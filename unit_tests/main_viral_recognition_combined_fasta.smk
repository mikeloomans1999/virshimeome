
#!/usr/bin/env python

'''
Selecting the right contigs and prediction viral sequences. 
Author: Mike Loomans
'''
include: 'main_config_parser.smk'

###############
##  Imports  ##
###############
from os import path, mkdir
import glob
from itertools import groupby

#################
##  Functions  ##
#################

def get_subdirectories(main_dir):
    # Get first layer subdirectories. 
    return [sample_dir for sample_dir in glob.glob(path.join(main_dir,"*")) if os.path.isdir(sample_dir)]

def make_dirs(abs_dir_list):
    [os.mkdir(directory) for directory in abs_dir_list if not os.path.isdir(directory)]

#################
##  Variables  ##
#################

# Input # 
project_name=config["PROJECT"]
virshimeome_dir=config["virshimeome_dir"]
hq_reads_dir=directory(config["hq_reads_dir"]) # This will be subdivided into numerous dirs and a things will have to be setup to iterate through them :(.
main_contig_dir=config["contig_dir"]
checkv_db=path.join(virshimeome_dir, "data", "checkv_db")
script_dir=path.join(virshimeome_dir, "scripts")
fasta_contig_file_paths = glob.glob(path.join(main_contig_dir, "**", "*.fasta"), recursive=True) # have an additional check in place for look for the fasta files here. 
sample_ids = [prospective_dir.split("/")[-1] for prospective_dir in glob.glob(os.path.join(main_contig_dir, '*')) if os.path.isdir(prospective_dir)]
path_dvf_script = path.join(config["deep_vir_finder_dir"], "dvf.py")

# Params #
# Resources
memory=config["memory"]
threads = config["threads"]
cores=config["cores"]

# Contig selection (custom script & virsorter2)
circular = config["circular"]
min_contig_length = config["min_contig_length"]
max_contig_length = config["max_contig_length"]
min_score_vir_recognition = config["min_score_vir_recognition"]

# Read to contig alignment (BWA)
algorithm_bwa=config["algorithm_bwa"]

# Relative abundaance (samtools & msamtools)
min_length_alignment_abundance=config["min_length_alignment_abundance"]
percentage_identity_abundance=config["percentage_identity_abundance"]
percentage_read_aligned_abundance=config["percentage_read_aligned_abundance"]
# circular  yet to be determined if used. 

# Output #

# Global output directory, sample specific output directory. 
output_dir = config["output_dir"]
make_dirs(path.join(output_dir, step_dir) for step_dir in ["1_1_vs","1_2_checkv"])

#############
##  Rules  ##
#############
rule all:
    input:
        expand("{main_dir}/combined_contigs.fasta", 
            main_dir=output_dir, 
            ),
        
        # Viral contig prediction (VS2)
        expand("{main_dir}/1_1_vs/final-viral-combined.fa", 
            main_dir=output_dir, 
            ),

        # Viral contig prediction (DVF)
        expand("{main_dir}/1_1_dvf/combined_contigs.fasta_gt1bp_dvfpred.txt",
            main_dir=output_dir
            ),
        
        # CheckV 
        expand("{main_dir}/1_2_checkv/quality_summary.tsv", 
            main_dir=output_dir,
            )

###########################
##  fasta concatenation  ##
###########################
rule combine_fasta_files:
    """
    Runs a python instance that selects fasta files with a minimum amount of bases present in each file 
    assuming ech file contains one contig. 
    """
    input:
        fasta_files = lambda wildcards: glob.glob(
            os.path.join(main_contig_dir, "**", "*.fna")
            )
    output:
        combined_contig_file = "{output_dir}/combined_contigs.fasta"
    threads:
        1 
    params:
        max_contig_length=max_contig_length,
        min_contig_length=min_contig_length,
        circular=circular,
        sample_contig_dir=main_contig_dir
    run:        
        def get_filtered_fasta(fasta_file, max_length, min_length, circular="", description=""):
            # https://www.biostars.org/p/710/#383479
            with open(fasta_file, 'rb') as fasta_file_read:
                faiter = (x[1] for x in groupby(fasta_file_read, lambda line: str(line, 'utf-8')[0] == ">"))
                for header in faiter:
                    seq_id = str(header.__next__(), 'utf-8').replace("\n", "") + "_MAG_" + description + "\n"
                    seq_id_list = seq_id.split("_")
                    length = int(seq_id_list[seq_id_list.index("length") + 1])
                    seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__()) + "\n"
                    if not (circular in seq_id and min_length <= length <= max_length):
                        seq_id = seq = ""
                    yield seq_id + seq

        # Open concat file in append mode and start adding sequences per MAG.
        with open(output.combined_contig_file, "a") as combined_fasta_file:
            for seq_file in input.fasta_files:
                fasta_to_write = get_filtered_fasta(
                    fasta_file=seq_file,
                    max_length=params.max_contig_length,
                    min_length=params.min_contig_length,
                    circular=params.circular,
                    description=seq_file.split(os.sep)[-2]
                )
                combined_fasta_file.writelines(fasta_to_write)


#######################
##  vir recognition  ##
#######################
rule vir_sorter:
    """
    Runs Virsorter2 to predict viral contigs.
    """
    input:
        fasta_files = "{output_dir}/combined_contigs.fasta"
    output:
        viral_combined = "{output_dir}/1_1_vs/final-viral-combined.fa",
        viral_score = "{output_dir}/1_1_vs/final-viral-score.tsv",
        viral_boundary = "{output_dir}/1_1_vs/final-viral-boundary.tsv",  
    params:
        min_score = min_score_vir_recognition
    threads:
        threads
    conda:
        path.join(virshimeome_dir, "envs", "virshimeome_base.yml")
    shell:
        """
        outdir=$(dirname {output.viral_combined})        
        mkdir -p $outdir
        echo virsorter run \
            -w $outdir \
            -i {input.fasta_files} \
            -j {threads} \
            all \
            --scheduler greedy

        virsorter run \
            -w $outdir \
            -i {input.fasta_files} \
            --min-score {params.min_score} \
            -j {threads} \
            all \
            --scheduler greedy
        """

rule deep_vir_finder:
    """
    Runs DeepVirfinder to predict viral contigs.
    """
    input:
        fasta_files = "{output_dir}/combined_contigs.fasta"
    output:
        dvf_summary = "{output_dir}/1_1_dvf/combined_contigs.fasta_gt1bp_dvfpred.txt"
    params:
        dvf_script = path_dvf_script
    threads:
        threads
    conda:
        path.join(virshimeome_dir, "envs", "deepvirfinder.yml")
    shell:
        """
        outdir=$(dirname {output.dvf_summary})        
        mkdir -p $outdir
        python {params.dvf_script} \
            -i {input.fasta_files} \
            -o $outdir \
            -c {threads}

        """

# ##############
# ##  checkv  ##
# ##############
# TODO check usage without contamination step. 
rule checkv:
    """
    Runs CheckV which constructs quality reports of the viral prediction step, that classified the contigs. 
    """
    input:
        viral_combined = "{output_dir}/1_1_vs/final-viral-combined.fa",
        checkv_db = checkv_db
    output:
        contig_quality_summary = "{output_dir}/1_2_checkv/quality_summary.tsv",
    threads:
        threads
    conda:
        path.join(virshimeome_dir, "envs", "virshimeome_base.yml")
    shell:
        """
        outdir=$(dirname {output.contig_quality_summary})        
        mkdir -p $outdir
        echo checkv end_to_end \
            {input.viral_combined} \
            $outdir \
            -t {threads} \
            -d {input.checkv_db} 
             
        checkv end_to_end \
            {input.viral_combined} \
            $outdir \
            -t {threads} \
            -d {input.checkv_db} 

        """
