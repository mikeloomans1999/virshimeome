
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
# from re import sub
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
# working_dir = config["working_dir"]
# local_dir = config["local_dir"]
virshimeome_dir=config["virshimeome_dir"]
hq_reads_dir=directory(config["hq_reads_dir"]) # This will be subdivided into numerous dirs and a things will have to be setup to iterate through them :(.
main_contig_dir=config["contig_dir"]
checkv_db=path.join(virshimeome_dir, "data", "checkv_db")
script_dir=path.join(virshimeome_dir, "scripts")
sample_ids = [abs_dir_path.split(os.sep)[-1] for abs_dir_path in get_subdirectories(main_contig_dir)]
fasta_contig_file_paths = glob.glob(path.join(main_contig_dir, "**", "*.fasta"), recursive=True) # have an additional check in place for look for the fasta files here. 
sample_ids = [prospective_dir.split("/")[-1] for prospective_dir in glob.glob(os.path.join(main_contig_dir, '*')) if os.path.isdir(prospective_dir)]

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

# Individual output directories, sample specific. 
sample_output_dirs = [os.path.join(output_dir,sample_id) for sample_id in sample_ids]
make_dirs(sample_output_dirs)

# preprocessing. 
contig_concatenation_file_path = expand("{sample_output_dirs}/contigs_combined.fasta",
                    sample_output_dirs=sample_output_dirs,max_contig_length=max_contig_length)

# Viral prediction output
viral_combined = expand("{sample_output_dir}/final-viral-combined.fa",
                    sample_output_dir = sample_output_dirs)
viral_score = expand("{sample_output_dir}/final-viral-score.tsv", 
                    sample_output_dir = sample_output_dirs)
viral_boundary = expand("{sample_output_dir}/final-viral-boundary.tsv", 
                    sample_output_dir = sample_output_dirs)
contig_quality_summary= expand("{sample_output_dir}/quality_summary.tsv", 
                    sample_output_dir=sample_output_dirs)

#############
##  Rules  ##
#############
rule all:
    input:
        expand("{main_dir}/{sample_dir}/combined_contigs.fasta", 
            main_dir=output_dir, 
            sample_dir=sample_ids, 
            ),
        
        # viral_combined,
        expand("{main_dir}/{sample_dir}/combined_contigs.fasta", 
            main_dir=output_dir, 
            sample_dir=sample_ids, 
            ),
        
        expand("{main_dir}/{sample_dir}/combined_contigs.fasta", 
            main_dir=output_dir, 
            sample_dir=sample_ids, 
            ),

        expand("{main_dir}/{sample_dir}/combined_contigs.fasta", 
            main_dir=output_dir, 
            sample_dir=sample_ids, 
            )
        # expand("{sample_output_dir}/final-viral-combined.fa", 
        #         sample_output_dir = sample_output_dirs),
        # contig_quality_summary,

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
            os.path.join(main_contig_dir, wildcards.sample_dir, "*.fna")
            )
    output:
        combined_contig_file = "{output_dir}/{sample_dir}/combined_contigs.fasta"
    params:
        max_contig_length=max_contig_length,
        min_contig_length=min_contig_length,
        circular=circular,
        sample_contig_dir=main_contig_dir
    run:
        print("Start combination script\n")
        
        def get_filtered_fasta(fasta_file, max_length, min_length, circular=""):
            """ returns sequence ID and sequence in list format if its header matches the correct values, 
                so every entry contains an ID and a sequence seperated by '\n' 
                which can be printed with a simple loop. """
            with open (fasta_file, 'rb') as fasta_file:
                # https://www.biostars.org/p/710/#383479
                faiter = (x[1] for x in groupby(fasta_file, lambda line: str(line, 'utf-8')[0] == ">"))
                for header in faiter:
                    seq_id = str(header.__next__(), 'utf-8')
                    seq_id_list = seq_id.split("_")
                    length = int(seq_id_list[seq_id_list.index("length")+1]) # TODO Move this step to an if function so it isn't excecuted every time.      
                    seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())+ "\n"
                    if not (circular in seq_id and length <= max_length and length >= min_length):
                        seq_id = "" # Change to null or None
                        seq = ""
                    yield (seq_id + seq)

        # Write new combined file.
        combined_fasta_filepath = output.combined_contig_file
        open(combined_fasta_filepath, mode='w').close()

        # get fasta files.
        fasta_files = [stringed_fasta_files[0] for stringed_fasta_files in {input.fasta_files}]

        # Open concat file in append mode and start adding sequences per MAG.
        with open(combined_fasta_filepath, "a") as combined_fasta_file:
            for seq_file in fasta_files:
                fasta_to_write = get_filtered_fasta(
                        fasta_file=seq_file, 
                        max_length=params.max_contig_length, 
                        min_length=params.min_contig_length, 
                        circular=params.circular
                        )
                for line in fasta_to_write:
                    combined_fasta_file.write(line) # Check if there is no double line stuff. 

#######################
##  vir recognition  ##
#######################
rule vir_recognition:
    """
    Runs Virsort2 to predict viral contigs for each sample. 
    """
    input:
        #fasta_files = contig_concatenation_file_path
        fasta_files = "{sample_output_dir}/vir_recognition_contigs_combined.fasta"
    output:
        # viral_combined = viral_combined,
        viral_combined = "{sample_output_dir}/final-viral-combined.fa",
        viral_score = "{sample_output_dir}/final-viral-score.tsv",
        # expand("{sample_output_dir}/final-viral-score.tsv", sample_output_dir = sample_output_dirs)
        viral_boundary = "{sample_output_dir}/final-viral-boundary.tsv",  
        # expand("{sample_output_dir}/final-viral-boundary.tsv", sample_output_dir = sample_output_dirs)
    params:
        min_contig_length = min_contig_length,
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
            --min-length {params.min_contig_length} \
            --min-score {params.min_score} \
            -j {threads} \
            all \
            --scheduler greedy
        """

# ##############
# ##  checkv  ##
# ##############
# # TODO check usage without contamination step. 
# rule checkv:
#     """
#     Runs CheckV which constructs quality reports of the viral prediction step, that classified the contigs. 
#     """
#     input:
#         viral_combined = viral_combined,
#         checkv_db = checkv_db
#     output:
#         contig_quality_summary = contig_quality_summary,
#     threads:
#         threads
#     conda:
#         path.join(virshimeome_dir, "envs", "virshimeome_base.yml")
#     shell:
#         """
#         outdir=$(dirname {output.contig_quality_summary})        
#         mkdir -p $outdir
#         echo checkv end_to_end \
#             -d {input.checkv_db}
#             {input.viral_combined} \
#             $outdir \
#             -t {threads} \

#         """
