

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
import glob, csv
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
checkm_db=path.join(virshimeome_dir, "data", "CheckM2_Database")
script_dir=path.join(virshimeome_dir, "scripts")
fasta_contig_file_paths = glob.glob(path.join(main_contig_dir, "**", "*.fasta"), recursive=True) # have an additional check in place for look for the fasta files here. 
sample_ids = [prospective_dir.split("/")[-1] for prospective_dir in glob.glob(os.path.join(main_contig_dir, '*')) if os.path.isdir(prospective_dir)]
path_dvf_script = path.join(config["deep_vir_finder_dir"], "dvf.py")

# Params #
# Resources
memory=config["memory"]
available_threads = workflow.cores

# Contig selection (custom script & virsorter2)
vs_db_dir = path.join(virshimeome_dir, "data", "vs2_db")
circular = config["circular"]
min_contig_length = config["min_contig_length"]
max_contig_length = config["max_contig_length"]
min_score_vir_recognition = config["min_score_vir_recognition"]

# dvf contig selection
min_score_dvf = config["min_score_dvf"]
max_pval_dvf = config["max_pval_dvf"]

# Checkv
checkv_db_dir = path.join(virshimeome_dir, "data", "checkv_db")

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
subdirs = ["1_1_vs", "1_1_dvf", "0_filtered_sequences"]
all_dirs = subdirs + ["2_checkv", "4_alignment", "3_checkm", "5_1_contig_lengths","5_2_viral_otus", "6_0_gene_calls", "6_1_all_against_all_search", "7_0_clustering", "8_0_distances", "data_visualization"]
make_dirs(path.join(output_dir, step_dir) for step_dir in all_dirs)
graph_filenames = ["circular_quality_absolute_bar.png", "contig_length_frequency.png", "circular_quality_percentage_bar.png", "contig_quality_boxplot.png"]

#############
##  Rules  ##
#############
rule all:
    input:
        # Split fasta
        expand("{main_dir}/{type}/seq_ids.tsv",
            main_dir=output_dir,
            type="1_1_dvf"
            ),
        
        # Sketch reference. 
        expand("{main_dir}/{type}/final-viral-combined.fa.msh",
            main_dir=output_dir,
            type="1_1_dvf"
            ),

        expand("{main_dir}/8_0_distances/{type}_vOTUs.mat",
            main_dir=output_dir,
            type="1_1_dvf"
            ),

        expand("{main_dir}/8_0_distances/{type}_vOTUs_ANDI.mat",
            main_dir=output_dir,
            type="1_1_dvf"
            )

        # Visualization
        # expand("{main_dir}/data_visualization/{files}.tsv", 
        #     main_dir=output_dir,
        #     files = graph_filenames
        #     )


rule fasta_split:
    input:
        complete_fasta_file = "{main_dir}/{type}/final-viral-combined.fa"
    output:
        seq_id_file = "{main_dir}/{type}/seq_ids.tsv"
    params:
        output_dir = "{main_dir}/{type}/split"
    run:
        def split_fasta(input_file, fasta_output_dir, seq_id_file):
            with open(input_file, 'r') as fasta_file:
                seq_id = ''
                sequence = ''
                for line in fasta_file:
                    line = line.strip()
                    if line.startswith('>'):
                        if seq_id and sequence:
                            fasta_output_file = f"{fasta_output_dir}/{seq_id}.fa"
                            write_fasta_file(seq_id, sequence, fasta_output_file)
                            write_seq_path(fasta_output_file, seq_id_file)
                        seq_id = line[1:]
                        sequence = ''
                    else:
                        sequence += line
                if seq_id and sequence:
                    write_fasta_file(seq_id, sequence, fasta_output_file)

        def write_seq_path(seq_id_path, output_file):
            with open(output_file, "a") as seq_id_file:
                seq_id_file.write(seq_id_path + "\n")

        def write_fasta_file(seq_id, sequence, output_file):
            sequences_formatted = "\n".join([sequence[x:x+70] for x in range(0, len(sequence), 70)])
            with open(output_file, 'w') as fasta_file:
                fasta_file.write(f">{seq_id}\n{sequences_formatted}\n")

        if not os.path.exists(params.output_dir):
            os.makedirs(params.output_dir)

        # function calling 
        split_fasta(input.complete_fasta_file, params.output_dir, output.seq_id_file)

rule mash_sketching:
    input:
        complete_seq_file = "{main_dir}/{type}/final-viral-combined.fa"
    output:
        v_contigs_sketches = "{main_dir}/{type}/final-viral-combined.fa.msh"
    params:
        mash_bin = path.join(".", virshimeome_dir, "mash", "mash")
    threads:
        10
    shell:
        """
            mash sketch -p {threads} -o {output.v_contigs_sketches} -i {input.complete_seq_file}
        """

rule aggregate_protein_similarity:
    input:
        v_contigs_sketches = "{main_dir}/{type}/final-viral-combined.fa.msh",
        seq_id_file = "{main_dir}/{type}/seq_ids.tsv"
    output:
        v_otu_distance_matrix = "{main_dir}/8_0_distances/{type}_vOTUs.mat"
    threads:
        10
    shell:
        """
           mash dist -p {threads} {input.v_contigs_sketches} -l {input.seq_id_file} > {output.v_otu_distance_matrix}
        """

rule distance_calc:
    input:
        v_contigs_sketches = "{main_dir}/{type}/final-viral-combined.fa.msh",
        seq_id_file = "{main_dir}/{type}/seq_ids.tsv"
    output:
        v_otu_distance_matrix = "{main_dir}/8_0_distances/{type}_vOTUs_ANDI.mat"
    threads:
        10
    shell:
        """
            outdir=$(dirname {input}.seq_id_file})        
            cd $outdir 
            cd split/
            cat * \
            | andi \
            > {output.v_otu_distance_matrix}

        """