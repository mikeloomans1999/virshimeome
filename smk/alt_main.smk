

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
from csv import reader 
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

# Iut # 
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
            type="0_filtered_sequences"
            ),
        
        # fastANI
        expand("{main_dir}/{type}/raw_all_against_all.out",
            main_dir=output_dir,
            type="0_filtered_sequences"
            ),
        
        # Filtering sequences based on similarity
        expand("{main_dir}/{type}/revised_seq_ids.tsv",
            main_dir = output_dir,
            type="0_filtered_sequences"
            )


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
            with open(output_file, "a") as output_object:
                output_object.write(seq_id_path + "\n")

        def write_fasta_file(seq_id, sequence, output_file):
            sequences_formatted = "\n".join([sequence[x:x+70] for x in range(0, len(sequence), 70)])
            with open(output_file, 'w') as fasta_file:
                fasta_file.write(f">{seq_id}\n{sequences_formatted}\n")

        if not os.path.exists(params.output_dir):
            os.makedirs(params.output_dir)

        # function calling 
        split_fasta(input.complete_fasta_file, params.output_dir, output.seq_id_file)

rule fastani:
    input:
        seq_id_file = "{main_dir}/{type}/seq_ids.tsv"
    output:
        fastani_out = "{main_dir}/{type}/raw_all_against_all.out",
        fastani_out_matrix = "{main_dir}/{type}/raw_all_against_all.out.matrix"
    params:
        contig_lenght = 2000
    threads:
        10 
    conda:
        path.join(virshimeome_dir, "envs", "fastani.yml")
    shell:
        """
        fastANI \
            --rl {input.seq_id_file} \
            --ql {input.seq_id_file} \
            -t {threads} \
            --fragLen {params.contig_lenght} \
            --matrix \
            -o {output.fastani_out} \
            > fastANI.log

        """

rule ani_filter:
    input:
        seq_id_file = "{main_dir}/{type}/seq_ids.tsv",
        fastani_out = "{main_dir}/{type}/raw_all_against_all.out",
        fastani_out_matrix = "{main_dir}/{type}/raw_all_against_all.out.matrix"
    output:
        revised_seq_ids = "{main_dir}/{type}/revised_seq_ids.tsv"
    params:
        ani_threshold = 99.9
    run:

        def convert_set_to_file(set_input, output_file):
            with open(output_file, "a") as new_file:
                for seq_id in set_input:
                    new_file.write(f"{seq_id}\n")

        # Get all ids
        with open(input.seq_id_file, "r") as seq_id_no_filter:
            all_ids = set(seq_id_no_filter.read().splitlines())
            
        identical_ids = []

        # Obatain all 
        removed = 0
        with open(input.fastani_out, "r") as fastani_out:
            fastani_out_read = reader(fastani_out, delimiter="\t")
            for row in fastani_out_read:
                if float(row[2]) >= float(params.ani_threshold) and row[0] != row[1]:
                    # Check if the current identical set of sequences is already mentioned here.    
                    if row[0] in all_ids and row[1] in all_ids:
                        all_ids.remove(row[0])
                        all_ids.remove(row[1])
                        identical_ids.append({row[0], row[1]})
                        removed += 2
                    else:
                        # One of the ids is already present in the set and one has to be added. 
                        # If there has been a mistake of some kind this will just do nothing. 
                        index_of_interest = 0 if row[0] in all_ids else 1 # so if index 0 is not in all ids then the other must be. 
                        if row[0] in all_ids:
                            index_of_interest = 0
                        elif row[1] in all_ids: # Not in ids
                            index_of_interest = 1
                        else:
                            continue
                            # None in ids are presumably already correctly identified. 
                            
                        for id_set in identical_ids:
                            if row[index_of_interest] in id_set:
                                id_set.add(row[index_of_interest])
                                removed += 1
                                # Remove from all_ids so we're left with just singletons. 
                                all_ids.remove(row[index_of_interest])
                                break
        # Apply filter
        longest_identical_ids = set()
        for id_set in identical_ids:
            longest_contig = 0 
            for id in id_set:
                seq_id_list = id.split("_")
                contig_length = int(seq_id_list[seq_id_list.index("length") + 1])
                if contig_length > longest_contig:
                    longest_id = id
            
            longest_identical_ids.add(id)
        removed -= len(longest_identical_ids)
        revised_ids = all_ids.union(longest_identical_ids)
        print(removed)
        convert_set_to_file(revised_ids, output.revised_seq_ids)



# See a 100% match
# Check existing list with sets for either ID in each set.  I can do this in parallel. 
# If neither ID exists create a new set entry with both inside. 
# If one of the IDs is found in a set add both to that set. (the duplicate will not get added. )

