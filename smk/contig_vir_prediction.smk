
#!/usr/bin/env python

'''
Selecting the right contigs and prediction viral sequences. 
Author: Mike Loomans
'''
include: 'virome_config_parser.smk'

###############
##  Imports  ##
###############
from os import path, getcwd
import glob

#################
##  Variables  ##
#################

# Input
max_contig_length = config["max_contig_length"]
working_dir = config["working_dir"]
local_dir = config["local_dir"]
contig_dir = config["contig_dir"]
all_fasta_file_paths = glob.glob(path.join(contig_dir,"**", "*.fasta"), recursive=True)
# pipeline_dir = config["pipeline_dir"] if config["pipeline_dir"] != "" else getcwd()

# Params
min_contig_length = config["min_contig_length"]
min_score_vir_recognition = config["min_score_vir_recognition"]
threads = config["threads"]
# circular  yet to be determined if used. 

# Output
output_dir=config["output_dir"]
combined_fasta_file_path=path.join(output_dir, f"vir_recognition_contigs_combined_max_length_{max_contig_length}.fasta")
viral_combined = path.join(output_dir, "final-viral-combined.fa"),
final_score = path.join(output_dir, "final-viral-score.tsv"),
viral_boundary = path.join(output_dir, "final-viral-boundary.tsv")


#############
##  Rules  ##
#############
rule all:
    input:
        combined_fasta_file_path,
        output_dir,
        viral_combined,
        final_score,
        viral_boundary

rule combine_fasta_files:
    input:
        fasta_file_path = expand("{file}",file=all_fasta_file_paths)
    output:
        combined_fasta_file_path
    params:
        max_contig_length = max_contig_length,
        contig_dir = contig_dir,
        combined_fasta_file_path = combined_fasta_file_path
    run:
        def get_fasta_sequence(fasta_filename):
            """ Reads the sequence of a fasta file containing only one sequence. 
            The filename specifies the absolute location of the fasta file. """
            with open(fasta_filename, 'r') as fasta_file:
                #next(fasta_file)
                # Concatenate all sequence lines into a single string
                return ''.join(line.strip() for line in fasta_file)

        def get_fasta_sequence_length_filename(fasta_filename):
            """ Return the length listed in the filename """
            filename_list = fasta_filename.split("_")
            return int(filename_list[filename_list.index("length") + 1].replace(".fasta",""))
        
        def get_fasta_files_under(min_length, contig_dir):
            """Returns all the fasta files in a specific directory that have less than the minimum length. """
            return [fasta_filename for fasta_filename in glob.glob(path.join(contig_dir,"**", "*.fasta"), recursive=True) 
                    if get_fasta_sequence_length_filename(fasta_filename) <= int(min_length)]

        def get_incomplete_genomes(contig_dir):
            """ Returns all the fasta files in a specific directory that do not have a  complete genome. """
            incomplete_genome_fasta_absolute_paths = []
            # Work in progress, unused as of now # 

            all_fasta_file_paths = glob.glob(path.join(contig_dir,"**", "*.fasta"), recursive=True)
            return incomplete_genome_absolute_fasta_file_paths

        def write_concat_fasta_files(fasta_filepaths, combined_fasta_filepath):
            """ Writes combined fasta file from a given list of fasta file paths"""
            with open(combined_fasta_filepath, 'a') as combined_fasta_file:
                for fasta_filepath in fasta_filepaths:
                    combined_fasta_file.write(f">{fasta_filepath.split(os.sep)[-1]}\n{get_fasta_sequence(fasta_filepath)}\n")

        small_fasta_files = get_fasta_files_under(max_contig_length, contig_dir)
        circular_fasta_file_paths = [filename for filename in small_fasta_files if "circular" in filename] # Get small circular genomes.
        # Write combined fasta file 
        write_concat_fasta_files(circular_fasta_file_paths, combined_fasta_file_path)

rule vir_recognition:
    input:
        fasta_files = combined_fasta_file_path
    output:
        viral_combined = path.join(output_dir, "final-viral-combined.fa"),
        viral_score = path.join(output_dir, "final-viral-score.tsv"),
        viral_boundary = path.join(output_dir, "final-viral-boundary.tsv")
    params:
        min_contig_length = min_contig_length,
        min_score = min_score_vir_recognition
    threads:
        threads
    shell:
        """
        outdir=$(dirname {output.viral_combined})        
        mkdir -p $outdir
        echo virsorter run \
            -w $outdir \
            -i {input.fasta_files} \
            -j {threads} \
            all
        
        virsorter run \
            -w $outdir \
            -i {input.fasta_files} \
            --min-length {params.min_contig_length} \
            --min-score {params.min_score} \
            -j {threads} \
            all

        """

# Check using checkV


# Map using bwa



# Run relative abundance estimate in another snakemake file. 
# msamtools profile --multi=proportional --label=sample1 --unit=rel -o sample1.IGC.profile.txt.gz sample1.IGC.filtered.bam


# rule taxonomic_characterization:
