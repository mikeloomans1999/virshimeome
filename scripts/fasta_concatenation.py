# This is by far one of the worst python scripts I have written so far but it works and I will not deal with any more snakemake magic now. 

from os import path, sep
import os 
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, 
                    help='sample id ')
parser.add_argument('--max_contig_length', type=int, 
                    help='max_contig_length as an integer. ')
parser.add_argument('--output', type=str,
                    help='super directory of the sample output directories. ')
args = parser.parse_args()

main_contig_dir = args.input
max_contig_length = args.max_contig_length
output_dir = args.output

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
    
def get_fasta_with_length_filter(min_length, sample_contig_dir):
    """Returns all the fasta files in a specific directory that have less than the minimum length. """
    return [fasta_filename for fasta_filename in glob.glob(path.join(sample_contig_dir, "**", "*.fasta"), recursive=True) 
        if get_fasta_sequence_length_filename(fasta_filename) <= int(min_length)]

def write_concat_fasta_files(fasta_filepaths, combined_fasta_filepath):
    """ Writes combined fasta file from a given list of fasta file paths"""
    with open(combined_fasta_filepath, 'a') as combined_fasta_file:
        for fasta_filepath in fasta_filepaths:
            combined_fasta_file.write(f"{fasta_filepath.split(os.sep)[-1]}\n{get_fasta_sequence(fasta_filepath)}\n")

samples = [dir.replace(main_contig_dir,"").replace("/", "") for dir in glob.glob(path.join(main_contig_dir, "*/"))]
for sample in samples:
    sample_contig_dir = path.join(main_contig_dir, sample)
    concat_output_file_fasta  = path.join(output_dir, sample, f"vir_recognition_contigs_combined_max_length_{max_contig_length}.fasta")
    print(concat_output_file_fasta)
    small_fasta_files = get_fasta_with_length_filter(max_contig_length, sample_contig_dir)
    circular_fasta_file_paths = [filename for filename in small_fasta_files if "circular" in filename] # Get small circular genomes.
    # Write combined fasta file 
    write_concat_fasta_files(circular_fasta_file_paths, concat_output_file_fasta)
