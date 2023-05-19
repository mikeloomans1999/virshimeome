
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
sample_ids = [abs_dir_path.split(os.sep)[-1] for abs_dir_path in get_subdirectories(hq_reads_dir)]
checkv_db=path.join(virshimeome_dir, "data", "checkv_db")
main_contig_dir=config["contig_dir"]
script_dir=path.join(virshimeome_dir, "scripts")
fasta_contig_file_paths = glob.glob(path.join(main_contig_dir, "**", "*.fasta"), recursive=True) # have an additional check in place for look for the fasta files here. 


# Params #
# Resources
memory=config["memory"]
threads = config["threads"]
cores=config["cores"]

# Contig selection (custom script & virsorter2)
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
sample_output_dirs = directory(expand("{output_dir}"+os.sep+"{sample_id}", 
                    output_dir=output_dir, sample_id=sample_ids))

# preprocessing. 
contig_concatenation_file_path = expand("{sample_output_dirs}/vir_recognition_contigs_combined_max_length.fasta",
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

make_dirs(sample_output_dirs)

#############
##  Rules  ##
#############
rule all:
    input:
        expand("{sample_output_dir}/vir_recognition_contigs_combined_max_length.fasta",
                sample_output_dir=sample_output_dirs, max_contig_length=max_contig_length),
        output_dir,
        # viral_combined,
        expand("{sample_output_dir}/final-viral-combined.fa", sample_output_dir = sample_output_dirs),
        contig_quality_summary,

###########################
##  fasta concatenation  ##
###########################
rule combine_fasta_files:
    """
    Runs a python instance that selects fasta files with a minimum amount of bases present in each file 
    assuming ech file contains one contig. 
    """
    input:
        fasta_contig_file_paths = fasta_contig_file_paths,
    output:
        contig_concatenation_file_path = "{sample_output_dir}/vir_recognition_contigs_combined_max_length.fasta", 
    params:
        output_dir = output_dir,
        main_contig_dir = main_contig_dir,
        max_contig_length = max_contig_length,
        script = path.join(script_dir, 'fasta_concatenation.py')
    run:
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

        print({params.sample_contig_dir})
        small_fasta_files = get_fasta_with_length_filter({params.max_contig_length}, {params.sample_contig_dir})
        
        circular_fasta_file_paths = [filename for filename in small_fasta_files if "circular" in filename] # Get small circular genomes.
        # Write combined fasta file 
        write_concat_fasta_files(circular_fasta_file_paths, {output.contig_concatenation_file_path})


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

##############
##  checkv  ##
##############
# TODO check usage without contamination step. 
rule checkv:
    """
    Runs CheckV which constructs quality reports of the viral prediction step, that classified the contigs. 
    """
    input:
        viral_combined = viral_combined,
        checkv_db = checkv_db
    output:
        contig_quality_summary = contig_quality_summary,
    threads:
        threads
    conda:
        path.join(virshimeome_dir, "envs", "virshimeome_base.yml")
    shell:
        """
        outdir=$(dirname {output.contig_quality_summary})        
        mkdir -p $outdir
        echo checkv end_to_end \
            -d {input.checkv_db}
            {input.viral_combined} \
            $outdir \
            -t {threads} \

        """
