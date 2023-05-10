
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
virshimeome_dir=config["virshimeome_dir"]
contig_dir = config["contig_dir"]
all_fasta_file_paths = glob.glob(path.join(contig_dir,"**", "*.fasta"), recursive=True)
hq_reads_dir=config["hq_reads_dir"] # This will be subdivided into numerous dirs and a things will have to be setup to iterate through them :(.

# Params
min_contig_length = config["min_contig_length"]
min_score_vir_recognition = config["min_score_vir_recognition"]
threads = config["threads"]
algorithm_bwa=config["algorithm_bwa"]
memory=config["memory"]
cores=config["cores"]
# circular  yet to be determined if used. 

# Output dir
output_dir=config["output_dir"]
# preprocessing. 
combined_fasta_file_path=path.join(output_dir, f"vir_recognition_contigs_combined_max_length_{max_contig_length}.fasta")
# Viral prediction
viral_combined = path.join(output_dir, "final-viral-combined.fa"),
viral_score = path.join(output_dir, "final-viral-score.tsv"),
viral_boundary = path.join(output_dir, "final-viral-boundary.tsv")
contig_quality_summary= path.join(output_dir,"quality_summary.tsv")

# Mapping to reads to contigs"
bam_file_alignment=path.join(output_dir, "reads_mapped_to_contigs.bam")
cram_file_alignment=path.join(output_dir, "reads_mapped_to_contigs.cram")

#############
##  Rules  ##
#############
rule all:
    input:
        combined_fasta_file_path,
        output_dir,
        viral_combined,
        viral_score,
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

#######################
##  vir recognition  ##
#######################
rule vir_recognition:
    input:
        fasta_files = combined_fasta_file_path
    output:
        viral_combined = viral_combined,
        viral_score = viral_score,
        viral_boundary = viral_boundary 
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

##############
##  checkv  ##
##############
# TODO check usage without contamination step. 
rule checkv:
    input:
        viral_combined=viral_combined,
    output:
        contig_quality_summary = contig_quality_summary,
    threads:
        threads
    shell:
        """
        outdir=$(dirname {output.contig_quality_summary})        
        mkdir -p $outdir
        echo checkv end_to_end \
            {viral_combined} \
            $outdir \
            -t {threads} \

        """

"""
The -x option in BWA-MEM specifies the algorithm and parameters to use for finding maximum exact matches (MEMs) during the seeding stage of read alignment. Specifically, the -x option takes a string argument that specifies the MEM algorithm and its parameters.

The available -x options are:

    -x ont2d: for Oxford Nanopore reads, which have high error rates, this option enables a more sensitive seeding algorithm that is better at handling the high error rate data.
    -x intractg: for long-read inter-chromosomal alignment, this option enables a more sensitive seeding algorithm that is better at handling inter-chromosomal gaps.
    -x sr: for short reads, this option uses the original BWA algorithm for seeding and is recommended for reads shorter than 70bp.
    -x spliced: for RNA-seq data, this option allows for the detection of spliced alignments. This mode uses a different seeding algorithm and allows for gap extension during the alignment process.

The default option is -x mem, which is recommended for most scenarios and uses a sensitive MEM algorithm suitable for short and long reads.
"""

# Map using bwa
# get contigs as reference and use the high quality paired-end reads to map to those reference reads. 
#################################
##  mappping to viral contigs  ##
#################################
rule relative_abundance_estimation:
    input:
        reads=path.join(hq_reads_dir, "sample_?.fg.gz"), # TODO reads
        viral_combined = viral_combined, # contigs
    output:
        bam_file_alignment=bam_file_alignment,
        cram_file_alignment=cram_file_alignment
    params:
        algorithm_bwa=algorithm_bwa, # mem is  default
        memory=memory,
        cores=cores
    shell:
        """
        echo bwa mem \
            -t {threads} \
            -a \
            -x {params.algorithm_bwa} \
            {input.viral_combined}  < (cat {input.reads}) \
        | \
        samtools view \
            -F4 \
            -b \
        | \
        samtools sort \
            -m {params.memory} -@ {params.cores} -T {output.bam_file_alignment} \
        | \
        samtools view \
            -C -T \
            {input.viral_combined} > {output.cram_file_alignment}

        """
##########################
##  relative abundance  ##
##########################
# Use msamtools later on to determine the relative abundance of viral species. 
rule:
    input:
        cram_file_alignment=cram_file_alignment,
    output:
        relative_abundance_profile_file=relative_abundance_profile_file #.txt.gz
    params:
        project_name=project_name,
        length_sequence_abundance=length_sequence_abundance,
        percentage_identity_abundance=percentage_identity_abundance,
        percentage_read_aligned_abundance=percentage_read_aligned_abundance,
    shell:
        # Normalize for sequence length. 
        # Calculate relative abundance. 
        """
        echo msamtools filter \
            -b -u \
            -l {params.length_sequence_abundance} \
            -p {params.percentage_identity_abundance} \
            -z {params.percentage_read_aligned_abundance} \
            --besthit {input.cram_file_alignment} \
        | \
        msamtools profile \
            --label={input.project_name} \
            -o {output.relative_abundance_profile_file}
                     
        """

# Run R code to make pretty graphs etc.. in another conda environment. 

# rule taxonomic_characterization:
