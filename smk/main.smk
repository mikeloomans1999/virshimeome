
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
contig_concatenation_file_path = expand("{sample_output_dirs}/vir_recognition_contigs_combined_max_length_{max_contig_length}.fasta",
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

# Mapping to reads to contigs
bam_file_alignment= expand("{sample_output_dir}/reads_mapped_to_contigs.bam", 
                    sample_output_dir = sample_output_dirs)
cram_file_alignment= expand("{sample_output_dir}/reads_mapped_to_contigs.cram", 
                    sample_output_dir = sample_output_dirs)

# Relative abundance
relative_abundance_profile_file= expand("{sample_output_dir}" + os.sep +  "relative_abundance_summary.txt.gz", 
                    sample_output_dir = sample_output_dirs)

make_dirs(sample_output_dirs)

#############
##  Rules  ##
#############
rule all:
    input:
        #contig_concatenation_file_path,
        expand("{sample_output_dirs}/vir_recognition_contigs_combined_max_length_{max_contig_length}.fasta",
                sample_output_dirs=sample_output_dirs,max_contig_length=max_contig_length),
        output_dir,
        expand("{sample_output_dir}/final-viral-combined.fa", sample_output_dir = sample_output_dirs),
        viral_score,
        viral_boundary,
        contig_quality_summary,
        cram_file_alignment,
        relative_abundance_profile_file

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
        contig_concatenation_file_path = contig_concatenation_file_path, 
    params:
        output_dir = output_dir,
        main_contig_dir = main_contig_dir,
        max_contig_length = max_contig_length,
        script = path.join(script_dir, 'fasta_concatenation.py')
    shell:    
        """
        python3 {params.script} \
            --input {params.main_contig_dir} \
            --max_contig_length {params.max_contig_length} \
            --output {params.output_dir}
        """

#######################
##  vir recognition  ##
#######################
rule vir_recognition:
    """
    Runs Virsort2 to predict viral contigs for each sample. 
    """
    input:
        #fasta_files = contig_concatenation_file_path
        expand("{sample_output_dirs}/vir_recognition_contigs_combined_max_length_{max_contig_length}.fasta", 
                sample_output_dirs=sample_output_dirs, 
                max_contig_length=max_contig_length)
    output:
        # viral_combined = viral_combined,
        viral_combined = expand("{sample_output_dir}/final-viral-combined.fa", sample_output_dir = sample_output_dirs),
        viral_score = viral_score,
        # expand("{sample_output_dir}/final-viral-score.tsv", sample_output_dir = sample_output_dirs)
        viral_boundary = viral_boundary 
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

# Map using bwa
# get contigs as reference and use the high quality paired-end reads to map to those reference reads. 
#################################
##  mappping to viral contigs  ##
#################################
rule bwa_mem_indexing:
    """
    Creates an index for each concatenated fasta files for each sample. 
    """
    input:
        viral_combined = viral_combined
    output:
        viral_combined_index = expand("{viral_combined}.bwt.2bit.64", 
                        viral_combined = viral_combined)
    threads:
        threads
    shell:
        """
        mkdir -p {output.viral_combined_index_dir}
        bwa-mem2 index \
            {input.viral_combined} 
        """

rule align_reads_to_contigs:
    """
    runs bwa-mem2 to align high quality reads to viral contigs and output the alignment file using samtools in CRAM format.
    """
    input:
        reads_1=expand("{hq_reads_dir}/{sample_id}/{sample_id}.1.fq.gz", hq_reads_dir=hq_reads_dir, sample_id=sample_ids),
        reads_2=expand("{hq_reads_dir}/{sample_id}/{sample_id}.2.fq.gz", hq_reads_dir=hq_reads_dir, sample_id=sample_ids),
        viral_combined = viral_combined, # Contigs
        viral_combined_index = expand("{viral_combined}.bwt.2bit.64", 
                        viral_combined = viral_combined)
    output:
        bam_file_alignment=bam_file_alignment,
        cram_file_alignment=cram_file_alignment
    params:
        algorithm_bwa=algorithm_bwa, # Default is mem
        memory=memory,
        cores=cores
    conda:
        path.join(virshimeome_dir, "conda_env", "virshimeome_base") # bwa & samtools
    shell:
        """
        echo ' \
        bwa-mem2 mem \
            -t {threads} \
            -x {params.algorithm_bwa} \
            {input.viral_combined} \
            {input.reads_1} {input.reads_2} \
        | \
        samtools view \
            -F 4 \
            -b \
        | \
        samtools sort \
            -m {params.memory} \
            -@ {params.cores} \
            -T {output.bam_file_alignment} \
        | \
        samtools view \
            -C \
            -T {input.viral_combined} \
            > {output.cram_file_alignment}
        '
        """
##########################
##  relative abundance  ##
##########################
# Use msamtools later on to determine the relative abundance of viral species. 
rule relative_abundance_estimation:
    """
    Estimates relative abundance of virus fragments found in each sample using msamtools and the previously created alignment file. 
    """
    input:
        cram_file_alignment=cram_file_alignment,
    output:
        relative_abundance_profile_file=relative_abundance_profile_file #.txt.gz
    params:
        project_name=project_name,
        min_length_alignment_abundance=min_length_alignment_abundance,
        percentage_identity_abundance=percentage_identity_abundance,
        percentage_read_aligned_abundance=percentage_read_aligned_abundance
    conda:
        path.join(virshimeome_dir, "conda_env", "virshimeome_base") # msamtools
    shell:
        # Normalize for sequence length. 
        # Calculate relative abundance. 
        """
        echo msamtools filter \
            -b -u \
            -l {params.min_length_alignment_abundance} \
            -p {params.percentage_identity_abundance} \
            -z {params.percentage_read_aligned_abundance} \
            --besthit {input.cram_file_alignment} \
        | \
        msamtools profile \
            --label={params.project_name} \
            -o {output.relative_abundance_profile_file}
                     
        """

# Run R code to make pretty graphs etc.. in another conda environment. 

# rule taxonomic_characterization:
