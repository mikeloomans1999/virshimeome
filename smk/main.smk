

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
visualization_script = path.join(virshimeome_dir, "scripts", "virshimeome_pipeline_visualization.py")

# Params #
# Resources
memory=config["memory"]
available_threads = workflow.cores

# Contig selection (custom script & virsorter2)
vs_db_dir = path.join(virshimeome_dir, "data", "vs_db")
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
all_dirs = subdirs + ["2_checkv", "4_alignment", "3_checkm", "5_1_contig_lengths","5_2_viral_otus", "6_0_gene_calls", "6_1_all_against_all_blastp", "7_0_clustering", "data_visualization"]
make_dirs(path.join(output_dir, step_dir) for step_dir in all_dirs)
graph_filenames = ["circular_quality_absolute_bar.png", "contig_length_frequency.png", "circular_quality_percentage_bar.png", "contig_quality_boxplot.png"]

#############
##  Rules  ##
#############
rule all:
    input:
        # Merge and filter fasta files across MAGs
        expand("{main_dir}/0_filtered_sequences/final-viral-combined.fa", 
            main_dir=output_dir, 
            ),
        
        # Viral contig prediction (VS2)
        expand("{main_dir}/1_1_vs/final-viral-combined.fa", 
            main_dir=output_dir, 
            ),

        # Viral contig prediction (DVF)
        expand("{main_dir}/1_1_dvf/final-viral-combined.fa_gt1bp_dvfpred.txt",
            main_dir=output_dir
            ),

        # DVF to sequences
        expand("{main_dir}/1_1_dvf/final-viral-combined.fa",
            main_dir=output_dir
            ),
        
        # CheckV 
        expand("{main_dir}/2_checkv/{type}/quality_summary.tsv", 
            main_dir=output_dir,
            type=subdirs
            ),

        # CheckM
        # expand("{main_dir}/3_checkm/{type}/quality_summary.tsv", 
        #     main_dir=output_dir,
        #     type=subdirs
        #     ),

        # Alignment
        expand("{main_dir}/4_alignment/{type}/contigs_aligned.blat", 
            main_dir=output_dir,
            type="1_1_dvf"
            ),
        
        # vOTU clutsering
        expand("{main_dir}/5_1_contig_lengths/{type}/contigs.all.lengths",
            main_dir=output_dir,
            type="1_1_dvf"
            ),
        expand("{main_dir}/5_2_viral_otus/{type}/vOTUs.tsv",
            main_dir=output_dir,
            type="1_1_dvf"
            ),

        # gene calling and comparison
        expand("{main_dir}/6_0_gene_calls/{type}/vOTUs.gbk",
            main_dir = output_dir,
            type="1_1_dvf"
            ),

        expand("{main_dir}/6_1_all_against_all_blastp/{type}/vOTUs.fasta36",
            main_dir = output_dir,
            type="1_1_dvf"
            )
            
        # Visualization
        # expand("{main_dir}/data_visualization/{files}.tsv", 
        #     main_dir=output_dir,
        #     files = graph_filenames
        #     )
        
###########################
##  fasta concatenation  ##
###########################
rule combine_sequence_files:
    """
    Runs a python instance that selects fasta files with a minimum amount of bases present in each file 
    assuming ech file contains one contig. 
    """
    input:
        sequence_files = lambda wildcards: glob.glob(
            os.path.join(main_contig_dir, "**", "*.fna")
            )
    output:
        combined_contig_file = "{main_dir}/0_filtered_sequences/final-viral-combined.fa"
    threads:
        1 
    params:
        max_contig_length=max_contig_length,
        min_contig_length=min_contig_length,
        circular=circular,
        sample_contig_dir=main_contig_dir
    run:        
        def read_fasta_file(fasta_file):
            # https://www.biostars.org/p/710/#383479
            with open(fasta_file, 'rb') as fasta_file_read:
                faiter = (x[1] for x in groupby(fasta_file_read, lambda line: str(line, 'utf-8')[0] == ">"))
                for header in faiter:
                    seq_id = str(header.__next__(), 'utf-8').replace("\n", "") 
                    seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__()) + "\n"
                    yield seq_id, seq

        # Open concat file in append mode and start adding sequences per MAG.
        with open(output.combined_contig_file, "a") as combined_fasta_file:
            for seq_file in input.sequence_files:
                for seq_id, seq in read_fasta_file(seq_file):
                    seq_id_list = seq_id.split("_")
                    length = int(seq_id_list[seq_id_list.index("length") + 1])
                    if params.circular in seq_id and params.min_contig_length <= length <= params.max_contig_length:
                        description = seq_file.split(os.sep)[-2] # MAG or sampleID
                        combined_fasta_file.write(f"{seq_id}_MAG_{description}\n{seq}")

#######################
##  vir recognition  ##
#######################
rule vir_sorter:
    """
    Runs Virsorter2 to predict viral contigs.
    """
    input:
        sequence_files = "{main_dir}/0_filtered_sequences/final-viral-combined.fa"
    output:
        viral_combined = "{main_dir}/1_1_vs/final-viral-combined.fa",
        viral_score = "{main_dir}/1_1_vs/final-viral-score.tsv",
        viral_boundary = "{main_dir}/1_1_vs/final-viral-boundary.tsv",  
    params:
        min_score = min_score_vir_recognition,
        vs_db_dir = vs_db_dir
    threads:
        10
    conda:
        path.join(virshimeome_dir, "envs", "virshimeome_base.yml") # vs2
    shell:
        """
        outdir=$(dirname {output.viral_combined})        
        mkdir -p $outdir
        echo virsorter run \
            -w $outdir \
            -i {input.sequence_files} \
            -j {threads} \
            all \
            --scheduler greedy

        virsorter config --init-source --db-dir={params.vs_db_dir}
        virsorter run \
            -w $outdir \
            -i {input.sequence_files} \
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
        sequence_files = "{main_dir}/0_filtered_sequences/final-viral-combined.fa"
    output:
        dvf_summary = "{main_dir}/1_1_dvf/final-viral-combined.fa_gt1bp_dvfpred.txt"
    params:
        dvf_script = path_dvf_script
    threads:
        10
    conda:
        path.join(virshimeome_dir, "envs", "deepvirfinder.yml")
    shell:
        """
        outdir=$(dirname {output.dvf_summary})        
        mkdir -p $outdir
        python {params.dvf_script} \
            -i {input.sequence_files} \
            -o $outdir \
            -c {threads}

        """

rule convert_dvf_results_to_sequences:
    input:
        dvf_summary = "{main_dir}/1_1_dvf/final-viral-combined.fa_gt1bp_dvfpred.txt",
        sequence_files = "{main_dir}/0_filtered_sequences/final-viral-combined.fa"
    output:
        viral_combined_dvf = "{main_dir}/1_1_dvf/final-viral-combined.fa"
    params:
        min_score_dvf=min_score_dvf,
        max_pval_dvf=max_pval_dvf
    run:
        def read_fasta_file(fasta_file):
            # https://www.biostars.org/p/710/#383479
            with open(fasta_file, 'rb') as fasta_file_read:
                faiter = (x[1] for x in groupby(fasta_file_read, lambda line: str(line, 'utf-8')[0] == ">"))
                print(type(faiter))
                print(fasta_file_read)
                for header in faiter:
                    seq_id = str(header.__next__(), 'utf-8').replace("\n", "")
                    seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__()) + "\n"
                    yield seq_id, seq

        # Get IDs of viral contigs matching filters. 
        viral_contig_ids = []
        with open(input.dvf_summary, "r") as dvf_file:
            dvf_summary_read = csv.reader(dvf_file, delimiter="\t")
            next(dvf_summary_read) # Skip header
            for row in dvf_summary_read:
                seq_id = row[0]
                score = float(row[2])
                pvalue = float(row[3])
                if score >= params.min_score_dvf and pvalue <= params.max_pval_dvf:
                    viral_contig_ids.append(">" + seq_id)

        # Write viral contigs to file.                 
        with open(output.viral_combined_dvf, "a") as dvf_viral_contig_file:
            for seq_id, seq in read_fasta_file(input.sequence_files):
                # Iterate through the entire object so the gc can clean up after.                
                if seq_id in viral_contig_ids:
                    dvf_viral_contig_file.write(seq_id + "\n" + seq)
        

# ##############
# ##  checkv  ##
# ##############
# TODO check usage without contamination step. 
rule checkv:
    """
    Runs CheckV which constructs quality reports of the viral prediction step, that classified the contigs. 
    """
    input:
        viral_combined = "{main_dir}/{type}/final-viral-combined.fa",
    output:
        contig_quality_summary = "{main_dir}/2_checkv/{type}/quality_summary.tsv",
        complete_genomes = "{main_dir}/2_checkv/{type}/complete_genomes.tsv",
        genome_completeness = "{main_dir}/2_checkv/{type}/completeness.tsv",
        checkv_viruses = "{main_dir}/2_checkv/{type}/viruses.fna",
        checkv_proviruses = "{main_dir}/2_checkv/{type}/proviruses.fna",
        checkv_contamination = "{main_dir}/2_checkv/{type}/contamination.tsv"
    params:
        checkv_db = checkv_db
    threads:
        10
    conda:
        path.join(virshimeome_dir, "envs", "virshimeome_base.yml") # checkv
    shell:
        """
        outdir=$(dirname {output.contig_quality_summary})        
        mkdir -p $outdir
        echo checkv end_to_end \
            {input.viral_combined} \
            $outdir \
            -t {threads} \
            -d {params.checkv_db} 
             
        checkv end_to_end \
            {input.viral_combined} \
            $outdir \
            -t {threads} \
            -d {params.checkv_db} 

        """

###############
##  checkm2  ##
###############
rule checkm:
    input:
        viral_combined = "{main_dir}/{type}/final-viral-combined.fa"
    output:
        contig_quality_summary = "{main_dir}/3_checkm/{type}/quality_summary.tsv",
    params:
        checkm_db = checkm_db
    conda:
        path.join(virshimeome_dir, "envs", "checkm2.yml")
    threads: 
        10
    shell:
        """
        outdir=$(dirname {output.contig_quality_summary})        
        mkdir -p $outdir

        tmp=$(mktemp -d)
        rm -rf $(dirname {output})
        time ( \

        echo checkm2 predict \
            --quiet \
            --database_path {params.checkm_db} \
            -x fna \
            --remove_intermediates \
            --threads {threads} \
            --input $(cat {input}) \
            --tmpdir $tmp \
            -o $(dirname $outdir) 

        checkm2 predict \
            --quiet \
            --database_path {params.checkm_db} \
            -x fna \
            --remove_intermediates \
            --threads {threads} \
            --input $(cat {input}) \
            --tmpdir $tmp \
            -o $(dirname $outdir) 
        rm -rf $tmp
        """


#################
##  alignment  ##
#################
rule alignment:
    input:
        contigs = "{main_dir}/{type}/final-viral-combined.fa"
    output:
        aligned_contigs = "{main_dir}/4_alignment/{type}/contigs_aligned.blat"
    threads:
        1 
    conda:
        path.join(virshimeome_dir, "envs", "virshimeome_base.yml") # BLAT
    shell:
        """
        blat \
            {input.contigs} \
            {input.contigs} \
            {output.aligned_contigs} \
            -out=blast8
        """

#######################
##  vOTU clustering  ##
#######################
rule contig_lengths:
    input:
        contigs = "{main_dir}/{type}/final-viral-combined.fa",
        aligned_contigs = "{main_dir}/4_alignment/{type}/contigs_aligned.blat"
    output:
        contig_lengths = "{main_dir}/5_1_contig_lengths/{type}/contigs.all.lengths"
    params:
        f2s_script = path.join(script_dir, "f2s"),
        seqlengths_script =path.join(script_dir, "seqlengths"),
        joincol_script = path.join(script_dir, "joincol"),
        hashsums_script = path.join(script_dir, "hashsums")
    threads:
        1
    shell:
        """
        cat {input.contigs} \
            | perl {params.f2s_script} \
            | perl {params.seqlengths_script} \
            | perl {params.joincol_script} <(cat {input.aligned_contigs} \
            | awk '{{if ($1 == $2) print $1 "\t" $12}}' \
            | perl {params.hashsums_script} \
            | tail -n +2) > {output.contig_lengths}
        """

rule viral_otu_clustering:
    input:
        contigs = "{main_dir}/{type}/final-viral-combined.fa",
        aligned_contigs = "{main_dir}/4_alignment/{type}/contigs_aligned.blat",
        contig_lengths = "{main_dir}/5_1_contig_lengths/{type}/contigs.all.lengths"
    output:
        viral_otus_table = "{main_dir}/5_2_viral_otus/{type}/vOTUs.tsv",
    params:
        joincol_script = path.join(script_dir, "joincol"),
        hashsums_script = path.join(script_dir, "hashsums")
    threads:
        1
    shell:
        """
        cut -f1,2,12 {input.aligned_contigs} \
            | perl {params.hashsums_script} \
            | tail -n +2 \
            | perl {params.joincol_script} {input.contig_lengths} 2 \
            | sort -k4,4nr -k1,1 \
            | awk '{{if ($3/$NF >= .90) print $1 "\t" $2}}' \
            | perl -lane 'unless (exists($clusters{{$F[1]}})) {{$clusters{{$F[1]}} = $F[0]; print "$F[1]\t$F[0]"}}' \
            > {output.viral_otus_table}
        """

rule viral_otu_cluster_sequences:
    input:
        contigs = "{main_dir}/{type}/final-viral-combined.fa",
        viral_otus_table = "{main_dir}/5_2_viral_otus/{type}/vOTUs.tsv",
    output:
        viral_otus_fna = "{main_dir}/5_2_viral_otus/{type}/vOTUs.fna" # Check for FASTA formatting. 
    params:
        f2s_script = path.join(script_dir, "f2s"),
        joincol_script = path.join(script_dir, "joincol")
    threads:
        1
    shell:
        """
        cat {input.contigs} \
            | perl {params.f2s_script} \
            | perl {params.joincol_script} <(cut -f2 {input.viral_otus_table}) \
            | awk '$NF == 1' \
            | awk -F'\t' '{printf "%s\n%s\n", $1, $2}' > s2f \
            && mv s2f {output.viral_otus_fna}
        """

rule gene_calling_protein_comparison:
    input:
        viral_otus_fna = "{main_dir}/5_2_viral_otus/{type}/vOTUs.fna"
    output:
        gene_calls = "{main_dir}/6_0_gene_calls/{type}/vOTUs.gbk",
        protein_translations = "{main_dir}/6_0_gene_calls/{type}/vOTUs.faa"
    conda:
        path.join(virshimeome_dir, "envs", "gene_calling.yml") # prodigal
    shell:
        """
        cat  {input.viral_otus_fna} | prodigal -a {output.protein_translations} -p meta > {output.gene_calls}
        """


rule all_against_all_fasta36:
    input:
        protein_translations = "{main_dir}/6_0_gene_calls/{type}/vOTUs.faa"
    output:
        all_against_all_fasta36 = "{main_dir}/6_1_all_against_all_blastp/{type}/vOTUs.fasta36"
    params:
        fasta_36_bin = path.join(virshimeome_dir,"FASTA", "bin", "fasta36")
    threads:
        12
    shell:
        """
        {params.fasta_36_bin} \
            -m {threads} 
            {input.protein_translations} \
            {input.protein_translations} \
            > {output.all_against_all_fasta36} 
        """

rule markov_chain_clustering:
    input:
        protein_translations = "{main_dir}/6_0_gene_calls/{type}/vOTUs.faa",
        all_against_all_fasta36 = "{main_dir}/6_1_all_against_all_blastp/{type}/vOTUs.fasta36"
    output:
        v_otu_lengths =  "{main_dir}/7_0_clustering/{type}/vOTUs.faa.lengths",
        v_otu_viral_orthologous_groups = "{main_dir}/7_0_clustering/{type}/vOTUs.faa"
    params:
        f2s_script = path.join(script_dir, "f2s"),
        seqlengths_script =path.join(script_dir, "seqlengths"),
        joincol_script = path.join(script_dir, "joincol"),
        hashsums_script = path.join(script_dir, "hashsums"),
        mcl_script = path.join(script_dir, "mcl")
    shell:
        """
        cat {input.protein_translations} \
            | perl {params.f2s_script} \
            | perl {params.seqlengths_script} \
            > {output.v_otu_lengths}

        cat {input.all_against_all_fasta36} \
            | perl {params.joincol_script} {output.v_otu_lengths} \
            | perl {params.joincol_script} {output.v_otu_lengths} 2 \
            | awk '{{print $1 "\t" $2 "\t" $11 "\t" $13/$14 "\t" ($8-$7)/(2*$13)+($10-$9)/(2*$14) "\t" ($7+$8-$9-$10)/($13+$14)}}' \
            | awk '{{if ($3 <= 0.05) print}}' \
            | awk '{{if ($5 >= 0.4) print}}' \
            | awk '{{if (sqrt(($4-1)^2) - (sqrt(sqrt($5))-.8) + sqrt($6^2) <= 0.1) print $1 "\t" $2}}' \
            | mcl - -o - --abc | awk '{{j++; for (i = 1; i <= NF; i++) {{print $i "\t" j}}}}' \
            > {output.v_otu_viral_orthologous_groups}
        """



# Map using bwa
# get contigs as reference and use the high quality paired-end reads to map to those reference reads. 
#################################
##  mappping to viral contigs  ##
#################################
# rule bwa_mem_indexing:
#     """
#     Creates an index for each concatenated fasta files for each sample. 
#     """
#     input:
#         viral_combined = viral_combined
#     output:
#         viral_combined_index = expand("{viral_combined}.bwt.2bit.64", 
#                         viral_combined = viral_combined)
#     threads:
#         threads
#     shell:
#         """
#         mkdir -p {output.viral_combined_index_dir}
#         bwa-mem2 index \
#             {input.viral_combined} 
#         """

# rule align_reads_to_contigs:
#     """
#     runs bwa-mem2 to align high quality reads to viral contigs and output the alignment file using samtools in CRAM format.
#     """
#     input:
#         reads_1=expand("{hq_reads_dir}/{sample_id}/{sample_id}.1.fq.gz", hq_reads_dir=hq_reads_dir, sample_id=sample_ids),
#         reads_2=expand("{hq_reads_dir}/{sample_id}/{sample_id}.2.fq.gz", hq_reads_dir=hq_reads_dir, sample_id=sample_ids),
#         viral_combined = viral_combined, # Contigs
#         viral_combined_index = expand("{viral_combined}.bwt.2bit.64", 
#                         viral_combined = viral_combined)
#     output:
#         bam_file_alignment=bam_file_alignment,
#         cram_file_alignment=cram_file_alignment
#     params:
#         algorithm_bwa=algorithm_bwa, # Default is mem
#         memory=memory,
#         cores=cores
#     conda:
#         path.join(virshimeome_dir, "conda_env", "virshimeome_base") # bwa & samtools
#     shell:
#         """
#         echo ' \
#         bwa-mem2 mem \
#             -t {threads} \
#             -x {params.algorithm_bwa} \
#             {input.viral_combined} \
#             {input.reads_1} {input.reads_2} \
#         | \
#         samtools view \
#             -F 4 \
#             -b \
#         | \
#         samtools sort \
#             -m {params.memory} \
#             -@ {params.cores} \
#             -T {output.bam_file_alignment} \
#         | \
#         samtools view \
#             -C \
#             -T {input.viral_combined} \
#             > {output.cram_file_alignment}
#         '
#         """
# ##########################
# ##  relative abundance  ##
# ##########################
# # Use msamtools later on to determine the relative abundance of viral species. 
# rule relative_abundance_estimation:
#     """
#     Estimates relative abundance of virus fragments found in each sample using msamtools and the previously created alignment file. 
#     """
#     input:
#         cram_file_alignment=cram_file_alignment,
#     output:
#         relative_abundance_profile_file=relative_abundance_profile_file #.txt.gz
#     params:
#         project_name=project_name,
#         min_length_alignment_abundance=min_length_alignment_abundance,
#         percentage_identity_abundance=percentage_identity_abundance,
#         percentage_read_aligned_abundance=percentage_read_aligned_abundance
#     conda:
#         path.join(virshimeome_dir, "conda_env", "virshimeome_base") # msamtools
#     shell:
#         # Normalize for sequence length. 
#         # Calculate relative abundance. 
#         """
#         echo msamtools filter \
#             -b -u \
#             -l {params.min_length_alignment_abundance} \
#             -p {params.percentage_identity_abundance} \
#             -z {params.percentage_read_aligned_abundance} \
#             --besthit {input.cram_file_alignment} \
#         | \
#         msamtools profile \
#             --label={params.project_name} \
#             -o {output.relative_abundance_profile_file}
                     
#         """

rule data_visualization:
    input:
        checkv_quality_summary = expand("{main_dir}/2_checkv/{type}/quality_summary.tsv",
            main_dir = output_dir,
            type=subdirs
            ),
        
        checkv_completeness_gneomes = expand("{main_dir}/2_checkv/{type}/complete_genomes.tsv",
            main_dir = output_dir,
            type=subdirs
            )
    output:
        outputfiles = "{main_dir}/data_visualization/{files}.tsv"
    params:
        visualizaton_script = visualization_script,
        virshimeome_output_dir = output_dir
    shell:
        """
        outdir=$(dirname {output.outputfiles})        
        mkdir -p $outdir
        python3 {params.visualizaton_script} \
            -i {params.virshimeome_output_dir} \
            -o $outdir

        """

