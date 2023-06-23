

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
min_contig_length = int(config["min_contig_length"])
max_contig_length = int(config["max_contig_length"])
min_score_vir_recognition = config["min_score_vir_recognition"]

# dvf contig selection
min_score_dvf = float(config["min_score_dvf"])
max_pval_dvf = float(config["max_pval_dvf"])

# Checkv
checkv_db_dir = path.join(virshimeome_dir, "data", "checkv_db")

# Read to contig alignment (BWA)
algorithm_bwa=config["algorithm_bwa"]

# Relative abundaance (samtools & msamtools)
min_length_alignment_abundance=config["min_length_alignment_abundance"]
ani_threshold=config["ani_threshold"]
percentage_read_aligned_abundance=config["percentage_read_aligned_abundance"]
# circular  yet to be determined if used. 

# Output #

# Global output directory, sample specific output directory. 
output_dir = config["output_dir"]
subdirs = ["1_1_vs", "1_1_dvf", "0_filtered_sequences"]
all_dirs = subdirs + ["2_checkv", "4_alignment", "3_checkm", "5_1_contig_lengths","5_2_viral_otus", "6_0_gene_calls", "6_1_all_against_all_search", "7_0_clustering", "8_0_identification", "9_0_distance", "data_visualization"]
make_dirs(path.join(output_dir, step_dir) for step_dir in all_dirs)
make_dirs([path.join(output_dir, "8_0_identification", sub_dir) for sub_dir in subdirs])
make_dirs([path.join(output_dir, "9_0_distance", sub_dir) for sub_dir in subdirs])
graph_filenames = ["circular_quality_absolute_bar.png", "contig_length_frequency.png", "circular_quality_percentage_bar.png", "contig_quality_boxplot.png"]

#############
##  Rules  ##
#############
rule all:
    input:
        # Merge and filter fasta files across MAGs
        expand("{main_dir}/0_filtered_sequences/combined_sequences.fa", 
            main_dir=output_dir
            ),

        # Split fasta
        expand("{main_dir}/0_filtered_sequences/seq_ids.tsv",
            main_dir=output_dir
            ),
        
        # fastANI
        expand("{main_dir}/9_0_distance/raw_all_against_all.out",
            main_dir=output_dir
            ),
        
        # Filtering sequences based on similarity
        expand("{main_dir}/0_filtered_sequences/revised_seq_ids.tsv",
            main_dir = output_dir
            ),
        
        # Apply duplicate filter
        expand("{main_dir}/0_filtered_sequences/final-viral-combined.fa", 
            main_dir=output_dir
            ),
        
        # Viral contig prediction (VS2)
        expand("{main_dir}/1_1_vs/final-viral-combined.fa", 
            main_dir=output_dir
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

        # CheckM takes a long time. 
        expand("{main_dir}/3_checkm/{type}/quality_summary.tsv", 
            main_dir=output_dir,
            type=subdirs
            ),
        
        # vOTU clutsering(script)
        # expand("{main_dir}/5_1_contig_lengths/{type}/contigs.all.lengths",
        #     main_dir=output_dir,
        #     type="1_1_dvf"
        #     ),

        # expand("{main_dir}/5_2_viral_otus/{type}/vOTUs.fna",
        #     main_dir=output_dir,
        #     type="1_1_dvf"
        #     ),

        # gene calling and comparison (prodigal)
        expand("{main_dir}/6_0_gene_calls/{type}/vOTUs.faa",
            main_dir = output_dir,
            type="1_1_dvf"
            ),

        # Sequence identification (phabox)
        expand("{main_dir}/8_0_identification/{type}/blast_results.tab",
            main_dir=output_dir,
            type="1_1_dvf"
            ),

        # Distance matrix (alfpy, euclid)
        expand("{main_dir}/9_0_distance/{type}/distance_euclid.mat",
            main_dir=output_dir,
            type="1_1_dvf"
            ),

        # Visualization
        # expand("{main_dir}/data_visualization/{files}.tsv", 
        #     main_dir=output_dir,
        #     files = graph_filenames
        #     )
        
###########################
##  fasta concatenation  ##
###########################
rule combine_sequence_files: # TODO This can be done in two in bash, but it works 
    """
    Runs a python instance that selects fasta files with a minimum amount of bases present in each file 
    assuming ech file contains one contig. 
    """
    input:
        sequence_files = lambda wildcards: glob.glob(
            os.path.join(main_contig_dir, "**", "*.fna")
            )
    output:
        combined_contig_file = "{main_dir}/0_filtered_sequences/combined_sequences.fa"
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

rule fasta_split:
    input:
        complete_fasta_file = "{main_dir}/0_filtered_sequences/combined_sequences.fa"
    output:
        seq_id_file = "{main_dir}/0_filtered_sequences/seq_ids.tsv"
    params:
        output_dir = "{main_dir}/0_filtered_sequences/split"
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
        seq_id_file = "{main_dir}/0_filtered_sequences/seq_ids.tsv"
    output:
        fastani_out = "{main_dir}/9_0_distance/raw_all_against_all.out"
    params:
        fastani_tenp = "{main_dir}/9_0_distance/raw_all_against_all_temp.out", 
        contig_length = min_contig_length, # There has already been a length trim but this value must exist. 
        ani_threshold = ani_threshold
    threads:
        20 
    conda:
        path.join(virshimeome_dir, "envs", "fastani.yml")
    shell:
        """
        fastANI \
            --rl {input.seq_id_file} \
            --ql {input.seq_id_file} \
            -t {threads} \
            --fragLen {params.contig_length} \
            -o {params.fastani_tenp} 
        
        echo fastANI \
            --rl {input.seq_id_file} \
            --ql {input.seq_id_file} \
            -t {threads} \
            --fragLen {params.contig_length} \
            -o {params.fastani_tenp} 

        awk -F'\t' '!($1 == $2) && $3 >= {params.ani_threshold}' {params.fastani_tenp} > {output.fastani_out}

        """

rule ani_filter:
    input:
        seq_id_file = "{main_dir}/0_filtered_sequences/seq_ids.tsv",
        fastani_out = "{main_dir}/9_0_distance/raw_all_against_all.out"
    output:
        revised_seq_ids = "{main_dir}/0_filtered_sequences/revised_seq_ids.tsv"
    params:
        ani_threshold = ani_threshold,
        ani_filter_script = path.join(script_dir, "fastani_filter.py")
    conda:
        path.join(virshimeome_dir, "envs", "py36_env.yml")
    shell:
        """
        python3 {params.ani_filter_script} \
            --seq_id_file {input.seq_id_file} \
            --fastani_out {input.fastani_out} \
            --revised_seq_ids {output.revised_seq_ids}
        """


rule concat_new_viral_sequences:
    input:
        revised_seq_ids = "{main_dir}/0_filtered_sequences/revised_seq_ids.tsv"
    output:
        output_fasta = "{main_dir}/0_filtered_sequences/final-viral-combined.fa"
    shell:
        """
        xargs -a {input.revised_seq_ids} -I{{}} cat {{}} > {output.output_fasta}
        """


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
        viral_boundary = "{main_dir}/1_1_vs/final-viral-boundary.tsv"
    params:
        min_score = min_score_vir_recognition,
        vs_db_dir = vs_db_dir
    threads:
        20
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
        20
    conda:
        path.join(virshimeome_dir, "envs", "deepvirfinder.yml")
    shell:
        """
        outdir=$(dirname {output.dvf_summary})        
        mkdir -p $outdir
        python {params.dvf_script} \
            -i {input.sequence_files} \
            -c {threads} \
            -o $outdir
            
        """

rule convert_dvf_results_to_sequences:
    input:
        dvf_summary = "{main_dir}/1_1_dvf/final-viral-combined.fa_gt1bp_dvfpred.txt",
        fasta_file = "{main_dir}/0_filtered_sequences/final-viral-combined.fa"
    output:
        viral_combined_dvf = "{main_dir}/1_1_dvf/final-viral-combined.fa"
    params:
        dvf_script_to_sequnces_script = path.join(script_dir, "dvf_to_sequences.py"),
        min_score_dvf=min_score_dvf,
        max_pval_dvf=max_pval_dvf
    shell:
        """
        python {params.dvf_script_to_sequnces_script} \
            --fasta_file {input.fasta_file} \
            --dvf_summary_file {input.dvf_summary}
            --viral_combined_dvf {output.viral_combined_dvf} \
            --min_score_dvf {params.min_score_dvf} \
            --max_pval_dvf {params.max_pval_dvf} 
        
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
        20
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
        20
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
        joincol_script = path.join(script_dir, "joincol"),
        fasta_fix_script = path.join(script_dir, "fasta_fix_script.py")
    threads:
        1
    shell:
        """
        cat {input.contigs} \
            | perl {params.f2s_script} \
            | perl {params.joincol_script} <(cut -f2 {input.viral_otus_table}) \
            | awk '$NF == 1' \
            | cut -f1,2 > s2f \
            > {output.viral_otus_fna}_temp
        
        cp {output.viral_otus_fna}_temp {output.viral_otus_fna}

        python3 {params.fasta_fix_script} --input {output.viral_otus_fna}_temp --output {output.viral_otus_fna}
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
        viral_otus_fna = "{main_dir}/5_2_viral_otus/{type}/vOTUs.fna"
    output:
        all_against_all_search = "{main_dir}/6_1_all_against_all_search/{type}/vOTUs.fasta36"
    params:
        fasta_36_bin = path.join(virshimeome_dir,"FASTA", "bin", "fasta36")
    threads:
        12
    shell:
        """
        outdir=$(dirname {output.all_against_all_search})        
        mkdir -p $outdir

        {params.fasta_36_bin} \
            -m {threads}  \
            {input.viral_otus_fna} \
            {input.viral_otus_fna} \
            > {output.all_against_all_search} 
        """

rule phabox_identification:
    input:
        contig_file = "{main_dir}/2_checkv/{type}/viruses.fna",
        protein_translations = "{main_dir}/6_0_gene_calls/{type}/vOTUs.faa"
    output:
        out_file = "{main_dir}/8_0_identification/{type}/blast_results.tab"
    params: 
        current_type = "{type}",
        outdir = "{main_dir}/8_0_identification/{type}",
        phabox_script = path.join(virshimeome_dir, "PhaBOX", "main.py"),
        database_dir = path.join(virshimeome_dir, "data", "phabox_db", "database"),
        parameter_dir = path.join(virshimeome_dir, "data", "phabox_db", "parameters")
    conda:
        path.join(virshimeome_dir, "envs", "phabox.yml") # Conda environment has a non-existent conflict with blast > 2.10.1, which causes problems
    threads:
        20
    shell:
        """
        mkdir -p {params.outdir}

        echo  python {params.phabox_script} \
            --threads {threads} \
            --contigs {input.contig_file} \
            --rootpth {params.outdir} \
            --out output/ \
            --dbdir {params.database_dir} \
            --parampth {params.parameter_dir} \
            --protein {input.protein_translations} \
            --reject 0 

         {params.phabox_script} \
            --threads {threads} \
            --contigs {input.contig_file} \
            --rootpth {params.outdir} \
            --out output/ \
            --dbdir {params.database_dir} \
            --parampth {params.parameter_dir} \
            --protein {input.protein_translations} \
            --reject 0 

        
        """

rule euclid_distance_calculation:
    input:
        fasta_file = "{main_dir}/2_checkv/{type}/final-viral-combined.fa"
    output:
        distance_matrix  = "{main_dir}/9_0_distance/{type}/distance_euclid.mat",
        seq_id_file =  "{main_dir}/9_0_distance/{type}/seq_ids.txt"
    params:
        distance_script = path.join(script_dir, "distance_calculation.py")
    threads:
        20
    conda:
        path.join(virshimeome_dir, "envs", "alfpy.yml")
    shell:
        """

        outdir=$(dirname {output.distance_matrix})        
        mkdir -p $outdir

        python {params.distance_script} \
            -input {input.fasta_file} \
            -output {output.distance_matrix} 

        echo python {params.distance_script} \
            -input {input.fasta_file} \
            -output {output.distance_matrix} 

        sed -i '1d' {output.distance_matrix}
        """



rule tree_creation:
    input:
        fasta_file = "{main_dir}/2_checkv/{type}/final-viral-combined.fa"
    output:
        tree_file = "{main_dir}/9_0_distance/{type}_tree.newick"
    params:
        spam_script = path.join(".", virshimeome_dir, "multi-SpaM", "bin", "multi-SpaM")
    threads:
        10 
    shell:
        """
        {params.spam_script} \
        -t {threads} \
        -i {input.fasta_file} \
        -o {output.tree_file}
        """

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
        visualizaton_script = path.join(script_dir, "pipeline_visualization.R"),
        method_comparison_visualization = path.join(script_dir, "viral_recognition_method_analysis.py") ,
        virshimeome_output_dir = output_dir
    shell:
        """
        outdir=$(dirname {output.outputfiles})        
        mkdir -p $outdir
        python3 {params.method_comparison_visualization} \
            -i {params.virshimeome_output_dir} \
            -o $outdir
        
        Rscript {params.visualizaton_script} 
        """