#!/usr/bin/env python

'''
Download and install dependencies

Authors: Carmen Saenz
'''

# configuration yaml file
from os import path, mkdir
import glob

config_path  =  'configuration yaml file' #args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file: ") #, config_path)
print(" *******************************")
print("  ")

# Settings
virshimeome_dir = config["virshimeome_dir"]
local_dir = config["local_dir"]
download_threads = config["download_threads"]
download_memory = config["download_memory"]

# Download virsorter2 database and profiles. 
data_dir = directory(path.join(virshimeome_dir,"data"))
vs2_db_dir = path.join(data_dir, "vs2_db")
check_v_db_dir = path.join(data_dir, "checkv_db")

#################
##  Functions  ##
#################
# Define a function to expand the list of Pfam HMM files
def vs_pfam_hmm_files():
    pfam_dir_db = path.join(vs2_db_dir,"hmm", "pfam")
    pfam_hmm_list  =  [f"Pfam-A-{hmm_type}.hmm" for hmm_type in ["Archaea", "Bacteria", "Eukaryota", "Mixed", "Viruses"]]
    return expand("{pfam_dir_db}" + os.sep + "{hmm_file}", 
            pfam_dir_db = pfam_dir_db,
            hmm_file = pfam_hmm_list
            )
# Define a function to expand the list of RBS files
def vs_rbs_files():
    rbs_dir_db = path.join(vs2_db_dir,"rbs")
    rbs_list = [f"rbs-catetory{note}.tsv" for note in ["-notes"]]
    return expand("{rbs_dir_db}" + os.sep + "{rbs_file}",
            rbs_dir_db = rbs_dir_db, 
            rbs_file = rbs_list
            )
# Define a function to expand the list of group files
def vs_group_files():
    group_dir_db = path.join(vs2_db_dir,"group")
    group_list = ["dsDNAphage","lavidaviridae","NCLDV","RNA", "ssDNA"]
    group_files = ["hallmark-gene.list", "model"]
    return expand("{group_dir_db}" + os.sep + "{group}"+ os.sep +"{file}",
            group_dir_db = group_dir_db, 
            group = group_list, 
            file = group_files
            )
# Define a function to expand the list of CheckV databases
def check_v_genome_files(): 
    genome_db_dir = path.join(check_v_db_dir,"genome_db")
    genome_file_types = ["error.tsv", "info.tsv", "reps.dmnd", "reps.faa", "reps.fna", "reps.log", "reps.tsv"]
    genome_db_files = expand("{genome_db_dir}" + os.sep + "checkv_{genome_file_type}",
            genome_db_dir = genome_db_dir,
            genome_file_type = genome_file_types
            )    
    return genome_db_files
def check_v_hmm_files():
    hmm_db_dir = path.join(check_v_db_dir, "hmm_db", "checkv_hmms")
    hmm_index_range = [*range(1,81)]
    hmm_db_files = expand("{hmm_db_dir}" + os.sep + "{hmm_index}.hmm", 
            hmm_db_dir = hmm_db_dir, 
            hmm_index = hmm_index_range
            )
    return hmm_db_files

#############
##  Rules  ##
#############
rule all:
    input:
        check_v_genome_files(),
        check_v_hmm_files(),
        # vs_pfam_hmm_files(),
        vs_rbs_files(),
        # vs_group_files(),
        # vs_viral_hmm = path.join(vs2_db_dir, "hmm", "viral", "combined.hmm"),
        # vs_hmm_congregator = path.join(vs2_db_dir, "hmm", "pfam", "Pfam-A-acc2desc.tsv")


# Construct virsorter2 database
rule vs_two_db:
    output:
        # vs_pfam_hmm_files(),
        vs_rbs_files(),
        vs_group_files(),
        # vs_viral_hmm = path.join(vs2_db_dir,"hmm", "viral", "combined.hmm"),
        # vs_hmm_congregator = path.join(vs2_db_dir, "hmm", "pfam", "Pfam-A-acc2desc.tsv")
    threads:
        download_threads
    params:
        data_dir = data_dir,
        vs2_db_dir = vs2_db_dir
    conda:
        "virshimeome" # vs2
    shell:
        """
        mkdir -p {params.data_dir}
        virsorter setup -d {params.vs2_db_dir} -j {threads} --scheduler greedy
        """
# Construct checkV database. 
rule checv_db:
    output:
        check_v_genome_files(),
        check_v_hmm_files(),
        check_v_db_dir = directory(check_v_db_dir)
    params:
        data_dir = data_dir,
    conda:
        "virshimeome" #checkv
    shell:
        """
        mkdir -p {params.data_dir}
        checkv download_database {params.data_dir}
        mv {params.data_dir}/checkv-db* {output.check_v_db_dir}
        export CHECKVDB = {output.check_v_db_dir}
        """

