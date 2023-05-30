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
env_dir = path.join(virshimeome_dir, "envs")
local_dir = config["local_dir"]
download_threads = config["download_threads"]
download_memory = config["download_memory"]

# Download virsorter2 database and profiles. 
data_dir = directory(path.join(virshimeome_dir, "data"))
vs2_db_dir = path.join(data_dir, "vs2_db")
checkv_db_dir = path.join(data_dir, "checkv_db")
checkm_db_dir = path.join(data_dir, "checkm_db")


#################
##  Functions  ##
#################
# Define a function to expand the list of Pfam HMM files

def vs_out()
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

    return vs_pfam_hmm_files + vs_rbs_files + vs_group_files


def checkv_out()
    # Define a function to expand the list of CheckV databases
    def check_v_genome_files(): 
        genome_db_dir = path.joincheck_v_db_dir,"genome_db")
        genome_file_types = ["error.tsv", "info.tsv", "reps.dmnd", "reps.faa", "reps.fna", "reps.log", "reps.tsv"]
        genome_db_files = expand("{genome_db_dir}" + os.sep + "checkv_{genome_file_type}",
                genome_db_dir = genome_db_dir,
                genome_file_type = genome_file_types
                )    
        return genome_db_files

    def checkv_hmm_files():
        hmm_db_dir = path.joincheck_v_db_dir, "hmm_db", "checkv_hmms")
        hmm_index_range = [*range(1,81)]
        hmm_db_files = expand("{hmm_db_dir}/{hmm_index}.hmm", 
                hmm_db_dir = hmm_db_dir, 
                hmm_index = hmm_index_range
                )
        return hmm_db_files
    return check_v_genome_files + checkv_hmm_files

#############
##  Rules  ##
#############
rule all:
    input:
        checkv_out(),
        vs_out(),
         "{data_dir}/CheckM2_database/uniref100.KO.1.dmnd"


# Construct virsorter2 database
rule vs_db:
    output:
        vs_out()
    threads:
        download_threads
    params:
        data_dir = data_dir,
        vs2_db_dir = vs2_db_dir
    conda:
        path.join(env_dir, "virshimeome_base.yml") # vs2
    shell:
        """
        mkdir -p {params.data_dir}
        virsorter setup -d {params.vs2_db_dir} -j {threads} --scheduler greedy
        """
# Construct checkV database. 
rule checkv_db:
    output:
        checkv_out()
    params:
        data_dir = data_dir,
    conda:
        path.join(env_dir, "virshimeome_base.yml") #checkv
    shell:
        """
        checkv download_database {params.data_dir}
        mv {params.data_dir}/checkv-db* {outputcheck_v_db_dir}
        export CHECKVDB={outputcheck_v_db_dir}
        """

rule checkm_db:
    output:
        "{data_dir}/CheckM2_database/uniref100.KO.1.dmnd"
    params:
        data_dir = data_dir
    conda:
        path.join(env_dir, "checkm.yml")
    shell:
        """
        checkm2 database \
            --download  \
            --path {params.data_dir}        

        """