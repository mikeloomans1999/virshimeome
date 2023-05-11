#!/usr/bin/env python

'''
Download and install dependencies

Authors: Carmen Saenz
'''

# configuration yaml file
from os import path
import glob

config_path = 'configuration yaml file' #args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file: ") #, config_path)
print(" *******************************")
print("  ")

virshimeome_dir=config["virshimeome_dir"]
local_dir=config["local_dir"]
download_threads=config["download_threads"]
download_memory=config["download_memory"]

# Download virsorter2 database and profiles. 
data_dir=directory(path.join(virshimeome_dir,"data"))
pfam_hmm_list = [f"Pfam-A-{hmm_type}.hmm" for hmm_type in ["acc2desc", "Archaea", "Bacteria", "Eukaryota", "Mixed", "Viruses"]]
rbs_list=[f"rbs-catetory{note}.tsv" for note in ["notes", ""]]
group_list=["dsDNAphage","lavidaviridae","NCLDV","RNA", "ssDNA"]
group_files=["hallmark-gene.list", "model"]

# db dirs
pfam_dir_db=path.join(virshimeome_dir, "data", "db", "pfam")
rbs_dir_db=path.join(virshimeome_dir,"data","db","rbs")
group_dir_db=path.join(virshimeome_dir,"data", "db","group")
checkv_dir=path.join(virshimeome_dir, "data", "checkv-db-v1.5")
rule vs_two_db:
    output:
        pfam_hmm_files=expand("{pfam_dir_db}" + os.sep + "{hmm_file}", pfam_dir_db=pfam_dir_db, hmm_file=pfam_hmm_list),
        viral_hmm=path.join(virshimeome_dir,"data","db","hmm", "viral", "combined.hmm"),
        rbs_files=expand("{rbs_dir_db}" + os.sep + "{rbs_file}",rbs_dir_db=rbs_dir_db, rbs_file=rbs_list),
        group_files=expand("{group_dir_db}" + os.sep + "{group}"+ os.sep +"{file}",group_dir_db=group_dir_db, group=group_list, file=group_files)
    params:
        data_dir=data_dir
    threads:
        download_threads
    conda:
        path.join(virshimeome_dir, "conda_env", "virshimeome_base")
    shell:
        """
        virsorter setup -d db -j {threads} {params.data_dir}
        """

# Download checkv database and profiles. 
rule checv_db:
    output:
        check_v_db=expand("{checkv_dir}" + os.sep + "{db_file}",checkv_dir=checkv_dir, db_file=["hmm_db","genome_db"])
    params:
        data_dir=data_dir
    conda:
        path.join(virshimeome_dir, "conda_env", "virshimeome_base")
    shell:
        """
        checkv download_database {params.data_dir}
        """

