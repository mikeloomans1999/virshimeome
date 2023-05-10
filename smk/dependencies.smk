#!/usr/bin/env python

'''
Download and install dependencies

Authors: Carmen Saenz
'''

# configuration yaml file
#import sys
import os.path
from os import path
import glob

#args = sys.argv
#print(args)
#args_idx = sys.argv.index('--configfile')
#print(args_idx)
config_path = 'configuration yaml file' #args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file: ") #, config_path)
print(" *******************************")
print("  ")

# Variables from configuration yaml file
if config['virshimeome_dir'] is None:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': virshimeome_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['virshimeome_dir']) is False:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': virshimeome_dir variable path does not exit. Please, complete ', config_path)
else:
    virshimeome_dir=config["virshimeome_dir"]

if config['local_dir'] is None:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': local_dir variable is empty. Please, complete ', config_path)
else:
    local_dir = config['local_dir']

if config['download_threads'] is None:
    print('ERROR in ', config_path, ': download_threads variable is empty. Please, complete ', config_path)
elif type(config['download_threads']) != int:
    print('ERROR in ', config_path, ': download_threads variable is not an integer. Please, complete ', config_path)
else:
    download_threads=config["download_threads"]


if config['download_memory'] is None:
    print('ERROR in ', config_path, ': download_memory variable is empty. Please, complete ', config_path)
elif type(config['download_memory']) != int:
    print('ERROR in ', config_path, ': download_memory variable is not an integer. Please, complete ', config_path)
else:
    download_memory=config["download_memory"]

# Download virsorter2 database and profiles. 
data_dir=directory(path.join(virshimeome_dir,"data"))
pfam_hmm_list = [f"Pfam-A-{hmm_type}.hmm" for hmm_type in ["acc2desc", "Archaea", "Bacteria", "Eukaryota", "Mixed", "Viruses"]]
rbs_list=[f"rbs-catetory{notes}.tsv" for note in ["notes", ""]]
group_list=["dsDNAphage","lavidaviridae","NCLDV","RNA", "ssDNA"]
group_files=["hallmark-gene.list", "model"]
rule vs_two_db:
    output:
        pfam_hmm_files=expand(path.join(virshimeome_dir, "data", "db","pfam",{hmm_file}), hmm_file=pfam),
        viral_hmm=path.join(virshimeome_dir,"data","db","hmm", "viral", "combined.hmm"),
        rbs_files=expand(path.join(virshimeome_dir,"data","db","rbs",{rbs_file}),rbs_file=rbs_list),
        group_files=expand(path.join(virshimeome_dir,"data", "db","group",{group},{file}), group=group_list, file=group_files)
    params:
        data_dir=data_dir
    threads:
        download_threads
    conda:
        path.join(virshimeome_dir, "envs", "virshimeome.yml")
    shell:
        """
        virsorter setup -d db -j {threads} {parmas.data_dir}
        """

# Download checkv database and profiles. 
rule checv_db:
    output:
        check_v_db=expand(path.join(virshimeome_dir, "data", "checkv-db-v1.5", {db_file}), db_file=["hmm_db","genome_db"])
    params:
        data_dir=data_dir
    conda:
        path.join(virshimeome_dir, "envs", "virshimeome.yml")
    shell:
        """
        checkv download_database {params.data_dir}
        """

