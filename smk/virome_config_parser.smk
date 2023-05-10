#!/usr/bin/env python

'''
Reading of beta and alpha diversity paramters. 
Authors: Carmen Saenz, Mani Arumugam
'''
from os import path

config_path = 'configuration yaml file' # + args[args_idx+1]
print(" *******************************")
print(" Reading configuration yaml file: ", config_path)
print(" *******************************")

def check_strings_in_list(str_list, main_str):
    for string in str_list:
        if string not in main_str:
            return False
    return True

# Variables from configuration yaml file
# PROJECT NAME
if config['PROJECT'] is None:
    print('ERROR in ', config_path, ': PROJECT variable is empty. Please, complete ', config_path)
else:
    project_id = config['PROJECT']

# WORKING DIR PATH
if config['working_dir'] is None:
    print('ERROR in ', config_path, ': working_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['working_dir']) is False:
    print('WARNING in ', config_path, ': working_dir path does not exit.')
    working_dir = config['working_dir']
else:
    working_dir = config['working_dir']
# VIRSHIMEOME DIR
if config['virshimeome_dir'] is None:
    print('ERROR in ', config_path, ': virshimeome_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['virshimeome_dir']) is False:
    print('WARNING in ', config_path, ': virshimeome_dir path does not exit.')
    virshimeome_dir = config['virshimeome_dir']
else:
    virshimeome_dir = config['virshimeome_dir']


# HIGH QUALITY READS DIR
if config['hq_reads_dir'] is None:
    print('ERROR in ', config_path, ': hq_reads_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['hq_reads_dir']) is False:
    print('WARNING in ', config_path, ': hq_reads_dir path does not exit.')
    hq_reads_dir = config['hq_reads_dir']
else:
    hq_reads_dir = config['hq_reads_dir']

# THREADS
if config['threads'] is None:
    print('ERROR in ', config_path, ': threads variable is empty. Please, complete ', config_path)
elif type(config['threads']) != int:
    print('ERROR in ', config_path, ': threads variable is not an integer. Please, complete ', config_path)

# LOCAL DIR
if config['local_dir'] is None:
    print('ERROR in ', config_path, ': local_dir variable is empty. Please, complete ', config_path)
else:
    local_dir = config['local_dir']

# MAX CONTIG LENGTH
if config['max_contig_length'] is None:
    print('ERROR in ', config_path, ': the max_contig_length variable is empty. Please, complete')
    max_contig_length = config['max_contig_length']
elif type(config['max_contig_length']) != int:
    print('ERROR in ', config_path, ': the max_contig_length is not an interger. Please, complete', config_path)
    max_contig_length = config['max_contig_length']
else:
    max_contig_length = config['max_contig_length']

# MIN CONTIG LENGTH
if config['min_contig_length'] is None:
    print('ERROR in ', config_path, ': the min_contig_length variable is empty. Please, complete')
    min_contig_length = config['min_contig_length']
elif type(config['min_contig_length']) != int:
    print('ERROR in ', config_path, ': the min_contig_length is not an interger. Please, complete', config_path)
    min_contig_length = config['min_contig_length']
else:
    min_contig_length = config['min_contig_length']

# MIN SCORE VIR RECOGNITION
if config['min_score_vir_recognition'] is None:
    print('ERROR in ', config_path, ': the min_score_vir_recognition variable is empty. Please, complete')
    min_score_vir_recognition = config['min_score_vir_recognition']
elif type(config['min_score_vir_recognition']) != int and type(config['min_score_vir_recognition']) != float:
    print('ERROR in ', config_path, ': the min_score_vir_recognition is not an interger or float. Please, complete', config_path)
    min_score_vir_recognition = config['min_score_vir_recognition']
else:
    min_score_vir_recognition = config['min_score_vir_recognition']

# CONTIG DIR
if config['contig_dir'] is None:
    print('ERROR in ', config_path, ': the contig_dir variable is empty. Please, complete')
    contig_dir = config['contig_dir']
elif type(config['contig_dir']) != str:
    print('ERROR in ', config_path, ': the contig_dir is not a string. Please, complete', config_path)
    contig_dir = config['contig_dir']
elif path.exists(config['contig_dir']) is False:
    print('WARNING in ', config_path, ': contig_dir path does not exit. The directory will be created.')
    contig_dir = config['contig_dir']
else:
    contig_dir = config['contig_dir']

# OUTPUT DIR
if config['output_dir'] is None:
    print('ERROR in ', config_path, ': the output_dir variable is empty. Please, complete')
    output_dir = config['output_dir']
elif type(config['output_dir']) != str:
    print('ERROR in ', config_path, ': the output_dir is not a string. Please, complete', config_path)
    contig_dir = config['output_dir']
elif path.exists(config['output_dir']) is False:
    print('WARNING in ', config_path, ': output_dir path does not exit. The directory will be created.')
    output_dir = config['output_dir']
else:
    output_dir = config['output_dir']

# BWA ALGORITHM
if config['algorithm_bwa'] is None:
    print('ERROR in ', config_path, ': the algorithm_bwa variable is empty. Please, complete')
    algorithm_bwa=config["algorithm_bwa"]
elif config["algorithm_bwa"] not in ["ont2d", "intractg","sr","spliced", "mem"]:
    print('ERROR in ', config_path, ': the algorithm_bwa variable does not contain an available algorithm choose from: "ont2d", "intractg","sr","spliced", "mem" ')
    algorithm_bwa=config["algorithm_bwa"]
else:
    algorithm_bwa=config["algorithm_bwa"]

# MEMORY ALLOCATION
if config["memory"] is None:
    print('ERROR in ', config_path, ': the memory variable is empty. Please, complete')
    memory=config["memory"]
elif type(config["memory"]) is int:
    print('ERROR in ', config_path, ': the memory is not an interger. Please, complete', config_path)
    memory=config["memory"]
else:
    memory=config["memory"]

# CORE ALLOCATION
if config["cores"] is None:
    print('ERROR in ', config_path, ': the cores variable is empty. Please, complete')
    cores=config["cores"]
elif type(config["cores"]) is int:
    print('ERROR in ', config_path, ': the cores is not an interger. Please, complete', config_path)
    cores=config["cores"]
else:
    cores=config["cores"]