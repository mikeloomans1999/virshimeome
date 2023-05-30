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
    print('WARNING in ', config_path, ': contig_dir path does not exit. ')
    contig_dir = config['contig_dir']
else:
    contig_dir = config['contig_dir']

# DEEP VIR FINDER DIR
if config['deep_vir_finder_dir'] is None:
    print('ERROR in ', config_path, ': the deep_vir_finder_dir variable is empty. Please, complete')
    deep_vir_finder_dir = config['deep_vir_finder_dir']
elif type(config['deep_vir_finder_dir']) != str:
    print('ERROR in ', config_path, ': the deep_vir_finder_dir is not a string. Please, complete', config_path)
    deep_vir_finder_dir = config['deep_vir_finder_dir']
elif path.exists(config['deep_vir_finder_dir']) is False:
    print('WARNING in ', config_path, ': deep_vir_finder_dir path does not exit.')
    deep_vir_finder_dir = config['deep_vir_finder_dir']
else:
    deep_vir_finder_dir = config['deep_vir_finder_dir']

# OUTPUT DIR
if config['output_dir'] is None:
    print('ERROR in ', config_path, ': the output_dir variable is empty. Please, complete')
    output_dir = config['output_dir']
elif type(config['output_dir']) != str:
    print('ERROR in ', config_path, ': the output_dir is not a string. Please, complete', config_path)
    contig_dir = config['output_dir']
elif path.exists(config['output_dir']) is False:
    print('WARNING in ', config_path, ': output_dir path does not exit. The directory will be created.')
    os.mkdir(config["output_dir"])
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
elif type(config["memory"]) != int:
    print('ERROR in ', config_path, ': the memory is not an interger. Please, complete', config_path)
    memory=config["memory"]
else:
    memory=config["memory"]

# RELATIVE ABUNDANCE MIN LENGTH READ ALIGNMENT
if config["min_length_alignment_abundance"] is None:
    print('ERROR in ', config_path, ': the min_length_alignment_abundance variable is empty. Please, complete')
    min_length_alignment_abundance=config["min_length_alignment_abundance"]
elif type(config["min_length_alignment_abundance"]) != int:
    print('ERROR in ', config_path, ': the min_length_alignment_abundance is not an interger. Please, complete', config_path)
    min_length_alignment_abundance=config["min_length_alignment_abundance"]
else:
    min_length_alignment_abundance=config["min_length_alignment_abundance"]

# PERCENTAGE IDENTITY FOR ABUNDANCE     
if config["percentage_identity_abundance"] is None:
    print('ERROR in ', config_path, ': the percentage_identity_abundance variable is empty. Please, complete')
    percentage_identity_abundance=config["percentage_identity_abundance"]
elif type(config["percentage_identity_abundance"]) != int and type(config["percentage_identity_abundance"]) != float :
    print('ERROR in ', config_path, ': the percentage_identity_abundance is not an interger or float. Please, complete', config_path)
    percentage_identity_abundance=config["percentage_identity_abundance"]
else:
    percentage_identity_abundance=config["percentage_identity_abundance"]

# PERCENTAGE OF READ ALIGNED FOR ABUNDANCE
if config["percentage_read_aligned_abundance"] is None:
    print('ERROR in ', config_path, ': the percentage_read_aligned_abundance variable is empty. Please, complete')
    percentage_read_aligned_abundance=config["percentage_read_aligned_abundance"]
elif type(config["percentage_read_aligned_abundance"]) != int and type(config["percentage_read_aligned_abundance"]) != float :
    print('ERROR in ', config_path, ': the percentage_read_aligned_abundance is not an interger or float. Please, complete', config_path)
    percentage_read_aligned_abundance=config["percentage_read_aligned_abundance"]
else:
    percentage_read_aligned_abundance=config["percentage_read_aligned_abundance"]

# CIRCULAR
if type(config["circular"]) != str:
    print('ERROR in ', config_path, ': the config is not string. Please, complete', config_path)
    circular=config["circular"]
else:
    circular=config["circular"]

# MIN SCORE DVF
if config["min_score_dvf"] is None:
    print('ERROR in ', config_path, ': the min_score_dvf variable is empty. Please, complete')
    min_score_dvf=config["min_score_dvf"]
elif type(config["min_score_dvf"]) != float :
    print('ERROR in ', config_path, ': the min_score_dvf is not an float. Please, complete', config_path)
    min_score_dvf=config["min_score_dvf"]
else:
    min_score_dvf=config["min_score_dvf"]

# MAX PVALUE DVF
if config["max_pval_dvf"] is None:
    print('ERROR in ', config_path, ': the max_pval_dvf variable is empty. Please, complete')
    max_pval_dvf=config["max_pval_dvf"]
elif type(config["max_pval_dvf"]) != float :
    print('ERROR in ', config_path, ': the max_pval_dvf is not an float. Please, complete', config_path)
    max_pval_dvf=config["max_pval_dvf"]
else:
    max_pval_dvf=config["max_pval_dvf"]