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

# VIRSHIMEOME DIR
if config['virshimeome_dir'] is None:
    print('ERROR in ', config_path, ': virshimeome_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['virshimeome_dir']) is False:
    print('WARNING in ', config_path, ': virshimeome_dir path does not exit.')
    virshimeome_dir = config['virshimeome_dir']
else:
    virshimeome_dir = config['virshimeome_dir']

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

# ANI THRESHOLD
if config["ani_threshold"] is None:
    print('ERROR in ', config_path, ': the ani_threshold variable is empty. Please, complete')
    ani_threshold=config["ani_threshold"]
elif type(config["ani_threshold"]) != int and type(config["ani_threshold"]) != float :
    print('ERROR in ', config_path, ': the ani_threshold is not an interger or float. Please, complete', config_path)
    ani_threshold=config["ani_threshold"]
else:
    ani_threshold=config["ani_threshold"]

# PHABOX REJECT
if config["prophage_reject"] is None:
    print('ERROR in ', config_path, ': the prophage_reject variable is empty. Please, complete')
    prophage_reject=config["prophage_reject"]
elif type(config["prophage_reject"]) != int and type(config["prophage_reject"]) != float :
    print('ERROR in ', config_path, ': the prophage_reject is not an interger or float. Please, complete', config_path)
    prophage_reject=config["prophage_reject"]
else:
    prophage_reject=config["prophage_reject"]

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