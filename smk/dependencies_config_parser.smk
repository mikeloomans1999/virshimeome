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
if config['virshimeome_dir'] is None:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': virshimeome_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['virshimeome_dir']) is False:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': virshimeome_dir variable path does not exit. Please, complete ', config_path)
else:
    virshimeome_dir=config["virshimeome_dir"]

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