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
if config['minto_dir'] is None:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': minto_dir variable is empty. Please, complete ', config_path)
elif path.exists(config['minto_dir']) is False:
    #print('ERROR in ')
    print('ERROR in ', config_path, ': minto_dir variable path does not exit. Please, complete ', config_path)
else:
    minto_dir=config["minto_dir"]

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

