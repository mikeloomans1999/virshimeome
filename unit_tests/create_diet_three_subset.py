import glob
import os
import random

# Variables
diet_three_dir = "/projects/arumugam/scratch/mjq180/Diet3_for_Mike/**/*"
outdir = "/projects/arumugam/scratch/mnc390/virome_testing/data/mag_subset/"

# Get all mags present. 
all_mags = glob.glob(diet_three_dir)
n_mags = len(all_mags)
percentage = 0.05
n_mags_chosen = int(n_mags * percentage)

print(n_mags_chosen)

# Get files
chosen_files = random.sample(all_mags, k=n_mags_chosen)

# Create new dir with all this stuff in symlinks.
for file_path in chosen_files:
    new_mag_dir = os.path.join(outdir, file_path.split(os.sep)[-2])
    
    if not os.path.exists(new_mag_dir):
        os.makedirs(new_mag_dir)
            
    filename = os.path.basename(file_path)
    sym_file = os.path.join(new_mag_dir, filename)
    os.symlink(file_path, sym_file)