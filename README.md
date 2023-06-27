
Setting up the initial conda environment:
```console
git clone https://github.com/mikeloomans1999/virshimeome
cd virshimeome
conda env create -f scripts/virshimeome.yml
```

Configure the config files.
Configuration for the dependencies config file, dependencies.yml.

```
######################
# General settings
######################
local_dir: # Local directory for storing the results. 
virshimeome_dir: # Local virshimeome directory

download_memory: # In gigatbytes
download_threads: # Number of threads available for downloads

```

Configuration for the pipeline execution, main.yml.
```
######################
## General settings ##
######################
# Dirs
PROJECT: # Project name
working_dir: /projects/arumugam/scratch/mnc390/virome_testing/data/test_data_virshimeome
local_dir: /projects/arumugam/scratch/mnc390/virome_testing/data/test_data_virshimeome
virshimeome_dir: # Local virshimeome directory 
hq_reads_dir: # Reads directory, not currently in use
contig_dir: # parent directory of contigs/MAGs
output_dir: # Your assigned output directory
deep_vir_finder_dir: # Local DeepVirFinder directory
# Resources 
memory: # In gigabytes

#####################
##  Tool settings  ##
#####################
# Contig selection (custom script & virsorter2)
circular: circular # Whether or not you want to filter for circular sequences or any other characteristics mentioned in the sequence identifier. 
max_contig_length: 500000
min_contig_length: 2000
min_score_vir_recognition: 0.5

# DVF
min_score_dvf: 0.7
max_pval_dvf: 0.05

# Read to contig alignment (BWA)
algorithm_bwa: mem # Currently not in use. 

# Relative abundance (samtools & msamtools)
min_length_alignment_abundance: 80  # Currently not in use. 
ani_threshold: 99.9 # Percentage not decimal  # Currently not in use. 
percentage_read_aligned_abundance: 80 # Percentage not decimal  # Currently not in use. 
```

Downloading all the databases
```console
snakemake --snakefile smk/dependencies.smk --configfile config/dependencies.yaml --use-conda
```

Run the pipeline, it is recommended to use the "-j" flag and specify multiple cores to use.
```console
snakemake --snakefile smk/main.smk --configfile config/main.yaml --use-conda
```
