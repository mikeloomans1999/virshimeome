
Setting up the initial conda environment:
```console
git clone https://github.com/mikeloomans1999/virshimeome
```
```console
cd virshimeome
```
```console
conda env create -f envs/virshimeome_base.yml
```

Configure the config files.
Configuration for the dependencies config file, config/dependencies.yml.

``` yaml
######################
# General settings
######################
local_dir: # Local directory for storing the results. 
virshimeome_dir: # Local virshimeome directory

download_memory: # In gigatbytes
download_threads: # Number of threads available for downloads

```

Configuration for the pipeline execution, main.yml.
``` yaml    
######################
## General settings ##
######################
# Dirs
PROJECT: # Project name
virshimeome_dir: # Virshimeome directory
contig_dir: # Parent directory of all contigs
output_dir: 
deep_vir_finder_dir: # Local DeepVirFinder repository

#####################
##  Tool settings  ##
#####################
# Contig selection (custom script & virsorter2) # Filtration parameter
circular:  # Special characters or words mentioned in sequence IDs
max_contig_length: 
min_contig_length: 

# VS2
min_score_vir_recognition: 0.5

# DVF
min_score_dvf: 0.7 # Combination of dvf_score and p-value used as a threshold for DeepVirFinder phage prediction.
max_pval_dvf: 0.05 

# phabox
prophage_reject: 0.2 # The reject threshold for the PhaBOX tools

# ANI
ani_threshold: 99.9 # Percentage not decimal
```

Downloading all the databases
```console
snakemake --snakefile smk/dependencies.smk --configfile config/dependencies.yaml --use-conda
```

Run the pipeline, it is recommended to use the "-j" flag and specify multiple cores to use.
```console
snakemake --snakefile smk/main.smk --configfile config/main.yaml --use-conda
```
Known issues:
There is a nonexistent dependency clash when installing blast+ >2.10.1, earlier version are unable to run  blastn in redhat distros past version 8.1 due to the absence of the libsnl library. 
To fix:
```console
conda activate phabox 
wget https://anaconda.org/bioconda/blast/2.14.0/download/linux-64/blast-2.14.0-h7d5a4b4_1.tar.bz2 ./
conda install blast-2.14.0-h7d5a4b4_1.tar.bz2
```
Run the pipeline in greedy mode if  there are any workload manager issues.
```console
snakemake --snakefile smk/main.smk --configfile config/main.yaml --use-conda --scheduler greedy
```

