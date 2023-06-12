# Load required packages
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
output_dir <- "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/"

#################
##  variables  ##
#################

mat <- scan('A.txt')
