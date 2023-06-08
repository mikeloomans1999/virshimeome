# Load required packages
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
librayr(phyloseq)
output_dir <- "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/"

# Define command-line arguments
# args <- commandArgs(trailingOnly=TRUE)
# if (length(args) != 3) {
#   stop("Please provide two arguments: input matrix file path and output directory path")
# }
# matrix_file <- args[1] # Location of the relative abundance matrix file
# output_dir <- args[2] # Code before is guaranteed to create the output directory, 
# beta_diversity_metric <- args[3] # Distance metrix for phyloseq all possible metrics are checked for in the config file. 


# Load the TSV matrix file as a data frame
# df <- read.table(matrix_file, header=TRUE, sep="\t", row.names="ID")
# row_names <- c("ID", rownames(df))
# sample_ids <- colnames(df)
# rownames(df) <- NULL # Remove rownames to avoid issues later on

# # Convert the data frame to a matrix
# mat <- data.matrix(df)

# # Create a phyloseq object from the matrix
# ps <- phyloseq(otu_table(mat, taxa_are_rows=TRUE))
#################
##  variables  ##
#################
