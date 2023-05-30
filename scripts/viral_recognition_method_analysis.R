# Load required packages
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)

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
combined_contigs_file <- file.path(output_dir, "combined_contigs.fasta")
checkv_dvf <- file.path(output_dir, "1_2_checkv", "1_1_dvf")
checkv_vs <- file.path(output_dir, "1_2_checkv", "1_1_vs")
checkv_quality_sum_dvf_file <- file.path(checkv_dvf, "quality_summary.tsv")
checkv_quality_sum_vs_file <- file.path(checkv_vs, "quality_summary.tsv")

pre_prediction_contigs <- read.table(combined_contigs_file, header=TRUE, sep="\t", row.names="contig_id")
fas <- readLines(combined_contigs_file)
fas <- fas[!grepl('^$', fas)] # remove empty lines 
(pre_prediction_contig_ids <- gsub('^.+\\|(\\w+)\\|.*$', '>\\1', fas)) # Get seq IDs

checkv_quality_sum_df_dvf <- read.table(checkv_quality_sum_dvf_file, header=TRUE, sep="\t", row.names="contig_id")
checkv_quality_sum_df_vs <- read.table(checkv_quality_sum_vs_file, header=TRUE, sep="\t", row.names="contig_id")

# Compare quality of circular to linear contigs. 


# Make a graph displaying all contig lengths for all. 

# Make a graph displayin contig length of phages for each MAG

# Make a tree displaying the phylogeny of each MAG that has been found. 

# Make a graph displaying the genome quality of each MAG and method. 