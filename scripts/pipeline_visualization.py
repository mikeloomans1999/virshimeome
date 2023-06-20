from Bio import Phylo
import argparse
from os import path
import matplotlib.pyplot as plt
from csv import reader

# function calling 
# "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/final-viral-combined.fa", "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/split", "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/split/hello.fa"


# Create a circular taxonomic tree. 
# Species discovered and virulence. 
# Perform accurate blast contigs as the family level blat is based on a limited local database. 
# 

# Virulence
phatyp_prediction = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/8_0_identification/0_filtered_sequences/output/phatyp_prediction.csv"
cherry_host = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/8_0_identification/0_filtered_sequences/output/cherry_prediction.csv"
phagcn_taxonomy = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/8_0_identification/0_filtered_sequences/output/phagcnx_prediction.csv"

outdir = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/data_visualization/"

def get_phage_lifestyle_number(phatyp_prediction_file):
    temperate = 0 
    virulent = 0 
    undefined = 0 
    lifestyle_dict = {} # dict
    with open(phatyp_prediction_file, "r") as phatype_file:
        phatype_file_read = reader(phatype_file, delimiter=",")
        next(phatype_file_read)
        for row in phatype_file_read:
            seq_id, prediction, score = row
            if float(score) > 0.5:
                lifestyle_dict[seq_id] = prediction
                if prediction == "virulent":
                    virulent += 1
                else:
                    temperate += 1
            else:
                lifestyle_dict[seq_id] = "undefined"
                undefined += 1 
    return temperate, virulent, undefined


def write_pie_lifestyle_plot(temperate, virulent, undefined, output_file):
    pie_lifestyle_fig, pie_lifestyle_ax = plt.subplots()

    if undefined > 0:
        pie_lifestyle_ax.pie([virulent, temperate, undefined], 
                            labels=["virulent", "temperate", "undefined"], 
                            autopct=lambda p : '{:.2f}%  ({:,.0f})'.format(p,p * sum([virulent, temperate, undefined])/100)
                            )
    else:
        pie_lifestyle_ax.pie([virulent, temperate], 
                            labels=["virulent", "temperate"], 
                            autopct=lambda p : '{:.2f}%  ({:,.0f})'.format(p,p * sum([virulent, temperate])/100)
                            )

    plt.savefig(output_file)

def convert_newick_seq_ids():
    with open(phagcn_taxonomy, "r") as taxonomy_file:
        taxonomy_file_read = reader(taxonomy_file, delimiter=",")
        next(taxonomy_file_read)
        for row in taxonomy_file_read:
            seq_id, _, taxnomy, score = row
        


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--phatyp_prediction_file', type=str, help='phatyp_prediction_file')
    parser.add_argument('--output_file', type=str, help='data_vizualization_dir')
    args = parser.parse_args()
    
    phatyp_prediction_file = args.phatyp_prediction_file
    output_file = args.output_file
    
    temperate, virulent, undefined = get_phage_lifestyle_number(phatyp_prediction_file)
    write_pie_lifestyle_plot(temperate, virulent, undefined, output_file)

