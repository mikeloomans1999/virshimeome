from Bio import Phylo
import os
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
outdir = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/data_visualization/"


temperate = 0 
undefined = 0 
virulent = 0 
lifestyle_dict = {} # dict
with open(phatyp_prediction, "r") as phatype_file:
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

plt.savefig(f"{outdir}phage_lifestyle_piechart.svg")
plt.savefig(f"{outdir}phage_lifestyle_piechart.png")

# Def convert newick seq_ids



tree = Phylo.read()

