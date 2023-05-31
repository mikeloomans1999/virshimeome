###############
##  imports  ##
###############
import matplotlib.pyplot as plt
import numpy as np
from os import path

# Phylogenetic tree
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqIO

#################
##  variables  ##
#################
output_dir = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/"
checkv_dvf_dir = path.join(output_dir, "1_2_checkv", "1_1_dvf")
checkv_vs_dir = path.join(output_dir, "1_2_checkv", "1_1_vs")

combined_contig_file = path.join(output_dir, "combined_contigs.fasta")
checkv_quality_sum_dvf_file = path.join(checkv_dvf_dir, "quality_summary.tsv")
checkv_quality_sum_vs_file = path.join(checkv_vs_dir, "quality_summary.tsv")

viral_fasta_file_dict = {
    'checkv_dvf': path.join(checkv_dvf_dir, "viruses.fna"),
    'checkv_vs':  path.join(checkv_vs_dir, "viruses.fna"),
    'dvf_prediction': path.join(output_dir, "1_1_dvf", "final-viral.fa"),
    'vs_prediction': path.join(output_dir, "1_1_vs", "final-viral.fa")
}

#################
##  load data  ##
#################
checkv_qs_df_dvf = np.genfromtxt(checkv_quality_sum_dvf_file, delimiter='\t', dtype=str, skip_header=1)
checkv_qs_df_vs = np.genfromtxt(checkv_quality_sum_vs_file, delimiter='\t', dtype=str, skip_header=1)
    
#################
##  functions  ##
#################
# Compare quality of circular to linear contigs. 
def quality_circular(checkv_qs_df):
    quality_circ = []
    quality_norm = []
    for row in checkv_qs_df:
        row_list = row.tolist()
        if "circular" in row_list[0]:
            quality_circ.append(row_list[7])
        else:
            quality_norm.append(row_list[7])
    return quality_circ, quality_norm

def calculate_percentages(quality_list, quality_types):
    total_count = len(quality_list)
    percentage_dict = {}
    for value in quality_types:
        percentage_dict[value] = (quality_list.count(value) / total_count) * 100
    return percentage_dict

def quality_by_contig_length(checkv_qs_df):
    quality_by_contig_length_dict = {}
    
    for row in checkv_qs_df:
        row_list = row.tolist()
        quality_value = row_list[7]
        contig_length = int(row_list[1])
        
        if quality_value not in quality_by_contig_length_dict:
            quality_by_contig_length_dict[quality_value] = []
        
        quality_by_contig_length_dict[quality_value].append(contig_length)
        
    return quality_by_contig_length_dict

def create_phylogenetic_tree(name, viral_fasta_file):
    # Read sequences from the FASTA file
    sequences = SeqIO.parse(viral_fasta_file, 'clustal')

    # Calculate distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(sequences)

    # Construct the tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # Color the tree nodes based on groups
    # for node in tree.get_terminals():
    #     node_name = node.name
        
    #     # Determine the group of the node based on its name
    #     # Adjust this logic according to your specific case
    #     if 'Group A' in node_name:
    #         node.color = group_colors['Group A']
    #     elif 'Group B' in node_name:
    #         node.color = group_colors['Group B']
    #     elif 'Group C' in node_name:
    #         node.color = group_colors['Group C']

    # Plot the tree
    fig, ax = plt.subplots(figsize=(8, 8))
    Phylo.draw(tree, axes=ax, do_show=False)

    # Add a legend for the group colors
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markersize=10)
        # for group, color in group_colors.items()
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    # Show the plot
    plt.savefig(name + "_phylogenetic_tree.png")


#####################
##  Contig length  ##
#####################
contig_lengths_dvf = sorted(checkv_qs_df_dvf[:,1].astype(int).tolist()) # Get all contig lengths 
contig_lengths_vs = sorted(checkv_qs_df_vs[:,1].astype(int).tolist())

# Plot #
contig_length_figure, contig_length_axes = plt.subplots()
contig_length_axes.plot(contig_lengths_dvf, range(0, len(contig_lengths_dvf)))
plt.savefig("contig_length_frequency.png")


########################
##  Percentage viral  ##
########################
len_combined_contig_sequences = 0
combined_contig_sequences = SeqIO.parse(combined_contig_file, 'fasta')
for record in  SeqIO.parse(combined_contig_file, 'fasta'):
    len_combined_contig_sequences += 1
print(len_combined_contig_sequences)
print(len(contig_lengths_vs) / len_combined_contig_sequences)




##############################
##  Quality by contig type  ##
##############################
# DVF
quality_circ_dvf, quality_norm_dvf = quality_circular(checkv_qs_df=checkv_qs_df_dvf)
# VS
quality_circ_vs, quality_norm_vs = quality_circular(checkv_qs_df=checkv_qs_df_vs)

names = ["dvf circular", "vs2 circular"] #, "dvf linear", "vs2 linear"]
quality_types = ['Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete']

###### ABSOLUTE ######

absolute_weight_dict = {}
for quality in quality_types:
    
    absolute_weight_dict[quality] = np.array([
            quality_circ_dvf.count(quality), 
            quality_circ_vs.count(quality)]) 
            #, norm_dvf_quality.count(quality)
            #, norm_vs_quality.count(quality)])

circular_absolute_figure, circular_absolute_axes = plt.subplots()
bottom = np.zeros(2)

for quality, weight_absolute in absolute_weight_dict.items():
    p = circular_absolute_axes.bar(names, weight_absolute, 0.5, label=quality, bottom=bottom)
    bottom += weight_absolute

plt.grid('on', color = 'gainsboro', axis = 'y', linestyle='--', zorder=3)
circular_absolute_axes.legend(loc="upper right")
plt.savefig("circular_quality_absolute_bar.png")

###### PERCENTAGES #######
# Circular
quality_circ_dvf_perc = calculate_percentages(quality_circ_dvf, quality_types)
quality_circ_vs_perc = calculate_percentages(quality_circ_vs, quality_types)
# Linear
# quality_norm_dvf_perc = calculate_percentages(quality_norm_dvf, quality_types)
# quality_norm_vs_perc = calculate_percentages(quality_norm_vs, quality_types)

weight_percentages = {}
for quality in quality_types:
    circ_dvf_quality_percentage =  quality_circ_dvf_perc.get(quality)
    circ_vs_quality_percentage = quality_circ_vs_perc.get(quality)

    # norm_dvf_quality = quality_norm_dvf_perc.get(quality)
    # norm_vs_quality = quality_norm_vs_perc.get(quality)
    
    weight_percentages[quality] = np.array([circ_dvf_quality_percentage, circ_vs_quality_percentage]) #, norm_dvf_quality, norm_vs_quality])

## Plotting ##
circular_percentage_figure, circular_percentage_axes = plt.subplots()
bottom = np.zeros(2)

for quality, weight_percentage in weight_percentages.items():
    p = circular_percentage_axes.bar(names, weight_percentage, 0.5, label=quality, bottom=bottom)
    bottom += weight_percentage

plt.grid(color = 'gainsboro', axis = 'y', linestyle='--', zorder=3)
circular_percentage_axes.legend(loc="upper right")
plt.savefig("circular_quality_percentage_bar.png")

################################
##  Quality by contig length  ##
################################
quality_by_contig_length_dvf_dict = quality_by_contig_length(checkv_qs_df_dvf)
quality_by_contig_length_vs_dict = quality_by_contig_length(checkv_qs_df_vs)

# Get average and stdev per quality type and plot those in barplots. 
values_dvf = [quality_by_contig_length_dvf_dict.get(quality_type) for quality_type in quality_types]
values_vs = [quality_by_contig_length_vs_dict.get(quality_type) for quality_type in quality_types]

contig_lengths_figure, contig_lengths_axes = plt.subplots(figsize=(10,10))
contig_lengths_axes.set_ylabel("contig length")
contig_lengths_axes.set_xlabel("contig quality")
contig_lengths_axes.boxplot(values_dvf)
contig_lengths_axes.set_xticklabels(quality_types)
plt.xticks(rotation=35)

plt.savefig("contig_quality_boxplot.png")

# Make a tree displaying the phylogeny of each MAG that has been found. 

#########################
##  Phylogenetic tree  ##
#########################

# for data_input_name, fasta_file_path in msa_file_dict.items():
#     create_phylogenetic_tree(data_input_name, fasta_file_path)
