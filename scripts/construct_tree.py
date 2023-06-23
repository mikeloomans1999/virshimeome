
import matplotlib.colors as mcolors
import random
import numpy as np
import pandas  as pd
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def get_seq_ids(seq_id_file):
    with open(seq_id_file, "r") as seq_id_open:
        lines = seq_id_open.readlines()
        return [line.strip() for line in lines]
    
def convert_two_arrays_to_dict(array_keys, array_values):
    dict_from_arrays = {}
    key_list = array_keys.tolist()
    for index, key in enumerate(key_list):
        dict_from_arrays[key] = array_values[index]
    return dict_from_arrays

def plot_dendrogram(seq_ids, distance_matrix, outfile):
    fig, ax = plt.subplots(figsize=(8, 8))
    linked_clustering = linkage(distance_matrix, method='complete')
    dendrogram(linked_clustering, orientation='top', labels=seq_ids, color_threshold=0, leaf_font_size=12)

    # Set the axis labels and title
    ax.set_xlabel('Species')
    ax.set_ylabel('Distance')
    ax.set_title('Dendrogram')

    # Adjust the position and rotation of the labels
    plt.xticks(rotation=90, ha='center')

    # Remove the spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Save the plot
    plt.savefig(outfile)
    
    
def plot_colorful_dendrogram(seq_ids, distance_matrix, outfile, feature_dict):
    
    # Add color to lineages
    lineage_colors = []
    unique_features = set(feature_dict.values())
    unique_features.add("unknown")
    unique_features = list(unique_features)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    linked_clustering = linkage(distance_matrix, method='complete')
    id_list_dendrogram = leaves_list(linked_clustering).tolist()
    
    # Get a list of color names
    
    # Convert colors to string representation
    color_names = list(mcolors.CSS4_COLORS.keys())    
    # Select the colors
    color_strings = color_names[:(len(seq_ids)+1)]
    
    colors = ["b"]*(2*len(seq_ids)-1)
    
    for index, seq_id in enumerate(seq_ids):
        # Get feature
        if seq_id in feature_dict:
            feature = feature_dict[seq_id]
        else:
            feature = "unknown"
        # Get color
        color = color_strings[unique_features.index(feature)]
        # Store color
        colors[id_list_dendrogram.index(index)] = color
    
    
    print(len(seq_ids), "/", len(lineage_colors))
    dendrogram(linked_clustering, orientation='top', labels=seq_ids,color_threshold=0.2, leaf_font_size=12, link_color_func=lambda k: colors[k])
    
    # Set the axis labels and title
    ax.set_xlabel('Species')
    ax.set_ylabel('Distance')
    ax.set_title('Dendrogram')

    # Adjust the position and rotation of the labels
    plt.xticks(rotation=90, ha='center')

    # Remove the spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Create a legend
    # handles = [plt.Line2D([0], [0], color=color, lw=2) for color in lineage_colors.values()]
    # labels = lineage_colors.keys()
    # ax.legend(handles, labels, loc='upper right')
    
    # Save the plot
    plt.savefig(outfile)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix_file', type=str, help='Matrix file in PHYLIP format')
    parser.add_argument('--outfile', type=str, help='File for final tree, filetype can be mentioned here ')
    parser.add_argument('--feature_file', type=str, help='File with all the sequence IDs in the "Accession" column and their respective characteristic in the "Pred" column in .csv format.  ', default="")
    args = parser.parse_args()

    matrix_file = args.matrix_file
    outfile = args.outfile
    feature_file = args.feature_file
    
    # Read distance matrix
    # distance_matrix = np.loadtxt(matrix_file)
    distance_matrix = np.array(np.loadtxt(matrix_file, delimiter = "\t", dtype=object)[:,1:], dtype=float)
    print(distance_matrix.shape)

    # Get sequence IDs
    seq_ids = np.loadtxt(matrix_file, usecols=0, dtype=str).tolist()
    print(len(seq_ids))
    
    condensed_matrix = squareform(distance_matrix)
    
    if feature_file != "":    
        pd_df = pd.read_csv(feature_file)
        feature_dict = convert_two_arrays_to_dict(pd_df["Accession"], pd_df["Pred"])
        plot_colorful_dendrogram(seq_ids, condensed_matrix, outfile, feature_dict)
    else:
        # Perform hierarchical clustering using linkage
        plot_dendrogram(seq_ids, condensed_matrix, outfile)


if __name__ == "__main__":
    main()
