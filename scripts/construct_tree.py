
import matplotlib.colors as mcolors
import random
import numpy as np
import pandas  as pd
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
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
    
def plot_colorful_dendrogram(df, outfile, feature_dict):
    # Create a list of colors for each feature
    unique_features = sorted(set(feature_dict.values()))
    color_map = plt.cm.get_cmap('tab10', len(unique_features))
    colors = [color_map(i) for i in range(len(unique_features))]
    feature_colors = {feature: color for feature, color in zip(unique_features, colors)}

    # Perform hierarchical clustering on the DataFrame
    df = df[list(feature_dict.keys())]
    linked = linkage(df.transpose(), method='average')

    # Create the dendrogram plot
    plt.figure(figsize=(8, 6))
    dendrogram(linked, orientation='top', labels=df.columns, color_threshold = 0, link_color_func = lambda k: 'black' )

    # Color the branches based on the feature_dict feature
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        sample = lbl.get_text()
        feature = feature_dict.get(sample)
        color = feature_colors.get(feature, 'gray')
        lbl.set_color(color)


    # Create a legend
    unique_families = sorted(set(feature_dict.values()))
    color_map = plt.cm.get_cmap('tab10', len(unique_families))
    family_colors = {family: color_map(i) for i, family in enumerate(unique_families)}
    unique_families_frequency = [f"{family} ({list((feature_dict.values())).count(family)})" for family in unique_families]
    handles = [plt.Line2D([0], [0], color=color, lw=2, linestyle='--') for color in family_colors.values()]
    labels = unique_families_frequency
    plt.legend(handles, labels, loc='upper right')

    # Set the axis labels and title
    plt.xlabel('MAGs')
    plt.ylabel('Distance')
    plt.title('Dendrogram with feature_dict-Based Coloring')

    # Show the plot
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

    if feature_file != "":    
        pd_df_matrix = pd.DataFrame(pd.read_csv(matrix_file, sep="\t", index_col = 0))
        pd_df = pd.read_csv(feature_file)
        feature_dict = convert_two_arrays_to_dict(pd_df["Accession"], pd_df["Pred"])
        plot_colorful_dendrogram(pd_df_matrix, outfile, feature_dict)
    else:
        # Read distance matrix
        distance_matrix = np.array(np.loadtxt(matrix_file, delimiter = "\t", dtype=object)[:,1:], dtype=float)
        print(distance_matrix.shape)

        # Get sequence IDs
        seq_ids = np.loadtxt(matrix_file, usecols=0, dtype=str).tolist()
        print(len(seq_ids))
        
        condensed_matrix = squareform(distance_matrix)
        # Perform hierarchical clustering using linkage
        plot_dendrogram(seq_ids, condensed_matrix, outfile)


if __name__ == "__main__":
    main()
