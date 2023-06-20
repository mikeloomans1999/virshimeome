import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import argparse


def plot_circular_dendrogram(seq_ids, distance_matrix, outfile):
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix_file', type=str, help='Matrix file in PHYLIP format')
    parser.add_argument('--outfile', type=str, help='File for final tree, filetype can be mentioned here ')
    args = parser.parse_args()

    matrix_file = args.matrix_file
    outfile = args.outfile

    # Read distance matrix
    distance_matrix = np.loadtxt(matrix_file)

    # Get sequence IDs
    seq_ids = 
    
    condensed_matrix = squareform(distance_matrix[:, 1:])
    # Perform hierarchical clustering using linkage
    plot_circular_dendrogram(seq_ids, condensed_matrix, outfile)


if __name__ == "__main__":
    main()
