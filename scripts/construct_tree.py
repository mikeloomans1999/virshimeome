import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

 
def get_seq_ids(fasta_file):
    seq_ids = set()
    with open(fasta_file, "r") as fasta_open:
        for line in fasta_open:
            if line.startswith(">"):
                seq_id = line.strip().lstrip(">")
                seq_ids.add(seq_id)
    return seq_ids
   

def plot_circular_dendrogram(seq_ids, linked_clustering, outfile):
    # Plot the dendrogram in circular form
    fig, ax = plt.subplots(figsize=(8, 8))
    dendrogram(linked_clustering, orientation='right', labels=seq_ids, color_threshold=0, leaf_font_size=12)

    # Set the axis labels and title
    ax.set_xlabel('Similarity')
    ax.set_ylabel('Species')
    ax.set_title('Circular Cladogram')

    # Adjust the position and rotation of the labels
    plt.yticks(rotation=90, va='center')

    # Set the aspect ratio to be equal to make it circular
    ax.set_aspect('equal', 'box')

    # Remove the ticks and spines
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Save the plot
    plt.savefig(outfile)


def main():
    matrix_file = "/projects/arumugam/scratch/mnc390/virome_testing/ani_virshimeome/9_0_distance/distance.mat"
    fasta_file = "/projects/arumugam/scratch/mnc390/virome_testing/ani_virshimeome/2_checkv/1_1_dvf/viruses.fna"
    outfile =  "/projects/arumugam/scratch/mnc390/virome_testing/ani_virshimeome/data_visualization/circular_cladogram.svg"
    # Sample similarity matrix
    similarity_matrix = np.genfromtxt(fname=matrix_file, delimiter="\t", skip_header=1, filling_values=1, )

    seq_ids = get_seq_ids(fasta_file)
    
    # Perform hierarchical clustering using linkage
    linked_clustering = linkage(similarity_matrix, method='complete')
    plot_circular_dendrogram(seq_ids, linked_clustering, outfile)
    

if __name__ == "__main__":
    main()