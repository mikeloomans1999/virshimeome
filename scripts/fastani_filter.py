import networkx as nx
import argparse

def get_ids(seq_id_file):
    # Get all ids
    with open(seq_id_file, "r") as seq_id_no_filter:
        return set(seq_id_no_filter.read().splitlines())

def cluster_duplicates(fastani_out):
    # Read the input file and create an undirected graph
    graph = nx.Graph()
    with open(fastani_out, "r") as file:
        for line in file:
            id_1, id_2, score,_,_ = line.split("\t")

            if float(score) >= 99.9:
                graph.add_edge(id_1, id_2, score=float(score))

    # Perform clustering based on connected components    
    return list(nx.connected_components(graph))

def get_duplicate_ids(fastani_out):
    possible_duplicate_ids = set()
    with open(fastani_out, "r") as file:
        for line in file:
            line_split = line.split("\t")
            possible_duplicate_ids.add(line_split[0])
            possible_duplicate_ids.add(line_split[1])
    return possible_duplicate_ids
            
def select_longest_contig(clusters):
    # Get longest sequence from each cluster
    longest_identical_ids = set()
    for cluster in clusters:
        longest_contig = 0
        for seq_id in cluster:
            seq_id_list = seq_id.split("_")
            contig_length = int(seq_id_list[seq_id_list.index("length") + 1])
            if contig_length > longest_contig:
                longest_id = seq_id
                longest_contig = contig_length
        longest_identical_ids.add(longest_id)
    return longest_identical_ids

def convert_set_to_file(set_input, output_file):
    with open(output_file, "a") as new_file:
        for seq_id in set_input:
            new_file.write(f"{seq_id}\n")

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--fastani_out', type=str, help='pairwise comparison average nucleotide identity ')
    parser.add_argument('--revised_seq_ids', type=str, help='output file with selected seq ids')
    parser.add_argument('--seq_id_file', type=str, help='seq id file with duplicates ')
    args = parser.parse_args()
    
    fastani_out = args.fastani_out
    revised_seq_ids = args.revised_seq_ids
    seq_id_file = args.seq_id_file
    
    
    all_ids = get_ids(seq_id_file)
    clusters = cluster_duplicates(fastani_out)
    possible_duplicate_ids = get_duplicate_ids(fastani_out)
    longest_identical_ids = select_longest_contig(clusters)
    
    all_ids -= possible_duplicate_ids
    revised_ids = all_ids.union(longest_identical_ids)
    convert_set_to_file(revised_ids, revised_seq_ids)

if __name__ == "__main__":
        main()