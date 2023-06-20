import numpy as np
from alfpy.utils import seqrecords
from alfpy import word_distance
from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distmatrix
import argparse


def get_seq_record(fasta_file):
    fasta_file_open = open(fasta_file)
    seq_records = seqrecords.read_fasta(fasta_file_open)
    fasta_file_open.close()
    return seq_records
    
def calc_matrix(seq_records, distance_metric):
    p = word_pattern.create(seq_records.seq_list, word_size=2)
    counts = word_vector.Counts(seq_records.length_list, p)
    dist = word_distance.Distance(counts, distance_metric)
    return distmatrix.create(seq_records.id_list, dist)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='Input fasta')
    parser.add_argument('--output', type=str, help='Output distance matrix')
    parser.add_argument('--seq_id_file', type=str, help='seq_id_file')
    parser.add_argument('--metric', type=str, help='Distance estimation method used. ', default="euclid_norm")
    args = parser.parse_args()
    
    #################
    ##  variables  ##
    #################
    fasta_file = args.input
    dist_mat_file = args.output
    distance_metric = args.metric
    seq_id_file = args.seq_id_file
    
    seq_records = get_seq_record(fasta_file)
    matrix = calc_matrix(seq_records, distance_metric)
    id_list = matrix.id_list
    norm_matrix = matrix.normalize()
    np_matrix = matrix.data
    np.savetxt(dist_mat_file, np_matrix)
    
    with open(seq_id_file, "w") as seq_id_open:
        for seq_id in id_list:
            seq_id_open.write(f"{seq_id}\n")

if __name__ == "__main__":
    main()