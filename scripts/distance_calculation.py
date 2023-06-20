from alfpy.utils import seqrecords
from alfpy import word_distance
from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distmatrix
import argparse

def calc_matrix(fasta_file, dist_mat_file, distance_metric):
    fasta_file_open = open(fasta_file)
    seq_records = seqrecords.read_fasta(fasta_file_open)
    fasta_file_open.close()

    p = word_pattern.create(seq_records.seq_list, word_size=2)
    counts = word_vector.Counts(seq_records.length_list, p)
    dist = word_distance.Distance(counts, distance_metric)
    matrix = distmatrix.create(seq_records.id_list, dist)

    out_file = open(dist_mat_file, 'w')
    matrix.write_to_file(out_file)
    out_file.close()

# t 10 --fasta /projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/0_filtered_sequences/combined_sequences.fa --out /projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/0_filtered_sequences/distane.mat

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='Input fasta', default="/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/0_filtered_sequences/combined_sequences.fa")
    parser.add_argument('--output', type=str, help='Output distance matrix', default= "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/0_filtered_sequences/combined_sequences_distances.mat")
    parser.add_argument('--metric', type=str, help='Distance estimation method used. ', default="euclid_norm")
    args = parser.parse_args()
    
    #################
    ##  variables  ##
    #################
    fasta_file = args.input
    dist_mat_file = args.output
    distance_metric = args.metric
    
    calc_matrix(fasta_file, dist_mat_file, distance_metric)
    
if __name__ == "__main__":
    main()