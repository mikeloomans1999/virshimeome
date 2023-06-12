
import os
from alfpy.utils import seqrecords, distmatrix
from alfpy.utils.data import subsmat
from alfpy import ncd, wmetric, word_distance, word_pattern, word_vector


if not os.path.exists("/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/split"):
    os.makedirs("/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/split")


fh = open("/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/final-viral-combined.fa")
seq_records = seqrecords.read_fasta(fh)
fh.close()  
print(seq_records)
dist = ncd.Distance(seq_records)
matrix = distmatrix.create(seq_records.id_list, dist)
matrix.display()

w_matrix = subsmat.get('blosum62')


p = word_pattern.create(seq_records.seq_list, word_size=2)
counts = word_vector.Counts(seq_records.length_list, p)
dist = word_distance.Distance(counts, 'euclid_norm')
matrix = distmatrix.create(seq_records.id_list, dist)
matrix.display()

# function calling 
# "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/final-viral-combined.fa", "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/split", "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/split/hello.fa"