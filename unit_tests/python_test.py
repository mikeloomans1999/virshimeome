from os import path, mkdir
import glob
from itertools import groupby

def read_fasta_file(fasta_file):
	with open(fasta_file, 'rb') as fasta_file_read:
		faiter = (x[1] for x in groupby(fasta_file_read, lambda line: str(line, 'utf-8')[0] == ">"))
		for header in faiter:
			seq_id = str(header.__next__(), 'utf-8').replace("\n", "") 
			seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__()) + "\n"
			yield seq_id, seq


file_path="/projects/arumugam/scratch/mnc390/virome_testing/data/vs2_test_data_contigs/test.fa"
result = ""
for seq_id, seq in read_fasta_file(file_path):
	result += seq_id + "\n" + seq

