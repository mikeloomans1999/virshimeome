from csv import reader
from itertools import groupby
import argparse


def get_matching_ids(dvf_summary_file, min_score_dvf, max_pval_dvf):
    # Get IDs of viral contigs matching filters. 
    viral_contig_ids = set()
    with open(dvf_summary_file, "r") as dvf_file:
        dvf_summary_read = reader(dvf_file, delimiter="\t")
        next(dvf_summary_read) # Skip header
        for row in dvf_summary_read:
            seq_id = row[0]
            score = float(row[2])
            pvalue = float(row[3])
            if score >= min_score_dvf and pvalue <= max_pval_dvf:
                viral_contig_ids.add(">" + seq_id)
    return viral_contig_ids

def combine__seqids_in_file(viral_combined_dvf,viral_contig_ids,fasta_file ):
    
    def read_fasta_file(fasta_file):
        # https://www.biostars.org/p/710/#383479
        with open(fasta_file, 'rb') as fasta_file_read:
            faiter = (x[1] for x in groupby(fasta_file_read, lambda line: str(line, 'utf-8')[0] == ">"))
            for header in faiter:
                seq_id = str(header.__next__(), 'utf-8').replace("\n", "")
                seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__()) + "\n"
                yield seq_id, seq
                
    # Write viral contigs to file.                 
    with open(viral_combined_dvf, "a") as dvf_viral_contig_file:
        for seq_id, seq in read_fasta_file(fasta_file):
            # Iterate through the entire object so the gc can clean up after.                
            if seq_id in viral_contig_ids:
                dvf_viral_contig_file.write(seq_id + "\n" + seq)
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_file', type=str, help='file with all sequence in fasta format. ')
    parser.add_argument('--dvf_summary_file', type=str, help='output of deepvirfinder summary ')
    parser.add_argument('--viral_combined_dvf', type=str, help='fasta file with all viral sequences output of this script')
    parser.add_argument('--min_score_dvf', type=float, help='minimum dvf score ')
    parser.add_argument('--max_pval_dvf', type=float, help='maximum pvalue for dvf hits ')
    
    args = parser.parse_args()
    
    
    fasta_file = args.fasta_file
    dvf_summary_file = args.dvf_summary_file
    min_score_dvf = args.min_score_dvf
    max_pval_dvf = args.max_pval_dvf
    viral_combined_dvf= args.viral_combined_dvf

    viral_contig_ids = get_matching_ids(dvf_summary_file, min_score_dvf, max_pval_dvf)
    combine__seqids_in_file(viral_combined_dvf, viral_contig_ids, fasta_file)

if __name__ == "__main__":
    main()