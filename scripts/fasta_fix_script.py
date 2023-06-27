import argparse

def parse_broken_fasta(file_path):
    fasta_dict = {}
    curr_seq = ""
    curr_header = ""

    with open(file_path, "r") as file:
        for line in file:
            line_split =line.split("	")
            fasta_dict[line_split[0]] = line_split[1]

        # Add the last sequence to the dictionary
        if curr_header:
            fasta_dict[curr_header] = curr_seq

    return fasta_dict

def write_fasta(file_path, fasta_dict):
    # Makes sense for our use case with snakemake to use append here. 
    with open(file_path, 'a') as file:
        for header, seq in fasta_dict.items():
            file.write(f">{header}\n{seq}\n")
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='input fasta')
    parser.add_argument('--output', type=str, help='output fasta')
    args = parser.parse_args()

    fasta_file_in = args.input
    fasta_file_out = args.output

    fasta_dict = parse_broken_fasta(fasta_file_in)
    write_fasta(fasta_file_out, fasta_dict)

if __name__ == "__main__":
    main()