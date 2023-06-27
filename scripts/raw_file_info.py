import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import argparse
from dvf_to_sequences import read_fasta_file
from statistics import median, mean

def plot_lengths(fasta_file, outfile):
    fasta_gen = read_fasta_file(fasta_file)
    start_position = 1
    sequence_lengths = []

    # Open the FASTA file
    for seq_id, sequence in fasta_gen:
        # Count the length of each sequence
        sequence_lengths.append(len(sequence))

    # Calculate the intervals
    num_intervals = 20
    interval_size = 20000

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(sequence_lengths, bins=num_intervals, color='skyblue', edgecolor='black')

    # Set plot title and labels
    ax.set_title('Sequence Length Distribution')
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Frequency (n)')

    # Set x-axis ticks at 20kb intervals with scientific notation
    x_ticks = range(start_position, max(sequence_lengths) + interval_size, interval_size)
    ax.set_xticks(x_ticks)
    ax.tick_params(axis='x', labelrotation=35)  
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1e'))
    plt.yscale("log")  
    fig.tight_layout()
    
    # Save the plot
    plt.savefig(outfile)

    
def analyze_fasta_file(fasta_file):
    fasta_gen = read_fasta_file(fasta_file)     
    # Initialize variables
    sequence_count = 0
    gc_count = 0
    total_bases = 0
    sequences_length = []
    # Open the FASTA file
    for seq_id, sequence in fasta_gen:
        sequences_length.append(len(sequence))
        
        # Calculate GC content for sequence lines
        gc_count += sequence.count('G') + sequence.count('C')

        
    total_bases = sum(sequences_length)
    # Calculate GC content percentage 

    data_out = {"gc_content": (gc_count / total_bases) * 100,
                "median": median(sequences_length),
                "average": mean(sequences_length),
                }
    # Print the analysis results
    return data_out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_file', type=str, help='fasta_file you want to analyze', default="")
    parser.add_argument('--output_file', type=str, help='data_visualization directory where the figures are dumped.')
    args = parser.parse_args()
    
    fasta_file = args.fasta_file
    output_file = args.output_file
    
    # I opted for two seperate files here because it is easier to have an over. 
    # If needed you can merge them. 
    print(analyze_fasta_file(fasta_file))
    plot_lengths(fasta_file, output_file)

if __name__ == "__main__":
        main()