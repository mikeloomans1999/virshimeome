import matplotlib.ticker as mticker
import argparse
from os import path
import matplotlib.pyplot as plt
from dvf_to_sequences import read_fasta_file

def plot_lengths(fasta_file, outfile):
    fasta_gen = read_fasta_file(fasta_file)
    start_position = 1
    sequence_lengths = []

    # Open the FASTA file
    for seq_id, sequence in fasta_gen:
        # Count the length of each sequence
        sequence_lengths.append(len(sequence))

    # Calculate the intervals
    num_intervals = 12
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
    
    # Save the plot
    plt.savefig(outfile)

    
def analyze_fasta_file(fasta_file):
    fasta_gen = read_fasta_file(fasta_file)     
    # Initialize variables
    sequence_count = 0
    gc_count = 0
    total_bases = 0

    # Open the FASTA file
    for seq_id, sequence in fasta_gen:
        sequence_count += 1
        # Calculate GC content for sequence lines
        total_bases += len(sequence)
        gc_count += sequence.count('G') + sequence.count('C')

    # Calculate GC content percentage
    gc_content = (gc_count / total_bases) * 100

    # Print the analysis results
    return sequence_count, (round(gc_content, 2), "%")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_file', type=str, help='fasta_file you want to analyze', default="")
    parser.add_argument('--output_dir', type=str, help='data_visualization directory where the figures are dumped.')
    args = parser.parse_args()
    
    fasta_file = args.fasta_file
    output_dir = args.output_dir
    
    # I opted for two seperate files here because it is easier to have an over. 
    # If needed you can merge them. 
    print(analyze_fasta_file(fasta_file))
    plot_lengths(fasta_file, path.join(output_dir, "sequence_length_by_frequency_raw.png"))

if __name__ == "__main__":
        main()