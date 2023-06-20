###############
##  imports  ##
###############
import matplotlib.pyplot as plt
import numpy as np
from os import path
from Bio import SeqIO
import argparse


    
#################
##  functions  ##
#################
# Compare quality of circular to linear contigs. 
def quality_circular(checkv_qs_df):
    quality_circ = []
    quality_norm = []
    for row in checkv_qs_df:
        row_list = row.tolist()
        if "circular" in row_list[0]:
            quality_circ.append(row_list[7])
        else:
            quality_norm.append(row_list[7])
    return quality_circ, quality_norm

def calculate_percentages(quality_list, quality_types):
    total_count = len(quality_list)
    percentage_dict = {}
    for value in quality_types:
        percentage_dict[value] = (quality_list.count(value) / total_count) * 100
    return percentage_dict

def quality_by_contig_length(checkv_qs_df):
    quality_by_contig_length_dict = {}
    
    for row in checkv_qs_df:
        row_list = row.tolist()
        quality_value = row_list[7]
        contig_length = int(row_list[1])
        
        if quality_value not in quality_by_contig_length_dict:
            quality_by_contig_length_dict[quality_value] = []
        
        quality_by_contig_length_dict[quality_value].append(contig_length)
        
    return quality_by_contig_length_dict


#################
##  load data  ##
#################
def load_data(checkv_quality_sum_dvf_file, checkv_quality_sum_vs_file, checkv_quality_sum_none_file):
    checkv_qs_df_dvf = np.genfromtxt(checkv_quality_sum_dvf_file, delimiter='\t', dtype=str, skip_header=1)
    checkv_qs_df_vs = np.genfromtxt(checkv_quality_sum_vs_file, delimiter='\t', dtype=str, skip_header=1)
    checkv_qs_df_none = np.genfromtxt(checkv_quality_sum_none_file, delimiter='\t', dtype=str, skip_header=1)
    return checkv_qs_df_dvf, checkv_qs_df_vs, checkv_qs_df_none

def plot_by_contig_length(checkv_qs_df, combined_contig_file, output_file):    
    #####################
    ##  Contig length  ##
    #####################
    contig_lengths = sorted(checkv_qs_df[:,1].astype(int).tolist()) # Get all contig lengths 

    # Plot #
    contig_length_figure, contig_length_axes = plt.subplots()
    contig_length_axes.plot(contig_lengths, range(0, len(contig_lengths)))
    plt.savefig(output_file)


    ########################
    ##  Percentage viral  ##
    ########################
    len_combined_contig_sequences = 0
    combined_contig_sequences = SeqIO.parse(combined_contig_file, 'fasta')
    for record in  SeqIO.parse(combined_contig_file, 'fasta'):
        len_combined_contig_sequences += 1
    print( "Number of contigs: ", len_combined_contig_sequences)
    print("Percentage of contigs that are of viral origin", len(contig_lengths) / len_combined_contig_sequences* 100, "%") 


def plot_quality_by_contig(checkv_qs_df_dvf, checkv_qs_df_vs, checkv_qs_df_none, outdir):
    ##############################
    ##  Quality by contig type  ##
    ##############################
    # DVF
    quality_circ_dvf, quality_norm_dvf = quality_circular(checkv_qs_df=checkv_qs_df_dvf)
    # VS
    quality_circ_vs, quality_norm_vs = quality_circular(checkv_qs_df=checkv_qs_df_vs)
    # None
    quality_circ_none, quality_norm_none = quality_circular(checkv_qs_df=checkv_qs_df_none)

    names = ["dvf circular", "vs2 circular", "raw"] #, "dvf linear", "vs2 linear"]
    quality_types = ['Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete']

    ###### ABSOLUTE ######

    absolute_weight_dict = {}
    for quality in quality_types:
        
        absolute_weight_dict[quality] = np.array([
                quality_circ_dvf.count(quality), 
                quality_circ_vs.count(quality),
                quality_circ_none.count(quality)] 
                )
                #, norm_dvf_quality.count(quality)
                #, norm_vs_quality.count(quality)])

    circular_absolute_figure, circular_absolute_axes = plt.subplots()
    bottom = np.zeros(3)

    for quality, weight_absolute in absolute_weight_dict.items():
        p = circular_absolute_axes.bar(names, weight_absolute, 0.5, label=quality, bottom=bottom)
        bottom += weight_absolute

    plt.grid('on', color = 'gainsboro', axis = 'y', linestyle='--', zorder=3)
     # Shrink current axis by 20%
    box = circular_absolute_axes.get_position()
    circular_absolute_axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    circular_absolute_axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.savefig(path.join(outdir,"circular_quality_absolute_bar.png"))



    ###### PERCENTAGES #######
    # Circular
    quality_circ_dvf_perc = calculate_percentages(quality_circ_dvf, quality_types)
    quality_circ_vs_perc = calculate_percentages(quality_circ_vs, quality_types)
    quality_circ_none_perc = calculate_percentages(quality_circ_none, quality_types)
    # Linear
    # quality_norm_dvf_perc = calculate_percentages(quality_norm_dvf, quality_types)
    # quality_norm_vs_perc = calculate_percentages(quality_norm_vs, quality_types)

    weight_percentages = {}
    for quality in quality_types:
        circ_dvf_quality_percentage =  quality_circ_dvf_perc.get(quality)
        circ_vs_quality_percentage = quality_circ_vs_perc.get(quality)
        circ_none_quality_percentage = quality_circ_none_perc.get(quality)
        # norm_dvf_quality = quality_norm_dvf_perc.get(quality)
        # norm_vs_quality = quality_norm_vs_perc.get(quality)
        
        weight_percentages[quality] = np.array([circ_dvf_quality_percentage, circ_vs_quality_percentage, circ_none_quality_percentage]) #, norm_dvf_quality, norm_vs_quality])

    ## Plotting ##
    circular_percentage_figure, circular_percentage_axes = plt.subplots()
    bottom = np.zeros(3)

    for quality, weight_percentage in weight_percentages.items():
        p = circular_percentage_axes.bar(names, weight_percentage, 0.5, label=quality, bottom=bottom)
        bottom += weight_percentage

    plt.grid(color = 'gainsboro', axis = 'y', linestyle='--', zorder=3)
    
    # Shrink current axis by 20%
    box = circular_percentage_axes.get_position()
    circular_percentage_axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    circular_percentage_axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.savefig(path.join(outdir,"circular_quality_percentage_bar.png"))


def quality_by_length(checkv_qs_df, outfile):
    ################################
    ##  Quality by contig length  ##
    ################################
    quality_types = ['Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete']
    
    quality_by_contig_length_dict = quality_by_contig_length(checkv_qs_df)

    # Get average and stdev per quality type and plot those in barplots. 
    values = [quality_by_contig_length_dict.get(quality_type) for quality_type in quality_types]
    values = values = [value if value is not None else 0 for value in values]
    
    contig_lengths_figure, contig_lengths_axes = plt.subplots(figsize=(10,10))
    contig_lengths_axes.set_ylabel("contig length")
    contig_lengths_axes.set_xlabel("contig quality")
    contig_lengths_axes.boxplot(values)
    contig_lengths_axes.set_xticklabels(quality_types)
    plt.xticks(rotation=35)
    plt.savefig(outfile)

def main():
    
    #################
    ##  variables  ##
    #################

    parser = argparse.ArgumentParser()
    parser.add_argument('--pipeline_output_dir', type=str, help='global pipeline output directory ')
    args = parser.parse_args()
    pipeline_output_dir = args.pipeline_output_dir


    checkv_dvf_dir = path.join(pipeline_output_dir, "2_checkv", "1_1_dvf")
    checkv_vs_dir = path.join(pipeline_output_dir, "2_checkv", "1_1_vs")
    check_v_none_dir = path.join(pipeline_output_dir, "2_checkv", "0_filtered_sequences")
    visualization_output_dir = path.join(pipeline_output_dir, "data_visualization")

    combined_contig_file = path.join(pipeline_output_dir, "0_filtered_sequences", "final-viral-combined.fa")
    checkv_quality_sum_dvf_file = path.join(checkv_dvf_dir, "quality_summary.tsv")
    checkv_quality_sum_vs_file = path.join(checkv_vs_dir, "quality_summary.tsv")
    checkv_quality_sum_none_file =  path.join(check_v_none_dir, "quality_summary.tsv")

    viral_fasta_file_dict = {
        'checkv_dvf': path.join(checkv_dvf_dir, "viruses.fna"),
        'checkv_vs':  path.join(checkv_vs_dir, "viruses.fna"),
        'dvf_prediction': path.join(pipeline_output_dir, "1_1_dvf", "final-viral-combined.fa"),
        'vs_prediction': path.join(pipeline_output_dir, "1_1_vs", "final-viral-combined.fa")
    }
    
    #################
    ##  functions  ##
    #################
    
    # Load data
    checkv_qs_df_dvf, checkv_qs_df_vs, checkv_qs_df_none = load_data(checkv_quality_sum_dvf_file, checkv_quality_sum_vs_file, checkv_quality_sum_none_file)
    
    # Plot by length
    plot_by_contig_length(checkv_qs_df_dvf, combined_contig_file, path.join(visualization_output_dir, "contig_length_by_frequency_dvf.svg"))
    plot_by_contig_length(checkv_qs_df_vs, combined_contig_file, path.join(visualization_output_dir, "contig_length_by_frequency_vs.svg"))
    plot_by_contig_length(checkv_qs_df_none, combined_contig_file, path.join(visualization_output_dir, "contig_length_by_frequency_none.svg"))
    
    # Plot by contig
    plot_quality_by_contig(checkv_qs_df_dvf,checkv_qs_df_vs, checkv_qs_df_none, visualization_output_dir)
    
    # Quality by length
    quality_by_length(checkv_qs_df_dvf, path.join(visualization_output_dir, "dvf_contig_quality_boxplot.png"))
    quality_by_length(checkv_qs_df_vs, path.join(visualization_output_dir, "vs2_contig_quality_boxplot.png"))
    quality_by_length(checkv_qs_df_none, path.join(visualization_output_dir, "none_contig_quality_boxplot.png"))
    
if __name__ == "__main__":
    main()