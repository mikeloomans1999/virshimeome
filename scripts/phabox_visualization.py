from Bio import Phylo
import argparse
from os import path
import matplotlib.pyplot as plt
from csv import reader, DictReader
from random import randint

# function calling 
# "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/final-viral-combined.fa", "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/split", "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/1_1_dvf/split/hello.fa"


# Create a circular taxonomic tree. 
# Species discovered and virulence. 
# Perform accurate blast contigs as the family level blat is based on a limited local database. 
# 

# Virulence
phatyp_prediction = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/8_0_identification/0_filtered_sequences/output/phatyp_prediction.csv"
cherry_host = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/8_0_identification/0_filtered_sequences/output/cherry_prediction.csv"
phagcn_taxonomy = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/8_0_identification/0_filtered_sequences/output/phagcnx_prediction.csv"

outdir = "/projects/arumugam/scratch/mnc390/virome_testing/virshimeome_pipeline_output/data_visualization/"

def get_phage_lifestyle_number(phatyp_prediction_file):
    temperate = 0 
    virulent = 0 
    undefined = 0 
    lifestyle_dict = {} # dict
    with open(phatyp_prediction_file, "r") as phatype_file:
        phatype_file_read = reader(phatype_file, delimiter=",")
        next(phatype_file_read)
        for row in phatype_file_read:
            seq_id, prediction, score = row
            if float(score) > 0.5:
                lifestyle_dict[seq_id] = prediction
                if prediction == "virulent":
                    virulent += 1
                else:
                    temperate += 1
            else:
                lifestyle_dict[seq_id] = "undefined"
                undefined += 1 
    return temperate, virulent, undefined


def write_pie_lifestyle_plot(temperate, virulent, undefined, output_file):
    pie_lifestyle_fig, pie_lifestyle_ax = plt.subplots()

    if undefined > 0:
        pie_lifestyle_ax.pie([virulent, temperate, undefined], 
                            labels=["virulent", "temperate", "undefined"], 
                            autopct=lambda p : '{:.2f}%  ({:,.0f})'.format(p,p * sum([virulent, temperate, undefined])/100)
                            )
    else:
        pie_lifestyle_ax.pie([virulent, temperate], 
                            labels=["virulent", "temperate"], 
                            autopct=lambda p : '{:.2f}%  ({:,.0f})'.format(p,p * sum([virulent, temperate])/100)
                            )

    plt.savefig(output_file)

def get_feature_frequency_dict(prediction_file):
    feature_frequency_dict = {}
    
    with open(prediction_file, 'r') as file:
        reader = DictReader(file)
        for row in reader:
            # Get the feature from the "Pred" column
            feature = row["Pred"]

            # If the feature is already in the dictionary, increment its frequency
            if feature in feature_frequency_dict:
                feature_frequency_dict[feature] += 1
            # If the feature is not in the dictionary, add it with a frequency of 1
            else:
                feature_frequency_dict[feature] = 1

    return feature_frequency_dict

def calculate_other(feature_names, frequencies):
    total = sum(frequencies)
    other_frequency = 0
    pop_n = 0 
    other_frequencies = frequencies.copy()
    freq_threshold = total * 0.02
    print(freq_threshold)
    for index, frequency in enumerate(frequencies):
        if frequency < freq_threshold:
            feature_names.pop(index - pop_n)
            other_frequencies.pop(index - pop_n)
            pop_n += 1
            other_frequency += frequency
    return feature_names, other_frequencies, other_frequency
    
def write_pie_multi(dict_feature_frequency,output_file):
    # Get the feature names and frequencies from the dictionary
    feature_names = list(dict_feature_frequency.keys())
    frequencies = list(dict_feature_frequency.values())

    # Merge uncommon features into "other"
    feature_names, frequencies, other_frequency = calculate_other(feature_names, frequencies)
    if other_frequency > 0:
        feature_names.append("other")
        frequencies.append(other_frequency)
    
    feature_names_absolute_values_comment = [feature_name + f"({frequencies[index]})" for index, feature_name in enumerate(feature_names)]
    
    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(10, 10))
    
    colors = ['#%06X' % randint(0, 0xFFFFFF) for _ in range(len(dict_feature_frequency))]
    # Generate the pie chart
    ax.pie(frequencies, labels=feature_names_absolute_values_comment, autopct='%1.1f%%',  colors=colors,  labeldistance=1.1, pctdistance=0.85)
    fig.tight_layout(pad=5)
    # Set the aspect ratio to be equal so that pie is drawn as a circle
    ax.axis('equal')

    # Save the pie chart to the output file
    plt.savefig(output_file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phatype_prediction_file', type=str, help='phatyp_prediction_file', default="")
    parser.add_argument('--host_prediction_file', type=str, help='host_prediction_file')
    parser.add_argument('--output_file', type=str, help='data_vizualization_dir')
    args = parser.parse_args()
    
    output_file = args.output_file
    host_prediction_file = args.host_prediction_file
    
    phatype_prediction_file = args.phatype_prediction_file
    if phatype_prediction_file != "":
        temperate, virulent, undefined = get_phage_lifestyle_number(phatype_prediction_file)
        write_pie_lifestyle_plot(temperate, virulent, undefined, output_file)
    if host_prediction_file != "":
        dict_feature_frequency = get_feature_frequency_dict(host_prediction_file)
        write_pie_multi(dict_feature_frequency, output_file)
if __name__ == "__main__":
        main()