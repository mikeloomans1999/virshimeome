The virshimeome pipeline is a bioinformatics tool that enables the identification and analysis of viral contigs in metagenomic datasets. 
The pipeline starts with the identification of viral contigs using the virsorter2 tool, which classifies contigs based on their similarity to viral genomes.
The resulting viral contigs are then validated using checkV to ensure that they are of high quality and accurately represent viral genomes. 
The next step involves mapping high quality reads to the identified viral contigs using samtools and msamtools, which allows for the calculation of the relative abundance of viral sequences in the metagenomic dataset. 
The output of the pipeline provides valuable insights into the viral diversity and abundance in the studied metagenomic samples, which is essential for understanding the role of viruses in various environmental and clinical contexts.


# TODO LIST
1. Get list of fasta files from contigs less than 500kb in length and circular. 
2. Use a vir predictor method to get contigs that match. 
3. Calculate relative abundance. 
