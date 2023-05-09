# TODO LIST
1. Get list of fasta files from contigs less than 500kb in length and circular. 
2. Use a vir predictor method to get contigs that match. 
3. Calculate relative abundance. 


            # # When you have to read it all. 
            # all_fasta_filenames = glob.glob(path.join(contig_dir,"**", "*.fasta"), recursive=True)
            # for fasta_filename in all_fasta_filenames:
            #     fasta_length = len(get_fasta_sequence(fasta_filename))
            #     if fasta_length <= min_length:
            #         accepted_fasta_absolute_filenames.append(fasta_filename)
