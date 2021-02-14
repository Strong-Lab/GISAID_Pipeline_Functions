#Biopython package
#The path to the biopython packages was not automatically added to anaconda. To allow import of the package, a manual path must be specified.
import os
#Absolute location: 

#os.environ['CONDA_PREFIX']=/home/wish1832/anaconda3
#Code to append biopython location to path
os.environ['PATH'] = os.environ['PATH']+';'+os.environ['CONDA_PREFIX']+"/pkgs/biopython-1.78-py38h7b6447c_0/lib/python3.8/site-packages"

from Bio import AlignIO as msa_in
from Bio.Align import MultipleSeqAlignment as msa
import Bio.SeqRecord
from Bio import SeqIO
from scipy import stats
import numpy as np
import pandas as pd

#Functions in this module depend on the seq_tools module in same folder. 
from Global_Scripts.Seq_Tools import X_percent,percentage_table
    
def sequence_filter(input_directory,suffix,X_cutoff,max_length_diff,output_directory=None):
    """
    sequence_filter() will filter all of the FASTA files in a specified directory based on the criteria passed to the function. The function will ignore the initial FASTA file dowloaded from GISAID (named 'allprot'). All other files in the directory should be output files from the protein file splitter function.
    
    Arguments
    input_directory: directory of FASTA files for filtering
    suffix: the date of data download, which will be added to the filenames of the output
    X_cutoff: the maximum content of ambiguous "X" codons desired. Should be entered as a fractional value.
    max_length_diff: the maximum deviation in length from the reference sequence for the protein desired. All sequences with a length less than the reference length minus this value, or greater than the reference length plus this value will be excluded.
    output directory: optional. May specify a directory for storing the filtered FASTA files.
    """
    directory=input_directory
    #Count how many files included for progress reporting
    n_files=0
    for filename in os.listdir(directory):
        #Include only FASTA files for each individual protein (exclude file for all proteins)
        if ".fasta" in filename and "allprot" not in filename:
            n_files+=1

    #Sequence filter for each file
    #Create counter for progress reporting
    current_file=0
    for filename in os.listdir(directory):
        if ".fasta" in filename and "allprot" not in filename:
            print("-"*80)
            current_file+=1
            print("Loading file {} of {}".format(current_file,n_files))
            #Retreive protein name
            protein=filename.split("_")[0]

            #For Nsp proteins, filename may be in uppercase. If this is the case, convert to the desired title case.
            if "NSP" in protein:
                #Convert to title case
                protein=protein.title()

            #Load sequences
            sequences=list(SeqIO.parse(directory+filename,"fasta"))
            print("Sequence Information for {}".format(protein))

            #Find the reference sequence WIV04 and determine its length for sequence filtering
            for sequence in sequences:
                if "WIV04" in sequence.description:
                    print("Reference sequence")
                    print(sequence)
                    #Store the length of the reference sequence when found, print the length, and exit the loop
                    ref_length=len(sequence.seq)
                    print("\nLength of Reference Sequence: {}".format(ref_length))
                    print("Length cutoff entered: +-{}".format(max_length_diff))
                    print("Lower bound of sequence length: {}-{}={}".format(ref_length,max_length_diff,ref_length-max_length_diff))
                    print("Upper bound of sequence length: {}+{}={}".format(ref_length,max_length_diff,ref_length+max_length_diff))
                    break

            #Print a Summary of Percentage Sequence Ambiguity
            Xprot=X_percent(sequences)
            print("\nNumber of Sequences by Percent Ambiguity")
            percentage_table(Xprot,bin_list=[0,0.1,1,10])
            print("Ambiguity_Cutoff: {0:.1f} percent".format(X_cutoff*100))

            #Filter Sequences and Write Out
            print("\nSequence filtering for {}".format(protein))
            #Designate output filename
            if output_directory==None:
                outfile=directory+protein+"_"+suffix+"_filtered.fasta"
            else:
                outfile=output_directory+protein+"_"+suffix+"_filtered.fasta"
            #Define lower and upper bounds based on the maximum deviation passed to the function
            lower_bound=ref_length-max_length_diff
            upper_bound=ref_length+max_length_diff
            #individual_filter() will filter the sequences from the current file and write to output
            individual_filter(sequences,outfile,X_cutoff,lower_bound,upper_bound)  
    
def individual_filter(sequences,outfile,cutoff,lower_bound,upper_bound):
    """
    Operates on FASTA genomes, removing any genomes that have a percent ambiguity content greater than 
    a defined threshold, and a length that is not within the lower and upper bounds.
    
    infile: the path of the input file to read genomes from. Must be in FASTA format.
    outfile: path of the output file to write the filtered genomes.
    X_cutoff: a fractional value for the maximum X content desired. (e.g. 0.01 for 1 percent)
    lower_bound: the lower bound for sequence length selection. All sequences with length less than this value will be excluded from the final output file. 
    upper_bound: the upper bound for sequence length selection. All sequences with length greater than or equal to this value will be excluded from the final output file.
    """
    #Protection against overwriting: if the output file already exists, return an error
    if os.path.exists(outfile):
        raise FileExistsError(outfile+" already exists.")
    
    #Print number of Sequences
    n_seq=len(sequences)
    print("Number of sequences: {}".format(len(sequences)))
    
    #Sequences will be stored in an empty list.
    output_sequences=[]
    
    #Create a counter for each exclusion scenario (sequence length above or below cutoff, sequence ambiguity above cutoff)
    too_short=0
    too_long=0
    too_ambiguous=0
    #Updated: counter for sequences from non-human hosts
    non_human=0
    
    #Iterate through the sequences in the given file, testing if the percent ambiguity is below the defined threshold
    for sequence in sequences:
        #First, test if the sequence was taken from a human host. The host information is in field seven (index=6).
        host=sequence.description.split("|")[6]
        if host=="Human":
            #Sequence length: convert the sequence of the SeqRecord object to a string and measure the length of the sequence.
            #Sequences that are within the defined length limits and have a X content below the threshold will be added
            if len(sequence.seq)>=lower_bound and len(sequence)<upper_bound:
                #Conditional for sequence ambiguity
                if (str(sequence.seq).count("X"))/len(sequence.seq)<=cutoff:
                    output_sequences.append(sequence)
                #If the seqeunce has a fractional ambiguity greater than the cutoff entered, exclude the sequence and
                #increase the count of sequences with ambiguity greater than the cutoff by one.
                else:
                    too_ambiguous+=1
            #If the sequence length is less than the range defined by the lower and upper bound, increase the count of 
            #sequences that are too short by one.
            elif len(sequence.seq)<lower_bound:
                too_short+=1
            #If the sequence length is greater than the range defined by the lower and upper bound, increase the count of 
            #sequences that are too long by one.
            elif len(sequence.seq)>=upper_bound:
                too_long+=1
        #If the sequence was not collected from a human, increase the count of non-human sequences by one.
        else:
            non_human+=1
    
    print("{} sequences were removed.".format(n_seq-len(output_sequences)))
    print("Excluded {} sequences taken from non-human hosts.".format(non_human))
    print("Excluded {} sequences of length less than {}".format(too_short,lower_bound))
    print("Excluded {} sequences of length greater than or equal to {}".format(too_long,upper_bound))
    print("Of the sequences within the defined range for length, {} sequences \nwith a percent X content greater than {:.0f} percent were removed.".format(too_ambiguous,cutoff*100))
    
    #Write filtered sequences to output file
    SeqIO.write(output_sequences,outfile,"fasta")
    print("{} sequences written to {}.".format(len(output_sequences),outfile))