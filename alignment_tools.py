#Import Modules and functions
import pandas as pd
import numpy as np
import os,sys

#Import biopython modules
#os.environ['CONDA_PREFIX']=/home/wish1832/anaconda3
#Code to append biopython location to path
os.environ['PATH'] = os.environ['PATH']+';'+os.environ['CONDA_PREFIX']+"/pkgs/biopython-1.78-py38h7b6447c_0/lib/python3.8/site-packages"

from Bio import AlignIO as msa_in
from Bio.Align import MultipleSeqAlignment as msa
import Bio.SeqRecord
from Bio import SeqIO

#Import Reference Cluster Finder
from MSA_Reader import find_ref_cluster_id

def read_seq(path):
    """
    Reads a multiple sequence alignment in FASTA format and stores the sequences in a format readable by the sort_seqs() and manual_correction() functions.
    """
    return list(SeqIO.parse(path,"fasta"))

def sort_seqs(seq_list):
    """
    Takes a sequence alignment produced by read_seq() and sorts in ascending order by the cluster number (number appearing after "Uniq" in the cluster ID).
    """
    return sorted(seq_list,key=lambda x:int(x.id.split(";")[0].split("Uniq")[1]))

def manual_correction(seq_list,msa_start,msa_end,erroneous_alignment,correct_alignment):
    """
    Takes a sequence alignment produced by read_seq() and makes a specified correction to the alignment.
    
    Arguments
    -----------
    seq_list: a list of BioPython SeqRecord objects created by read_seq()
    
    msa_start: location of the beginning of the region to be corrected. Index is based on the residue postion, starting at 1 (i.e. to start at residue 70, enter 70 instead of 69). 
    
    msa_end: location of the end of the region to be corrected. Index is based on the residue postion, starting at 1 (i.e. to end at residue 80, enter 80 instead of 79).
    
    erroneous_alignment: a string giving the motif of the erroneous alignment within the start and end positions specified. Use capital letters, and use blanks exactly as they appear in the alignment. 
    
    correct_alignment: enter the motif of the correct alignment within the start and end positions specified. Use capital letters, and use blanks exactly as they appear in the correct alignment. 
    
    Returns:
    -----------
    seq_list: a list of BioPython SeqRecord objects with the correction made.
    
    Example
    ----------
    To correct the alignment below, enter 72 and 76 as the msa_start and msa_end values, 'A-I--' for erroneous_alignment, and '---AI' for correct_alignment.
    
    72 ---AI 76
       A-I--
       ---AI
       ---AI
    72 ---AI 76
    """
    #Count sequences needing correction, and store information on sequences for display (if number of sequences to correct is below 25)
    n_seq_correct=0
    correction_info=[]
    #Identity Erroneous sequences: print the sequence along with the cluster ID in which it appears
    for i in range(len(seq_list)):
        if seq_list[i].seq[(msa_start-1):msa_end]==erroneous_alignment:
            #Increase count of sequences to correct by one
            n_seq_correct+=1
            cluster_ID=seq_list[i].id.split(";")[0]
            correction_string=f"{cluster_ID} : {msa_start} {seq_list[i].seq[msa_start-1:msa_end]} {msa_end}"
            correction_info.append(correction_string)
    
    #Print information on the sequences to be corrected, depending on the number of erroneous sequences found
    if n_seq_correct<=25:
        #Print all sequences to be corrected if the number of sequences to be corrected is less than or equal to 25
        print(f"Found {n_seq_correct} sequences with the specified motif {erroneous_alignment}.")
        print(*(cluster for cluster in correction_info),sep="\n")
    elif n_seq_correct==0:
        print("No sequences found with the specified erroneous motif.")
    else:
        #If the number of sequences found is greater than 25, print the amount of sequences found and the first 25 elements
        print(f"Found {n_seq_correct} sequences with the specified motif {erroneous_alignment}.")
        print(*(cluster for cluster in correction_info[0:25]),"...",sep="\n")
        
    #Edit the sequence alignment
    for i in range(len(seq_list)):
        if seq_list[i].seq[(msa_start-1):msa_end]==erroneous_alignment:
            #Convert SeqRecord Object to MutableSeq object to make changes
            seq_list[i].seq=seq_list[i].seq.tomutable()
            #Modify the string within the region specified to be equal to the correct alingment
            seq_list[i].seq[msa_start-1:msa_end]=correct_alignment
            #Convert the sequence back into an immutable Seq object
            seq_list[i].seq=seq_list[i].seq.toseq()
    
    print(f"Complete. {n_seq_correct} sequences changed to {correct_alignment}.",end="\n\n")
    return seq_list

def write_sequences(seq_list,path):
    """
    Writes the list of SeqRecord objects to a FASTA file at the defined path.
    """
    SeqIO.write(seq_list,path,"fasta")
    
def sort_FASTA(in_path,out_path):
    """
    Loads a multiple sequence alignment at in_path, sorts sequences in ascending order by their cluster ID, and writes them as a FASTA file to out_path. The sequence headers must be in the format produced by USEARCH for this function to work. 
    """
    sequences=read_seq(in_path,out_path)
    sequences=sort_seqs(sequences)
    write_sequences(sequences,out_path)