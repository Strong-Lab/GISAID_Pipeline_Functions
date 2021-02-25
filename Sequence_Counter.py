#The path to the biopython packages was not automatically added to anaconda. To allow import of the package, a manual path must be specified.
import os,sys
#Absolute location: 

#os.environ['CONDA_PREFIX']=/home/wish1832/anaconda3
#Code to append biopython location to path
os.environ['PATH'] = os.environ['PATH']+';'+os.environ['CONDA_PREFIX']+"/pkgs/biopython-1.78-py38h7b6447c_0/lib/python3.8/site-packages"

from Bio import AlignIO as msa_in
from Bio.Align import MultipleSeqAlignment as msa
import Bio.SeqRecord
from Bio import SeqIO
import numpy as np
import pandas as pd
import warnings

def pipeline_seq_counts(main_directory,proteins="get-all",remove_start=0,remove_end=0):
    """
    Returns a dataframe giving the number of sequences included at each of the following steps in the GISAID analysis pipeline
    
    1. Initial Sequences
    2. Sequences Passing Filter Criteria
    3. Sequences Clustered and Aligned
    4. Clustered Sequences Linked to Metadata
    5. Metadata-Sequence Pairs Included in Time Series Analysis
    
    Arguments
    -----------
    main_directory (string): path to the directory with the 'batch' of proteins (batch is a single download from GISAID). The path should point to the master directory which contains the subdirectories Raw_FASTA, Filtered_Sequences, Seq_Derep, etc.

    proteins (list): a list giving the protein names for which to generate sequence counts. May alternately be the string "get-all" to compute counts for all proteins in the directory.
    
    remove_start (integer, default=0):  a positive integer or zero giving the number of weeks to exclude from the start of the time series.
    
    remove_end (integer, default=0): a positive integer or zero giving the number of weeks to exclude from the end of the time series.
    
    Returns
    -----------
    counts_df: a DataFrame giving the number of sequences included at each step in the analysis pipeline.
    """
    #Create dictionary for storing lists of sequence counts for each step
    counts_by_step=dict()
    
    #If proteins is set to get-all, build a list of all proteins processed through all five steps.    
    if proteins=="get-all":
        #Search through time_series directory (proteins that were processed through all five steps will be in this directory)
        proteins=[]
        for file in os.listdir(main_directory+"Time_Series_Tables/Time_Series_Percent/"):
            #Only read filenames for the global seuquence counts to avoid repeats of the same protein
            if "Global" in file:
                #File format: <protein_name>_Percent_Global.csv (protein name is the value to the left of the first underscore)
                protein_name=file.split("_")[0]
                #In case of multiple files, check to see if the protein name already exists in the list before adding it
                if protein_name not in proteins:
                    proteins.append(protein_name)
                    
    #Sequence counts for step 1 (# of initial sequences in each file)
    step="Sequences Downloaded"
    directory="{}Raw_FASTA/".format(main_directory)
    #Progress meter
    progress=1
    print("Step {0} of {1}: {2}".format(progress,5,step))
    sequence_counts=count_FASTA(directory,proteins)
    progress+=1
    #Add list of counts to dictionary under the step
    counts_by_step[step]=sequence_counts

    #Sequence counts for step 2 (# of sequences remaining after filtering)
    step="Sequences Passing Filter Criteria"
    directory="{}Filtered_Sequences/".format(main_directory)
    print("Step {0} of {1}: {2}".format(progress,5,step))
    sequence_counts=count_FASTA(directory,proteins)
    #Add list of counts to dictionary under the step
    counts_by_step[step]=sequence_counts
    #Progress meter: advance by one
    progress+=1

    #Step 3 (# of sequences remaining after clustering and excluding singletons)
    step="Sequences Clustered and Aligned"
    directory="{}Seq_Derep/Cluster_FASTA/".format(main_directory)
    print("Step {0} of {1}: {2}".format(progress,5,step))
    sequence_counts=count_clustered_sequences(directory,proteins)
    #Add list of counts to dictionary under the step
    counts_by_step[step]=sequence_counts
    #Progress meter: advance by one
    progress+=1

    #Step 4 (# of sequences clustered and linked to metadata)
    step="Clustered Sequences Linked to Metadata"
    directory="{}Time_Series_Tables/Metadata_with_Variants/".format(main_directory)
    print("Step {0} of {1}: {2}".format(progress,5,step))
    sequence_counts=seq_in_meta(directory,proteins)
    #Add list of counts to dictionary under the step
    counts_by_step[step]=sequence_counts
    #Progress meter: advance by one
    progress+=1

    #Step 5 (# of sequences in the final time series analysis)
    step="Metadata-Sequence Pairs Included in Time Series Analysis"
    directory="{}Time_Series_Tables/Time_Series_Frequency/".format(main_directory)
    print("Step {0} of {1}: {2}".format(progress,5,step))
    sequence_counts=count_TS(directory,proteins,remove_end=2)
    #Add list of counts to dictionary under the step
    counts_by_step[step]=sequence_counts
    #Progress meter: advance by one
    progress+=1

    #Create a Pandas DataFrame from the results
    #Form a DataFrame, then transpose (each dictionary key represents a row, and default is for each key to represent a column)
    counts_df=pd.DataFrame(counts_by_step).transpose()
    counts_df.columns=proteins
    
    return counts_df

def count_FASTA(directory,proteins):
    """
    Returns the number of sequences present in the FASTA files of the defined proteins, in the target directory.
    """
    #Use sequence_counts to store the number of sequences for each protein
    sequence_counts=[]

    for i in range(len(proteins)):
        print("Counting: ({0} of {1})\r".format(i+1,len(proteins)),end="")
        #Load sequences and count the number of sequences present
        sequences=load_FASTA(directory,proteins[i])
        n_seq=len(sequences)
        #Add the number of sequences to the list
        sequence_counts.append(n_seq)
    print("Counting Complete.")
    
    return sequence_counts

def count_clustered_sequences(directory,proteins):
    """
    Returns a list of the number of sequences in the clustered FASTA file for each protein entered.
    
    arguments
    ----------
    directory: path to the directory containing the clustered FASTA files.
    
    proteins (list): a list of the desired proteins for which to obtain sequence counts
    
    returns
    ----------
    sequence_counts: a list of sequence counts for each protein, in the same order as the protein list passed.
    """
    #Use sequence_counts to store the number of sequences for each protein
    sequence_counts=[]

    #Load FASTA file with clusters and examine headers (headers contain sample size)
    for i in range(len(proteins)):
        #Update progress meter
        print("Counting: ({0} of {1})\r".format(i+1,len(proteins)),end="")
        #Search for the file matching the protein name in the directory
        sequences=load_FASTA(directory,proteins[i])
        #Total number of sequences is equal to the sum of all cluster sizes. Create a counter for storing the sum
        total_sequences=0
        #Get the cluster size for each cluster in the file
        for j in range(len(sequences)):
            #Fetch the header data stored in the attribute .id, and extract the cluster size (between the two semicolons)
            cluster_size=sequences[j].id.split(";")[1]
            #returned string is size=<size>. Extract the number after the equals sign
            cluster_size=cluster_size.split("=")[1]
            #Add cluster size to the total
            total_sequences+=int(cluster_size)
        sequence_counts.append(total_sequences)
    
    print("Counting Complete.")
    return sequence_counts

def seq_in_meta(directory,proteins):
    """
    Returns a list of the number of sequences in the metadata with variants file for each protein entered.
    
    arguments
    ----------   
    directory (string): path to the directory containing the clustered FASTA files.
    
    proteins (list): a list of the desired proteins for which to obtain sequence counts
    
    returns
    ----------
    sequence_counts: a list of sequence counts for each protein, in the same order as the protein list passed.
    """
    #Use sequence_counts to store the number of sequences for each protein
    sequence_counts=[]

    for i in range(len(proteins)):
        #Update progress meter
        print("Counting: ({0} of {1})\r".format(i+1,len(proteins)),end="")
        #Find the path to the metadata with variants file for the protein 
        path=find_path(directory,proteins[i])
        #The number of rows in the metadata file is the number of sequences included
        n_seq=len(pd.read_csv(path,dtype="O",sep="\t"))
        #Add the sequence count for the protein to the list of counts
        sequence_counts.append(n_seq)

    print("Counting Complete.")
    return sequence_counts

def count_TS(directory,proteins,remove_start=0,remove_end=0):
    """
    Returns a list of the number of sequences in the final time series data for each protein entered.
    
    arguments
    ----------   
    directory: path to the directory containing the clustered FASTA files.
    
    proteins (list): a list of the desired proteins for which to obtain sequence counts
    
    remove_start (int): a positive integer or zero giving the number of weeks to exclude from the start of the time series.
    
    remove_end (int): a positive integer or zero giving the number of weeks to exclude from the end of the time series.
    
    returns
    ----------
    sequence_counts: a list of sequence counts for each protein, in the same order as the protein list passed.
    """
    #Use sequence_counts to store the number of sequences for each protein
    sequence_counts=[]

    for i in range(len(proteins)):
        #Update progress meter
        print("Counting: ({0} of {1})\r".format(i+1,len(proteins)),end="")
        #Sequences in Time Series Analysis: sum of 'total genomes' row of global TS frequency data
        path=find_TS_path(directory, proteins[i])
        import Plotting_Functions as plotf
        #Return an error if negative values are entered for remove_start or remove_end 
        if remove_start<0 or remove_end<0:
            raise ValueError("Negative value for remove_start and/or remove_end.")
        #Load time series, set code column as the index, and sum the values in the first row 
        if remove_end==0:
            n_seq=plotf.prepare_TS(path,reindex=True).iloc[0,remove_start:].sum()
        #If a nonzero number of weeks is specified to be removed from the end, slice the number of rows specified from the end with the "-" operator
        elif remove_end>0:
            n_seq=plotf.prepare_TS(path,reindex=True).iloc[0,remove_start:-remove_end].sum()
        #Add the sequence count for the protein to the list of counts
        sequence_counts.append(n_seq)
    
    print("Counting Complete.")
    return sequence_counts

def find_path(directory,protein):
    """
    Searches through the specified target directory for the file matching the given protein, and returns the path of the file.
    Arguments
    ----------
    directory (string): path to the desired directory.

    protein (string): name of the protein to search for in the target directory.

    Returns
    ----------
    path (string): path to the file in the target directory related to the protein
    """
    #Count files found to catch cases where multiple or zero files are found for a protein
    files_found=0
    #Iterate through all files in the directory
    for filename in os.listdir(directory):
        #Protein name may be in uppercase or title case depending the script that created the file
        if (protein.upper() in filename) or (protein in filename):
            files_found+=1
            path=directory+filename

    #Warn if zero files or multiple files are found for a protein
    if files_found==0:
        warnings.warn("File not found for protein {}.".format(protein))
    elif files_found>1:
        warnings.warn("Multiple files found for protein {}. Using the path that comes last alphabetically.".format(protein))
   
    return path

def load_FASTA(directory,protein):
    """
    Searches through the specified target directory for the FASTA file matching the given protein, and returns the sequences from that file.

    Arguments
    ----------
    protein (string): name of the protein to search for in the target directory.

    directory (string): path to the desired directory.

    Returns
    ----------
    sequences: a list of Biopython SeqRecord objects read from the discovered file.
    """
    #Count files found to catch cases where multiple or zero files are found for a protein
    files_found=0
    #Iterate through all files in the directory
    for filename in os.listdir(directory):
        #Protein name may be in uppercase or title case depending the script that created the file
        if (protein.upper() in filename) or (protein in filename):
            files_found+=1
            seq_path=directory+filename

    #Warn if zero files or multiple files are found for a protein
    if files_found==0:
        warnings.warn("File not found for protein {}.".format(protein))
    elif files_found>1:
        warnings.warn("Multiple files found for protein {}. Using the path that comes last alphabetically.".format(protein))
    
    #Load sequences and store as Biopython SeqRecord objects
    sequences=list(SeqIO.parse(seq_path,"fasta"))

    return sequences

def find_TS_path(directory,protein):
    """
    Searches through the specified target directory for the file containing global time series frequency data for the given protein, and returns the path of the file.
   
    Arguments
    ----------
    directory (string): path to the desired directory.
    
    protein (string): name of the protein to search for in the target directory.

    Returns
    ----------
    path (string): path to the file in the target directory related to the protein
    """
    #Count files found to catch cases where multiple or zero files are found for a protein
    files_found=0
    #Iterate through all files in the directory
    for filename in os.listdir(directory):
        #Protein name may be in uppercase or title case depending the script that created the file
        if ((protein.upper() in filename) or (protein in filename))&("Global" in filename):
            files_found+=1
            path=directory+filename

    #Warn if zero files or multiple files are found for a protein
    if files_found==0:
        warnings.warn("File not found for protein {}.".format(protein))
    elif files_found>1:
        warnings.warn("Multiple files found for protein {}. Using the path that comes last alphabetically.".format(protein))
   
    return path