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

def X_percent(sequences):
    """
    X_percent() will take a list of Biopython sequences and return a NumPy array containing the percentage X content for each sequence.
    """
    #Create an empty list to store the X content of each sequence
    x_content=[]
    #For each sequence:
    for i in range(len(sequences)):
        #Convert Sequence to String
        seq=str(sequences[i].seq)
        #Count X's in sequence
        x_freq=seq.count("X")
        #Compute percent X in the sequence using the length of the sequence
        cent_x=(x_freq/len(seq))*100
        #Add the X content percentage to the list
        x_content.append(cent_x)
    return np.array(x_content)

def X_freq():
    """
    X_Freq() will take a list of Biopython sequences and return a list containing the 
    number of X codons in each sequence.
    
    The X content is returned as a NumPy array.
    """
    #Create an empty list to store the X content of each sequence
    x_content=[]
    #For each sequence:
    for i in range(len(sequences)):
        #Convert Sequence to String
        seq=str(sequences[i].seq)
        #Count X's in sequence
        x_freq=seq.count("X")
        #Compute percent X in the sequence using the length of the sequence
        cent_x=(x_freq/len(seq))*100
        #Add the X content percentage to the list
        x_content.append(cent_x)
    return np.array(x_content)

def Drange(start,end,step):
    """
    Base range() function can only compute on integers. 

    The way around this is to perform a list comprehension producing the desired range in whole numbers, 
    and then dividing the range by a power of ten to produce the range sequence desired. The division power 
    of ten required to produce the desired range sequence depends on the desired step: if the step is between
    0.1 and 1, a division by 10^1 is required. If the step is between 0.01 and 0.1, a division by 10^2 is 
    required, and so on.
    """
    power_raise=0
    #If the step is less than one, calculate the power of ten required to produce a range sequence with the desired step
    while step<1:
        step=step*10
        power_raise+=1

    #List comprehension formula to return the correct range sequence
    return [x/(10**power_raise) for x in range(int(start*10**power_raise),int(end*10**power_raise),int(step))]

def sequence_lengths(sequences):
    """
    Given a list of biopython seqrecord objects, sequence_lengths will calculate the length of each sequence and return a list with the length of each sequence. 
    """
    #Create a variable to store sequence lengths
    lengths=[]
    #Iterate through all of the sequences, calculating the sequence length for each sequence
    for i in range(len(sequences)):
        seq_length=len(sequences[i].seq)
        #Append the sequence length to the list of sequence lengths
        lengths.append(seq_length)
    lengths=np.array(lengths)
    return lengths
            
def sequence_statistics(sequences,protein):
    """
    sequence statistics()
    
    sequences: a list of biopython seq objects read from a sequence file
    protein: name of the protein, which will be printed upon output (string).
    """
    print("Sequence Ambiguity for ",protein,sep="")
    #Compute X Percent of all sequences
    Xprot=X_Percent(sequences)
    #Compute sequence lengths
    lengths=sequence_lengths(sequences)
    #Display table of sequence frequency by percentage X 
    Percentage_Table(0,5.1,.1,list(Xprot))
    print()
    print("Statistics for {}".format(protein))
    print(stats.describe(lengths))
    print(stats.mode(lengths),end="\n\n")
    lower_bound=int(stats.mode(lengths).mode-10)
    print("The sequence length cutoff (reference sequence-10), is {}.".format(lower_bound))
    print("Sequences with length less than {}: {}".format(lower_bound,len(lengths[lengths<(lower_bound)])))
    print("Sequence mode equals reference length: {}".format(int(stats.mode(lengths).mode)==prot_ref_length[protein]))

def percentage_table(data,bin_list,csv_out=None):
    """
    Percentage_Table() will print a table grouping the elements of a list or numpy array into specified bins.
    
    Arguments:
        data: a list or numpy array of values to sort into bins
        bin list: a list consisting of the boundaries for the desired bins. Bins will be drawn between each value listed: for example, [0,0.1,0.2] will draw two bins: the first from 0 to 0.1, and the second from 0.1 to 0.2.
        Bins may have even or uneven spacing. Even bins can be created using the built-in range() function for integer bins. For bin values with decimals, use the Drange() function provided in this module.
        csv_out: optional. A CSV file with the table printed will be created at the specified path.
    """
    #A dictionary is used to count how many values fall into each bin 
    bin_dict={}
    #Dictionary must be set up beforehand with each bin range, all initially mapped to zero (zero values in that range).
    for j in range(len(bin_list)-1):
        #Bins are set up in the format of "<number> to <next_number>"
        bin_dict["{} to {}".format(bin_list[j],bin_list[j+1])]=0
    #Create a special bin for values greater than the values contained in the last bin
    bin_dict["Greater than {}".format(bin_list[-1])]=0

    for i in range(len(data)):
        for j in range(len(bin_list)-1):      
            #Special case or the last bin, placement in the bin includes both the lower bound and the upper bound
            if j==len(bin_list)-2:
                if data[i]>=bin_list[j] and data[i]<=bin_list[j+1]:
                    #Increment the number of items in the last bin by one if the item is inside the range of the bin
                    bin_dict["{} to {}".format(bin_list[j],bin_list[j+1])]+=1
                #If the value is greater than the upper bound in the last bin, store it in the "Greater than" bin
                elif data[i]>bin_list[j+1]:
                    bin_dict["Greater than {}".format(bin_list[-1])]+=1
            #For all other bins, inclusion in the bin is defined as 
            #Example: if bin_list[0]=0.0 and bin_list[1]=0.1, then any value between 0.0 and 0.1, not including 0.1,
            #Will be placed into the first bin
            else:
                if data[i]>=bin_list[j] and data[i]<bin_list[j+1]:
                    #If the value falls into the range of the bin, increment the number of items in that bin by one
                    bin_dict["{} to {}".format(bin_list[j],bin_list[j+1])]+=1
    #Print table
    for item in bin_dict:
        #Print table information to console if no output file is given, otherwise print to the output file.
        if csv_out:
            print("{},{}".format(item,bin_dict[item]),end="\n",file=csv_out)
        else:
            print("[{:<8}: {:>3} ({:>2%})]".format(item,bin_dict[item],bin_dict[item]/len(data)),end="\n")