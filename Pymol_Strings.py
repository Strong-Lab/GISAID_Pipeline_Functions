import pandas as pd
import numpy as np
import os
from scipy import stats

#Biopython path
#os.environ['CONDA_PREFIX']=/home/wish1832/anaconda3
#Code to append biopython location to path
os.environ['PATH'] = os.environ['PATH']+';'+os.environ['CONDA_PREFIX']+"/pkgs/biopython-1.78-py38h7b6447c_0/lib/python3.8/site-packages"

from Bio import SeqIO
import Bio.SeqRecord

def make_pymol_string(subset):
    #Store all of the residue numbers in the subset in a list
    values=[]
    for value in subset["Residue_Number"]:
        values.append(str(value))
    #Collapse list into a string with format value_1+value_1+...value_n
    delimeter="+"
    return delimeter.join(values)

def pymol_from_list(L):
    delimeter="+"
    return delimeter.join(L)

def pymol_strings_10K(path):
    """
    Given the path to the variant finder reports for all proteins, pymol_strings_10K will return the selection strings to color residues in pymol
    based on the number of genomes wtih variants in each position.
    Selection groups will be based on a log-10 scale, with all positions mutated in more than 10,000 genomes in the same group.
    
    Entering strings into pymol: Use the format below.
        select <selection_name>, <optional:selection or chain> and resi <output from this function>
        
        selection name: desired name for selection (recommended to use the protein name and the number of genomes mutated in this selection group)
        selection or chain: when pymol structures contain proteins other than the one for which selections are being defined,
            they should be specified here to avoid selecting residues from other proteins. For more information on defining selections,
            see https://pymolwiki.org/index.php/Selection_Algebra
    
    The correct path should contain one directory for each protein, with each of those directories containing the 'variant_counts.csv' file.
    This file structure will be created by Auto_Parse() in the MSA_Reader module at the specified output directory.
    
    The output directory path specified when calling Auto_Parse() should be the input path for this function.
    """
    for directory in os.listdir(path):
        #Files are stored in separate directories for each protein. Iteration will proceed for every directory in the 1114 proteins folder. 
        #load the csv file for variant counts by position
        filename=path+directory+"/"+directory+"_variant_counts.csv"

        #Print header
        print("-"*80)
        print(filename)
        print(directory,end="\n\n")

        #Load csv file
        Prot=pd.read_csv(filename)

        #Form subset dataframes based on criteria
        Prot_ten_thousand=Prot[Prot["Total_Variants"]>=10000]
        Prot_thousands=Prot[(Prot["Total_Variants"]>=1000) & (Prot["Total_Variants"]<10000)]
        Prot_hundreds=Prot[(Prot["Total_Variants"]>=100)&(Prot["Total_Variants"]<1000)]
        Prot_tens=Prot[(Prot["Total_Variants"]>=10)&(Prot["Total_Variants"]<100)]
        Prot_ones=Prot[(Prot["Total_Variants"]<10) & (Prot["Total_Variants"]>=2)]
        Prot_zero=Prot[Prot["Total_Variants"]==0]

        #Check for errors: ensure all residues are included in the subsets above
        if (len(Prot_zero)+len(Prot_ones)+len(Prot_tens)+len(Prot_hundreds)+len(Prot_thousands)+len(Prot_ten_thousand))==len(Prot):
            #If no errors, print the pymol strings for the protein
            print("Pymol string to select residues with more than 10,000 variants")
            print(make_pymol_string(Prot_ten_thousand))
            print("\nPymol string to select residues with 1000-10000 variants")
            print(make_pymol_string(Prot_thousands))
            print("\nPymol string to select residues with 100-1000 variants")
            print(make_pymol_string(Prot_hundreds))
            print("\nPymol string to select residues with 10-100 variants")
            print(make_pymol_string(Prot_tens))
            print("\nPymol string to select residues with 2-10 variants")
            print(make_pymol_string(Prot_ones))
            print("\nPymol string to select residues with zero variants")
            print(make_pymol_string(Prot_zero))
            print()
        else:
            raise ValueError("Some residues in protein {} do not match any of the criteria provided.".format(directory))