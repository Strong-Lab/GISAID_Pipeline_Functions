#Protein File Splitter
#The main function is split_proteins()

#Import libraries
import os

#Function to open FASTA file
def load_file(name):
    with open(name) as file:
        #Read contents and store in a list with one element for each line
        lines=file.readlines()
        #Remove newline characters from input
        for i in range(len(lines)):
            lines[i]=lines[i].strip()
        #Return the contents of the file as a list
        return lines

#Function to write lines to FASTA file, appending to the end if the file exists
def writeline_FASTA(name, line):
    with open(name, 'a') as file:
        file.write(line+"\n")

#The origional file contains amino acid sequences for all proteins together from every strain.
#The fuction below will separate the data such that there is one file for each protein.
def split_proteins(input_file,directory,date):
    """
    split_proteins() will split the sequences in the FASTA downloaded from GISAID into separate files based on the name of the protein in the metadata entry for each sequence.
    
    arguments
    -----------
    input_file: the raw FASTA file downloaded from GISAID. Name begins with "allprot".
    
    directory: the output directory for writing the individual protein files.
    
    date: the date of the download (used for naming the files). Can be in any format that creates a legal filename.
    """
    print("Opening FASTA file. This may take several minutes.")
    #Load contents of FASTA files
    contents=load_file(input_file)
    #The total number of sequences equals the length of the data divided by 2
    print("Number of sequences read: ",len(contents)/2,sep="",end="\n\n")
    #Iterating through every even index will read header data
    print("Sorting FASTA entries into files by protein.")
    for i in range(0,len(contents),2):
        if i%5000==0 and i!=0:
            #Every 5,000 headers, a progress report will be printed detailing the number of headers sorted out of all headers
            #Half of the data consists of headers (hence the division by 2).
            fraction_complete=(i/2)/(len(contents)/2)
            print("\rProgress: {0}/{1} ({2:.2%})".format(i/2,len(contents)/2,fraction_complete),end="")
        #Read name of protein associated with FASTA header and store name for file creation
        name=protein_name(contents[i])
        #Protein name in header data will determine which file to sort the header and sequence into
        filename=directory+name+"_"+date+".fasta"
        #Append header to FASTA
        writeline_FASTA(filename,contents[i])
        #Append AA sequence to FASTA (sequence is on the line after the header)
        writeline_FASTA(filename,contents[i+1])
    print("\rProgress: {0}/{1} ({2:.2%})".format(len(contents)/2,len(contents)/2,1.00))
    print("Complete. Files created at {}".format(os.path.abspath(directory)))
	
#The function below reads the header line and retrieves the protein name corresponding to the header and the sequence
def protein_name(header):
    #The protein name is listed between the ">" and "|" characters in the FASTA file.
    #The for loop below will iterate through the header until the ">" and "|" characters are found
    for j in range(0,len(header)):
        #When the ">" is found, the index of it is noted. 
        if header[j] == ">":
            #The substring will begin one character after the ">"
            sub_start = j+1
            continue
        #The substring should slice the header right before the "|" character to yield the protein name.
        elif header[j] == "|":
            #Slicing at the index of the "|" character will achieve this.
            sub_end = j
            break
    substring=header[sub_start:sub_end]
    return substring