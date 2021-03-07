#Protein File Splitter
#The main function is split_proteins()
#Import libraries
import resource, os, sys, time

def split_proteins(input_path,directory,date):
    """
    Takes the amino acid FASTA file with all sequences downloaded from GISAID (default name "allprot<download_date>.fasta") and sorts sequences into separate files based on protein.
    
    Arguments
    ----------
    input_path (string): path to the "allprot" FASTA file downloaded from GISAID. 
    
    directory (string): relative path to the Raw_FASTA directory. Files for each protein will be created in this directory.
    
    date (string): The download date suffix that appears in the input file. 
    """
    #Before running protein splitter function: check if there are any FASTA files in the directory besides the input 'allprot' file
    #Other FASTA files may exist from previous runs of the protein splitter function in the same directory
    print(f"Input file: {input_path}")
    found_files=False
    for file in os.listdir(directory):
        if ".fasta" in file and directory+file!=input_path:
            #The first time a file is found (when found_files is false, notify user that other files are found.)
            if found_files==False:
                print(f"Existing FASTA files besides {input_path} detected in directory:")
            #Set found_files to true
            found_files=True
            #Print each file found
            print(file)

    #Ask user to confirm deletion of above files
    if found_files==True:
        choice=ask_to_overwrite()  
        #Take action based on input
        if choice=="Y":
            print("Continuing to protein splitter function.")
        elif choice=="n":
            sys.exit()

    #MAIN

    #Determine Number of Sequnces for use with progress meter
    print()
    print("Step 1 of 2")
    print(f"Reading {input_path}")
    with open(input_path,"r") as file:
        #Initialize Counter
        total_sequences=0
        for line in file:
            #Progress Meter
            if total_sequences%5000==0:
                print(f"Sequences Read: {total_sequences}\r",end="")
            if ">" in line:
                #Read the protein name for each header and create a list of unique protein names
                #protein=protein_name(line)
                if protein_name(line) not in protein_names:
                    protein_names.append(protein_name(line))
                #Advance sequence count by one for each header read
                total_sequences+=1
    print(f"Complete. Sequences Read: {total_sequences}")
    print()

    #Define paths of output files and open them
    output_paths=define_all_files(protein_names,directory,date)
    output_files=open_all_files(protein_names,output_paths)

    #Sort Sequences
    print("Step 2 of 2")
    print(f"Sorting Sequences")
    tic=time.time()
    with open(input_path,"r") as file:
        #Initialize Counter
        sequences=0
        #Open all files for writing

        for line in file:
            #Time estimate: update every 100,000 sequences, excluding the first sequence
            if sequences%100_000==0 and sequences!=0:
                #Calculate time elapsed so far
                intermediate=time.time()
                time_elapsed=intermediate-tic
                #Average read speed: number of sequences read per second so far
                sequences_per_sec=sequences/time_elapsed
                #Calculate sequences remaining. Estimated time will be based on seq. remaining and the average reading speed.
                sequences_remaining=total_sequences-sequences
                time_estimate=sequences_remaining/sequences_per_sec
            #Progress meter: update every 5,000 sequences, excluding the first sequence 
            if sequences%5000==0 and sequences>=100_000:
                print(f"Progress: {(sequences/total_sequences)*100:.{4}} % ({sequences}/{total_sequences}). Estimated time remaining: {time_estimate//60:.0f} min {time_estimate%60:.0f} sec \r",end="")
            elif sequences%5000==0 and sequences<100_000:
                print(f"Progress: {(sequences/total_sequences)*100:.{4}} % ({sequences}/{total_sequences}).\r",end="")

            #If the line begins with ">", it represents a FASTA header. Protein identity will be determined from the FASTA header and saved
            if ">" in line:
                #Advance sequence count by one for each header read
                sequences+=1
                #Store variables after reading each header: these will apply to the sequence after each header and are reset with each new header
                #Store protein name from header.
                protein=protein_name(line)
                #Write the header to the target output file (based on protein)
                output_files[protein].write(line)

            #Reading sequences: sequences always come after the header
            elif ">" not in line:
                #Write the sequence to the FASTA file based on the protein in its header
                output_files[protein].write(line)

    #Close all files
    close_all_files(output_files)

    print(f"Progress: {(total_sequences/total_sequences)*100:.{4}} % ({total_sequences}/{total_sequences}). Estimated time remaining: 0 min 0 sec \r",end="")

    toc=time.time()

    print(f"\nComplete. Sequences Sorted: {sequences}")
    print(f"Time elapsed: {toc-tic} sec")
    print(f"Average sorting speed: {sequences/(toc-tic)} sequences/sec")
    print(f"Files created at {os.path.abspath(directory)}")

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

def protein_name(header_line):
    """
    Reads the header line and returns the name of the protein corresponding to the header and the sequence.
    """
    return header_line.split("|")[0].split(">")[1]

def define_all_files(protein_list, directory, date):
    """
    Given a list of protein names, define_all_files will create a dictionary with the path of the output files for each protein.
    """
    output_paths=dict()
    for protein in protein_list:
        output_paths[protein]=directory+protein+"_"+date+".fasta"
    return output_paths 
    
def open_all_files(protein_list,output_paths):
    """
    Taking the list of proteins and the list of output paths generated by define_all_files(), open_all_files() will create a dictionary of file objects with the proteins in protein_list as keys. The files must be closed by close_all_files once all output has been written.
    """
    #Create dictionary for storing file objects
    output_files=dict()
    #Iterate through protein names and open file objects at each protein. Store opened file objects in output_files
    for protein in protein_list:
        output_files[protein]=open(output_paths[protein],"w")
    return output_files

def close_all_files(output_files):
    """
    Closes all file objects created by open_all_files().
    """
    for file in output_files.values():
        file.close
        
def ask_to_overwrite():
    """
    Used at the beginning of the protein_splitter function. Asks the user to confirm overwrite of any FASTA files besides the input file. The function will repeat until a valid option is entered.
    
    Returns
    ----------
    choice: Y or n (whether or not to delete files: Y=yes, n=no, q=quit)
    """
    choice=input(f"Above files may be overwritten. Proceed? Y/n ")
    #Base case: if a valid option is entered, exit function and return choice
    if choice=="Y" or choice=="n":
        return choice
    #Recursive case: if an invalid option is entered, notify user and run funcation again.
    else:
        print("Invalid specification. Please type Y to continuw and n to quit")
        ask_to_overwrite()