#Protein File Splitter
#The main function is split_proteins()
#Import libraries
import resource, os, sys, time

def split_proteins(input_file,directory,date):
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
        print("Duplicate sequences may appear in final output unless these files are deleted.")
        choice=ask_to_delete()  
        #Take action based on input
        if choice=="Y":
            for file in os.listdir(directory):
                if ".fasta" in file and directory+file!=input_path:
                    os.remove(directory+file)
            print("Files Deleted. Continuing to protein splitter function.")
        elif choice=="n":
            print("Files kept. Continuing to protein splitter function.")
        elif choice=="q":
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
            #Advance sequence count by one for each header read
            if ">" in line:
                total_sequences+=1
    print(f"Complete. Sequences Read: {total_sequences}")
    print()

    #Sort Sequences
    print("Step 2 of 2")
    print(f"Sorting Sequences")
    tic=time.time()
    with open(input_path,"r") as file:
        #Initialize Counter
        sequences=0
        for line in file:
            #Time estimate: update every 5,000 sequences, excluding the first sequence
            if sequences%2000==0 and sequences!=0:
                #Calculate time elapsed so far
                intermediate=time.time()
                time_elapsed=intermediate-tic
                #Average read speed: number of sequences read per second so far
                sequences_per_sec=sequences/time_elapsed
                #Calculate sequences remaining. Estimated time will be based on seq. remaining and the average reading speed.
                sequences_remaining=total_sequences-sequences
                time_estimate=sequences_remaining/sequences_per_sec
            #Progress meter: update every 1,000 sequences, excluding the first sequence 
            if sequences%1000==0 and sequences>=2000:
                print(f"Progress: {(sequences/total_sequences)*100:.{4}} % ({sequences}/{total_sequences}). Estimated time remaining: {time_estimate//60:.0f} min {time_estimate%60:.0f} sec \r",end="")
            elif sequences%1000==0 and sequences<2000:
                print(f"Progress: {(sequences/total_sequences)*100:.{4}} % ({sequences}/{total_sequences}).\r",end="")

            #If the line begins with ">", it represents a FASTA header. Protein identity will be determined from the FASTA header and saved
            if ">" in line:
                #Advance sequence count by one for each header read
                sequences+=1
                #Store variables after reading each header: these will apply to the sequence after each header and are reset with each new header
                #Store protein name from header.
                protein=protein_name(line)
                #Set the output file based on the protein name
                output_path=directory+protein+"_"+date+".fasta"
                #Write the header to the target output file
                writeline_FASTA(output_path,line.strip())

            #Reading sequences: sequences always come after the header
            elif ">" not in line:
                #Write the sequence to the FASTA file based on the protein in its header
                writeline_FASTA(output_path,line.strip())

    print(f"Progress: {(total_sequences/total_sequences)*100:.{4}} % ({total_sequences}/{total_sequences}). Estimated time remaining: 0 min 0 sec \r",end="")
    toc=time.time()

    print(f"\nComplete. Sequences Sorted: {sequences}")
    print(f"Time elapsed: {toc-tic} sec")
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

#Function to write lines to FASTA file, appending to the end
def writeline_FASTA(name, line):
    with open(name, 'a') as file:
        file.write(line+"\n")
        
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

def ask_to_delete():
    """
    Used at the beginning of the protein_splitter function. Asks the user to confirm deletion of any FASTA files besides the input file. The function will repeat until a valid option is entered.
    
    Returns
    ----------
    choice: Y or n (whether or not to delete files: Y=yes, n=no, q=quit)
    """
    choice=input(f"Delete above files? Y/n/q ")
    #Base case: if a valid option is entered, exit function and return choice
    if choice=="Y" or choice=="n" or choice=="q":
        return choice
    #Recursive case: if an invalid option is entered, notify user and run funcation again.
    else:
        print("Invalid specification. Please type Y to delete files, n to keep files, or q to quit.")
        ask_to_delete()