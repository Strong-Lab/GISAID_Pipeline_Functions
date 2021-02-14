#!/usr/bin/env Rcript

#Command Line Arguments ---------
library(optparse)

option_list=list(
  make_option(c("-i", "--input"),type="character",default=NULL,help="Path to directory containing cluster information files created by USEARCH. Find_Reference_clusers.R will iterate through all cluster information files in the directory.",metavar="character"), 
  make_option(c("-o", "--output"),type="character",default="where_reference.tsv", help="Path and filename of the output file. The file must be in the same directory as the multiple sequence alignment FASTA files and end with 'where_reference.tsv' to be read by the multiple sequence alignment reader.", metavar="character"), 
  make_option(c("-r","--ref_id"),type="character",default=NULL, help="String with any segment of the FASTA header of the reference sequence that is unique to the reference. This can be the full isolate name ('hCoV-19/Wuhan/WIV04/2019'), the abbreviated isolate name ('WIV04'), or the GISAID accession ID ('EPI_ISL_402124').",metavar="character")
  );
#Convert command line input to variables to be passed to functions
#Read arguments and convert to format used by parse_args()
option_parse=OptionParser(option_list=option_list);
#Store arguments according to the long flag (double dashes) in the option list above
args=parse_args(option_parse);

# Functions ----------
#Master function: runs the functions below to record reference cluster IDs.
library(dplyr,quietly = TRUE,warn.conflicts = FALSE)
find_ref_clusters <-function(input_directory,output_path,reference_ID){
  df<-reference_cluster(input_directory,reference_ID,output_path)
  print(df)
}

reference_cluster <- function(directory,reference_ID,output_path){
  #Record number of files in directory
  n_files<-length(list.files(directory,pattern = "clusters.tsv"))
  #Set up object with filenames to iterate through
  files<-list.files(directory,pattern = "clusters.tsv")
  #Set up dataframe for output (2 columns and as many rows as there are files)
  df <- matrix(nrow = n_files,ncol = 2) 
  #Store number of files read for progress meter
  progress<-0
  for (i in 1:n_files){
    #Progress meter
    percent_progress<-(progress/n_files)*100
    cat(paste0("\rFiles Read: ",progress," of ",n_files,".(",round(percent_progress,digits = 1),"% complete)"),sep="\r")
    #Step 2A: Initialize input path and extract protein name
    #Define input path by adding directory path before the filename
    input_path=paste0(directory,files[i])
    #Extract the protein name (before the first underscore in the filename)
    protein<-strsplit(files[i],split="_")[[1]][1]
    
    #Step 2B: load cluster information file and store ID of sequence with the reference genome 
    #Load cluster information file and name columns according to --tabbedout output from USEARCH (see link below) 
    #https://drive5.com/usearch/manual/cmd_fastx_uniques.html
    cluster <- read.delim(input_path,sep="\t",header = FALSE,fill = TRUE)
    colnames(cluster) <- c("Input_ID","Cluster_Name","Cluster_num","Member_num","Cluster_Size","Target_Seq")
    
    #Search the Cluster ID column (cluster) for the FASTA header of the reference strain and store the cluster ID. 
    reference <- cluster %>% 
      filter(grepl(reference_ID,Input_ID)==TRUE) %>% 
      select(Cluster_Name)
    
    #Step 2C: write output based on conditionals
    if (nrow(reference)>1){
      #Conditional 2.C.1: number of rows in reference is greater than one (more than one entry of reference cluster found)
      #In the GISAID FASTA data, the reference cluster is printed twice and grouped in the same cluster. This has a minimal impact on 
      #analysis due to the number of genomes analyzed and the fact that the reference cluster by definition has no variants to be counted twice.
      #The function below verifies that the reference sequences identified are grouped in the same cluster.
      status <-check_equal(reference)
      #In the case that they are not, a warning is returned and "NA" is written out as the cluster ID.
      if (status==FALSE){
        df[i,1]<-protein
        df[i,2]<-"NA"
        print(paste("Warning: reference strand for ",protein," detected in different clusters. NA written to file.",collapse=""))
      }
      #In the case that all cluster IDs are the same, the first element of the cluster ID matrix is written as the cluster ID.
      else{
        df[i,1]<-protein
        df[i,2]<-as.character(reference[1,1])
      }
    #Conditional 2.C.2: Reference sequence is detected in exactly one cluster. Write cluster name to file in this case.
    }else if(nrow(reference)==1){
      df[i,1]<-protein
      df[i,2]<-as.character(reference[1,1])
    #Conditional 2.C.3: Reference sequence is not detected in any cluster. In this case, print a warning and write "NA" as cluster ID.
    }else if (nrow(reference)==0){
      df[i,1]<-protein
      df[i,2]<-"NA"
      print(paste("Warning: reference strand not found for ",protein,". NA written to file.\n",collapse=""))}
  progress<-progress+1
  }
  cat("\r",paste0("Files Read: ",n_files," of ",n_files," (",round(100,digits=1)," % complete)."),sep = "\n")
  #Step 3: Write matrix as a tab-separated file
  write.table(df,file=output_path,sep="\t",quote=F,row.names = F,col.names = F)
  cat(paste0("Output written to ",args$output))
}

check_equal<-function(reference){
  #Code below will catch cases where the reference appears in multiple clusters
  status<-TRUE
  for(i in 1:nrow(reference)){
    paste0("Testing i=",i)
    #Test if all elements are equal to the first
    result<-isTRUE(all.equal(as.character(reference[i,1]),as.character(reference[1,1])))
    if(result==FALSE){
      print("Warning: reference strand detected in two different clusters.")
      status<-FALSE}}
  return(status)
}

#Main----------
cat("Reference Cluster Finder",sep="\n")
cat(paste0("Input directory: ",args$input),sep="\n")
cat(paste0("Reference ID: ",args$ref_id),sep="\n")
cat(paste0("Output path: ",args$output),sep="\n")
cat("",sep = "\n")
find_ref_clusters(input_directory = args$input,reference_ID = args$ref_id,output_path = args$output)
