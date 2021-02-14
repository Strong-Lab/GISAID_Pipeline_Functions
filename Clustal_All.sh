#!/bin/bash

#Clustal_all will perform clustal omega on all FASTA files in a given input directory
#The input directory is provided as the first argument passed to the script.
#Directory must be specified with a forward slash at the end for the script to run properly.

INDIR=$1
OUTDIR=$2
SUFFIX=$3

for FILE in $INDIR*
do
	#Slice filename between "/" and "_" to yield protein name
	PROTEIN=${FILE##*/}
	PROTEIN=${PROTEIN%%_*}
	echo "Protein: ${PROTEIN}" 	
	echo "Input file: $FILE"
	#Forming new filename: use the output directory, the protein name, and the suffix (download date)
	OUTPUT="${OUTDIR}${PROTEIN}_${SUFFIX}_msa.fasta"
	echo "Output file: $OUTPUT"
	echo ""
	#Perform Clustal Omega with three iterations
	echo "Performing Clustal Omega on ${PROTEIN}:"
	time clustalo -i $FILE -t Protein --infmt fa -o $OUTPUT --outfmt fa -v --iter=3
	echo "--------------------------------------------------------------------------------"
done
