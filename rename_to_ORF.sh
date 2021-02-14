#!/bin/bash

#As of December 28th, 2020, FASTA files downloaded from GISAID
#Give the protein names of unnamed open reading frames (ORFS) as
#Ns<number> instead of ORF<number>.

#Given the directory containing the FASTA files as input,
#Rename_to_ORF will rename all files with the pattern NS<digit>
#To ORF<digit>, ignoring all files with the pattern NSP<digit>.

DIR=$1
for FILE in $DIR*
do
	#Slice FILE after the last '/' to remove the path
	FILENAME=${FILE##*/}
	#If the filename contains NS* (and not NSP*), rename to ORF*
	case $FILENAME in
		NS[^P]*) 
			#Use built-in sed to change 'NS' in filename to 'ORF'
			NEWFILE=${FILENAME/NS/ORF}
			#Rename using the path of the new file
			mv $FILE $DIR$NEWFILE
			printf "File: $FILE\nRenamed to: $DIR$NEWFILE\n\n";;
	esac
done
echo "Renaming Complete"
