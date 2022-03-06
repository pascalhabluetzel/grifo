#!/bin/bash

set -e # e=Exit immediately if a command exits with a non-zero exit status
start=$SECONDS

display_usage() {
   echo ""
   echo "
¦¦¦¦¦¦  ¦¦¦¦¦¦  ¦¦  ¦¦¦¦¦  ¦¦¦¦¦
¦¦   ¦¦ ¦¦   ¦¦ ¦¦ ¦¦     ¦¦   ¦¦
¦¦¦¦¦¦  ¦¦¦¦¦¦  ¦¦ ¦¦     ¦¦   ¦¦
¦¦   ¦¦ ¦¦   ¦¦ ¦¦ ¦¦     ¦¦   ¦¦
¦¦¦¦¦¦  ¦¦   ¦¦ ¦¦  ¦¦¦¦¦  ¦¦¦¦¦
   "
   echo "Description: Brico: a tool for taxonomic assignment of nanopore metabarcoding sequences."
   echo "Usage: bash $0 -s my_sequences.fasta"
   echo -e "\t -s \t Path to the sequence file. The file need to be in .fasta format."
   echo -e "\t -t \t Number of threads to be used."
   echo -e "\t -h \t How to use BRICO"
   echo -e "\n"
   echo -e "Example:"
   echo -e "./brico -s my_sequences.fasta"
}

while getopts "hs:t:" arg
    do
     case $arg in
        h)
          display_usage
          exit
          ;;
        s) 
          sample=${OPTARG}
          ;;
        t)
          threads=${OPTARG}
          ;;
     esac
done

shift $((OPTIND-1)) # Removes all the options that have been parsed by getopts from the parameters list

# Check whether input is in the right format
    if ! [[ "${sample}" = *.fasta ]]
    then
        printf "\nERROR: -s: BRICO needs sequences in .fasta format as input."
        display_usage
        exit 1
    fi

# Set default value for the number of threads to 1
    if [ -z "${threads}" ]
    then
        threads=1
    fi

# Load the configuration file
source ./brico.cfg


# If brico has been used before in the same directory, the previous result will be removed
#[[ -f "${sample}".gz ]] && rm "${sample}".gz

gzip "${sample}"

end=$SECONDS
duration=$(( end - start ))
echo "All good! It took BRICO $duration seconds to complete the job."