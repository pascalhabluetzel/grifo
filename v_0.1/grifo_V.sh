#!/bin/bash

#*********************************************************************************#
#                                                                                 #
# GRIFO is a workflow to process ONT data from begin to end. It takes fast5 input #
# and produces OTU abundance tables                                               #
#                                                                                 #
# Developers and maintainers: Pascal Habluetzel and Els De Keyzer                 #
# e-mail: pascal.hablutzel@vliz.be                                                #
#         els.de.keyzer@uantwerpen.be                                             #
#                                                                                 #
#*********************************************************************************#

#***********************************#
#                                   #
# Instructions for Valya            #
#                                   #
#***********************************#

# Install scikit-learn, spoa, hbdscan, pandas, vsearch, minimap2, qcat and nanofilt in your environment.
# Place data in a folder called ./data/ within your working directory.
# Place the following scripts in your working directory and make them executable: ashure.py, bilge_pype.py, grifo_V.sh (the first two can be found here: https://github.com/BBaloglu/ASHURE/tree/master/src)
# Place the grifo.cfg configuration file in your working directory and adjust the minimum and maximum length according to your amplicon length.
# Also in the configuration file, set the number of iterations for the clustering algorithm to a useful value (e.g. 10). More is better but takes more time of course.
# Run the script as ./grifo_V.sh
# In the results folder, there will be some files called barcodeXY_nonchimeras.fasta (these are your ASVs) and barcodeXY_table.csv (these are your read counts per ASV).




#set -e # e=Exit immediately if a command exits with a non-zero exit status; CREST4 produces a non-zero exit. Command can therefore not be used.

start=$SECONDS

# Safe working and results directory
wkdir=$(pwd)
resultsdir="results-"$(date +%F-%R)""

display_usage() {
   echo ""
   echo "
 ¦¦¦¦¦  ¦¦¦¦¦¦  ¦¦ ¦¦¦¦¦¦  ¦¦¦¦¦
¦¦      ¦¦   ¦¦ ¦¦ ¦¦     ¦¦   ¦¦
¦¦  ¦¦¦ ¦¦¦¦¦¦  ¦¦ ¦¦¦¦¦  ¦¦   ¦¦
¦¦   ¦¦ ¦¦   ¦¦ ¦¦ ¦¦     ¦¦   ¦¦
 ¦¦¦¦¦  ¦¦   ¦¦ ¦¦ ¦¦      ¦¦¦¦¦
   "
   echo "Description: GRIFO: a workflow for taxonomic assignment of nanopore metabarcoding sequences."
   echo "Usage: ./grifo.sh"
   echo -e "\t -t \t Number of threads to be used."
   echo -e "\t -h \t How to use GRIFO"
   echo -e "\n"
   echo -e "Example:"
   echo -e "./grifo.sh"
}

while getopts "ht:" arg
    do
     case $arg in
        h)
          display_usage
          exit
          ;;
        t)
          threads=${OPTARG}
          ;;
     esac
done

shift $((OPTIND-1)) # Removes all the options that have been parsed by getopts from the parameters list

# Set default value for the number of threads to 1
    if [ -z "${threads}" ]
    then
        threads=1
    fi

# Load the configuration file
source ./*.cfg


#***********************************#
#                                   #
# Check the input data type         #
# Do basecalling and demultiplexing #
# if not already done               #
#                                   #
#***********************************#

# All data files must be in a folder called ./data/ in the working directory.

if [ $demultiplexed = yes ]
then
      echo "Sequences are already demultiplexed."
elif [ $demultiplexed != yes ]
then
      echo "Check input file type."
# Count files in fastq.gz or fast5 format and provide the names of the subfolders where these files are stored.
      fastq_count=$(find ./data/ -name "*.fastq.gz" | wc -l)
      fastq_dir=$(find ./data/ -type f -name "*.fastq.gz" | sed -r 's|/[^/]+$||' |sort |uniq)
      fast5_count=$(find ./data/ -name "*.fast5" | wc -l)
      fast5_dir=$(find ./data/ -type f -name "*.fast5" | sed -r 's|/[^/]+$||' |sort |uniq)
elif [ $fastq_count -gt 0 ]
then
      echo "found $fastq_count g-zipped fastq files in $fastq_dir, initiating filtering"
elif [ $fast5_count -gt 0 ]
then
      echo  "found $fast5_count fast5 files in $fast5_dir, initiating base calling"
      for val in $fast5_dir
      do
      mkdir -p ./$resultsdir/basecalling/$val
      /opt/ont-guppy-cpu_3.0.3/bin/guppy_basecaller -i .$val/ -s ./$resultsdir/basecalling/$val --flowcell $flow_cell --kit $library_kit
      done
else
      echo "\nERROR: GRIFO needs sequences in .fast5 or .fastq format as input."
fi


#***********************************#
#                                   #
# Filtering                         #
#                                   #
#***********************************#

echo "Start filtering..."
fastq_dir=$(find ./data/ -type f -name "*.fastq.gz" | sed -r 's|/[^/]+$||' |sort |uniq)
echo $fastq_dir
for dir in $fastq_dir
      do
      echo "Filtering ${dir##*/}..."
      mkdir -p ./$resultsdir/qc/"${dir##*/}"
      cd $dir
        for fastq in *.fastq.gz
            do
            echo $(pwd)/"$fastq"
            gunzip -c "$fastq" | NanoFilt --length $minimum_length --maxlength $maximum_length -q $qscore | gzip > $wkdir/$resultsdir/qc/"${dir##*/}"/"$fastq"
        done
      cat $wkdir/$resultsdir/qc/"${dir##*/}"/*.fastq.gz > $wkdir/$resultsdir/qc/"${dir##*/}"/"${dir##*/}"_concatenated.fastq.gz
      gunzip $wkdir/$resultsdir/qc/"${dir##*/}"/"${dir##*/}"_concatenated.fastq.gz
      sed -n '1~4s/^@/>/p;2~4p' $wkdir/$resultsdir/qc/"${dir##*/}"/"${dir##*/}"_concatenated.fastq > $wkdir/"${dir##*/}"_concatenated.fasta
      gzip $wkdir/$resultsdir/qc/"${dir##*/}"/"${dir##*/}"_concatenated.fastq
      cd $wkdir
done


#***********************************#
#                                   #
# Convert fasta to csv              #
#                                   #
#***********************************#

for dir in $fastq_dir
do
mv $wkdir/"${dir##*/}"_concatenated.fasta input.fasta
fasta_to_csv () {
    PYCMD=$(cat <<EOF
fasta = 'input.fasta'
output = 'output.csv'

out_lines = []
temp_line = ''
with open(fasta, 'r') as fp:
    for line in fp:
        if line.startswith('>'):
            out_lines.append(temp_line)
            temp_line = line.strip() + ','
        else:
            temp_line += line.strip()
out_lines.append(temp_line)

with open(output, 'w') as fp_out:
    fp_out.write('id,sequence' + '\n'.join(out_lines))
EOF
    )

    python3 -c "$PYCMD"
}
fasta_to_csv
mv output.csv $wkdir/"${dir##*/}"_concatenated.csv
done


#***********************************#
#                                   #
# Clustering                        #
#                                   #
#***********************************#

for dir in $fastq_dir
    do
    ./ashure.py clst -i "${dir##*/}"_concatenated.csv -o "${dir##*/}"_clusters.csv -iter $niter -r
done


#***********************************#
#                                   #
# Convert csv back to fasta         #
#                                   #
#***********************************#

for dir in $fastq_dir
do
cp $wkdir/"${dir##*/}"_clusters.csv input_clusters.csv
csv_to_fasta () {
    PYCMD=$(cat <<EOF
csv = 'input_clusters.csv'
output = 'output_clusters.fasta'

out_lines = []
temp_line = ''
with open(csv, 'r') as csv:
    for line in csv:
        cols = line.split(",")
        out_lines.append(temp_line)
        temp_line = ">" + cols[0] + "\n" + cols[1] + "\n"

out_lines.append(temp_line)

with open(output, 'w') as csv_out:
    csv_out.write(''.join(out_lines)[13:])
EOF
    )

    python3 -c "$PYCMD"
}
csv_to_fasta
mv output_clusters.fasta $wkdir/"${dir##*/}"_clusters.fasta
done 


#***********************************#
#                                   #
# Chimera detection                 #
#                                   #
#***********************************#

for dir in $fastq_dir
    do
    vsearch --uchime_denovo $wkdir/"${dir##*/}"_clusters.fasta --chimeras $wkdir/"${dir##*/}"_chimeras.fasta --nonchimeras $wkdir/"${dir##*/}"_nonchimeras.fasta
done


#***********************************#
#                                   #
# Convert fasta back to csv         #
#                                   #
#***********************************#

for dir in $fastq_dir
do
cp $wkdir/"${dir##*/}"_nonchimeras.fasta input.fasta
fasta_to_csv () {
    PYCMD=$(cat <<EOF
fasta = 'input.fasta'
output = 'output.csv'

out_lines = []
temp_line = ''
with open(fasta, 'r') as fp: 
    for line in fp:
        if line.startswith('>'):
            out_lines.append(temp_line)
            int_line = line.strip() + ','
            temp_line = int_line.strip('>')
        else:
            temp_line += line.strip()
out_lines.append(temp_line)

with open(output, 'w') as fp_out:
    fp_out.write('id,sequence' + '\n'.join(out_lines))
EOF
    )

    python3 -c "$PYCMD"
}
fasta_to_csv
mv output.csv $wkdir/"${dir##*/}"_nonchimeras.csv
done


#***********************************#
#                                   #
# ASV table                         #
#                                   #
#***********************************#

cd $wkdir

for dir in $fastq_dir
    do
    minimap2 -c --cs $wkdir/"${dir##*/}"_nonchimeras.fasta $wkdir/$resultsdir/qc/"${dir##*/}"/"${dir##*/}"_concatenated.fastq.gz > results.paf
    sed -i '1s/^/query id\tquery length\tquery start\tquery end\tstrand direction\ttarget id\ttarget length\ttarget start\ttarget end\tmatching bases\tlength match\tmapping quality\t\t\t\t\t\t\t\t\t\t\t\t\n/' results.paf
    
OTU_cleanup () {
    PYCMD=$(cat <<EOF

import pandas as pd

df_paf = pd.read_csv("results.paf", sep='\t')
pd.to_numeric(df_paf['matching bases'])
df_paf['percentage'] = df_paf['matching bases']/df_paf['length match']
sort = df_paf.sort_values(by=['percentage'], ascending=False)
wtduplicates = sort.drop_duplicates(subset=['query id'], keep='first')
counts = wtduplicates['target id'].value_counts()
clusters = pd.Index.to_series(counts.index)
df_counts = pd.concat([counts.rename('counts'), clusters.rename('clusters')], axis=1)    

df_counts.to_csv("table.csv", sep=',', index=False)

EOF
    )

    python3 -c "$PYCMD"
}
OTU_cleanup

    mv $wkdir/table.csv $resultsdir/"${dir##*/}"_table.csv
    
done


for dir in $fastq_dir
    do
    mv $wkdir/"${dir##*/}"_clusters.fasta $resultsdir/
    mv $wkdir/"${dir##*/}"_clusters.csv $resultsdir/
    mv $wkdir/"${dir##*/}"_concatenated.csv $resultsdir/
    mv $wkdir/"${dir##*/}"_chimeras.fasta $resultsdir/
    mv $wkdir/"${dir##*/}"_nonchimeras.fasta $resultsdir/
    mv $wkdir/"${dir##*/}"_nonchimeras.csv $resultsdir/
done    

    rm -r __pycache__
    mv ashure.log $resultsdir/
    rm $wkdir/input.clusters.csv


end=$SECONDS
duration=$(( end - start ))
echo "All good! It took GRIFO $(($duration/3600)) hours, $(($duration/60)) minutes and $(($duration%60)) seconds to complete the job."