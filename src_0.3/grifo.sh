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
# Taxonomic assignment              #
#                                   #
#***********************************#

for file in *_nonchimeras.fasta
    do
    base=$(basename ${file} _nonchimeras.fasta)
    echo "Running blastn on" ${base}
    blastn \
    -db $path_to_blastdb \
    -query ${file} \
    -task blastn \
    -dust no \
    -num_threads $numthreads
    -outfmt "7 delim=, sseqid stitle qacc sacc evalue bitscore length pident" \
    -max_target_seqs $mts \
    -perc_identity $pct_ident \
    -out ${base}_blastn.csv
done

#***********************************#
#                                   #
# ASV table                         #
#                                   #
#***********************************#

export path_to_otu_taxonomy_file

cd $wkdir

for dir in $fastq_dir
	do
# Use minimap2 to map the raw sequences on the cluster consensuses to obtain counts
    minimap2 -c --cs $wkdir/"${dir##*/}"_nonchimeras.fasta $wkdir/$resultsdir/qc/"${dir##*/}"/"${dir##*/}"_concatenated.fastq.gz > results.paf
    sed -i '1s/^/query id\tquery length\tquery start\tquery end\tstrand direction\ttarget id\ttarget length\ttarget start\ttarget end\tmatching bases\tlength match\tmapping quality\t\t\t\t\t\t\t\t\t\t\t\t\n/' $wkdir/results.paf
    sed '/^#/d' "${dir##*/}"_blastn.csv > input.csv
    echo "${dir##*/}" > "${dir##*/}"_name
    cp "${dir##*/}"_name name

clean_up () {
    PYCMD=$(cat <<EOF

import pandas as pd

# load file
df = pd.read_csv("input.csv", sep=',')

df_paf = pd.read_csv("results.paf", sep='\t', usecols=range(12))
pd.to_numeric(df_paf['matching bases'])
df_paf['percentage'] = df_paf['matching bases']/df_paf['length match']
sort = df_paf.sort_values(by=['percentage'], ascending=False)
wtduplicates = sort.drop_duplicates(subset=['query id'], keep='first')
counts = wtduplicates['target id'].value_counts()
clusters = pd.Index.to_series(counts.index)
df_counts = pd.concat([counts.rename('counts'), clusters.rename('clusters')], axis=1)

# add column names
df.columns=['accession', 'taxonomic_annotation', 'cluster', 'accession2', 'evalue', 'bitscore', 'alignment_length', 'percentage_identity']

# select only rows with alignment length >= 1200 bp
df2 = df[df['alignment_length'] >= 1200]

# arrange rows by match percentage
df3 = df2.sort_values(by=['percentage_identity'], ascending=False)

# keep only first row of each ASV
df4 = df3.drop_duplicates(subset=['cluster'], keep='first', inplace=False, ignore_index=False)

# add sample name information
name = open("name","r")  
name = name.read().replace('\n','')  
df4['sample_name'] = name

df4['taxonomy'] = df4['taxonomic_annotation'].replace('"', '')

df5 = df4[['sample_name', 'cluster', 'accession', 'evalue', 'bitscore', 'alignment_length', 'percentage_identity', 'taxonomic_annotation']]

df_counts.columns=['counts', 'cluster']
df_counts_new = df_counts[['counts', 'cluster']]
df6 = pd.merge(df_counts_new, df5, on=['cluster'], how='outer')
df7 = df.drop(['accession', 'accession2'], axis=1)
df7.rename(columns={'taxonomic_annotation': 'accession'}, inplace=True)

df6.to_csv("output.csv", sep=';', index=False, header=False)

import os

otu_taxonomy_path = os.environ.get("path_to_otu_taxonomy_file")

df_tax = pd.read_csv(otu_taxonomy_path, sep='\t')
df_tax.columns=['accession', 'scientific_names']
df7 = pd.merge(df_tax, df6, on=['accession'], how='inner')
df7 = df7.drop(['taxonomic_annotation'], axis=1)
df7[['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']] = df7['scientific_names'].str.split('; ', expand=True)
df7[['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']] = df7[['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']].apply(lambda x: x.str[3:])

df7.to_csv("output.csv", sep=';', index=False, header=True)

EOF
    )

    python3 -c "$PYCMD"
}
clean_up


#***********************************#
#                                   #
# Cleaning up working environment   #
#                                   #
#***********************************#

    mv output.csv ./$resultsdir/"${dir##*/}"_asv_table.csv
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


#***********************************#
#                                   #
# Darwin Core Archive output        #
#                                   #
#***********************************#
'''
mkdir -p $resultsdir/src
    cp conversion_code.py ./$resultsdir/src/
    cp WoRMS.py ./$resultsdir/src/
mkdir -p $resultsdir/raw
    cp metadata_table.csv ./$resultsdir/raw/
mkdir -p $resultsdir/processed
for dir in $fastq_dir
    do
    cp ./$resultsdir/"${dir##*/}"_asv_table.csv ./$resultsdir/raw/asv_table.csv
    cp ./$resultsdir/"${dir##*/}"_taxa_table.csv ./$resultsdir/raw/taxa_table.csv
    cd ./$resultsdir/src/
    python ./conversion_code.py
    cd $wkdir
    cp ./$resultsdir/processed/dna_extension.csv ./$resultsdir/processed/"${dir##*/}"_dna_extension.csv
    cp ./$resultsdir/processed/occurrence.csv ./$resultsdir/processed/"${dir##*/}"_occurrence.csv     
done    
'''    

#***********************************#
#                                   #
# Finishing                         #
#                                   #
#***********************************#

end=$SECONDS
duration=$(( end - start ))
echo "All good! It took GRIFO $(($duration/3600)) hours, $(($duration%3600/60)) minutes and $(($duration%60)) seconds to complete the job."
