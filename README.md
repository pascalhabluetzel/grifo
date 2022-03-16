# GRIFO
GRIFO is a versatile pipeline for environmental metabarcoding using nanopore sequences with density based clustering for error correction and other utilities.
The name is derived of the mythical creature <em>grifo</em> (Esperanto) which is partly eagle, lion and deer. The name describes how the pipeline integrates existing tools into a new workflow.

---
**NOTE**

It works with almost all markdown flavours (the below blank line matters).

---

## Dependencies

GRIFO uses the clustering module of <em>ashure</em>, which itself relies on the following tools:

``` bash
pip install pandas          # for organizing underlying data
pip install scikit-learn    # for clustering
pip install hdbscan         # for clustering
pip install spoa            # for clustering
```

## Configuration file

Parameters:<br>

``` bash
demultiplexed="yes" # 'yes' or 'no'
flow_cell="FLO-MIN106"
library_kit="SQK-PSK004"
minimum_length="600" # minimum sequence length
maximum_length="800" # maximum sequence length
qscore="8"
niter="3"
locus="COI" # 'COI', '18S' or '16S'
```

## Doing checks on the data and environment

``` bash
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
```

## Check input files

GRIFO looks for the input files only in the ./data/ folder to avoid that it finds the sequences in the results folder. Make sure the data is in a folder called ./data/

``` bash
if [ $demultiplexed = yes ]
then
      echo "Sequences are already demultiplexed."
elif [ $demultiplexed != yes ]
then
      echo "Check input file type."
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
```

## Quality control and filtering

Filtering with <em>nanofilt</em><br>
Input: unfiltered fastq files<br>
Output: filtered fastq files

``` bash
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
```

## Conversion to .csv

``` bash
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
 ```

## Clustering

Clustering with <em>ashure</em>

``` bash
for dir in $fastq_dir
    do
    ./ashure.py clst -i "${dir##*/}"_concatenated.csv -o "${dir##*/}"_clusters.csv -iter $niter -r
done
```

## Conversion to .fasta

``` bash
for dir in $fastq_dir
do
mv $wkdir/"${dir##*/}"_clusters.csv input_clusters.csv
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

rm -r __pycache__; rm input_clusters.csv; rm input.fasta
mv ashure.log ./$resultsdir/ashure.log
```

## Taxonomic assignment

Taxonomic assignement with <em>CREST4</em>

``` bash
if [ $locus == "COI" ]
then database="bold"
else
echo "No reference database for this locus."
fi

echo "Using $database database for taxonomic assignment."

for file in *_clusters.fasta
    do
    crest4 -f $file -d $database
done
```

## OTU table

Use minimap2 to align all sequences with the cluster centroids. Calculate a mapping percentage by dividing column 10 of the .paf file by column 11. If there are different matches for forward and reverse version of the sequence, remove the one with the smaller mapping percentage. Remove all sequences with a mapping percentage below a specific value (e.g. 0.9). Count the occurrence of each cluster in column 6. Merge this information with the best blast hit for each cluster. The resulting table has 4 columns: cluster number, number of reads, alignment percentage with reference database (SILVA, BOLD, etc.), taxonomic assignment.

``` bash
cd $wkdir

for dir in $fastq_dir
    do
    cp "${dir##*/}"_clusters.fasta.crest4/search.hits ./
    cp "${dir##*/}"_clusters.fasta.crest4/assignments.txt ./
    minimap2 -c --cs "${dir##*/}"_clusters.fasta $wkdir/$resultsdir/qc/"${dir##*/}"/"${dir##*/}"_concatenated.fastq.gz > results.paf

OTU_cleanup () {
    PYCMD=$(cat <<EOF

import pandas

df_assignments = pandas.read_csv("assignments.txt", header=None, sep='\t')
df_assignments.columns = ['k', 'taxon']
df_paf = pandas.read_csv("results.paf", header=None, sep='\t')
df_paf = df_paf.iloc[:, [0, 5, 9, 10]]
df_paf.columns = ['read', 'cluster', 'value1', 'value2']

df_paf['percentage'] = ( df_paf['value1'].T / df_paf['value2'] ).T

result0 = df_paf.sort_values(by=['percentage'], ascending=False)
result0.read.duplicated()
result = result0.drop_duplicates(subset=['read'], keep='first')
x = result['cluster'].value_counts()
print(x)
y = pandas.DataFrame(x)

y['count'] = y.index
y.columns=['cluster', 'k']



with open('search.hits', 'r') as reader:
    # Note: readlines doesn't trim the line endings
    text = reader.readlines()
with open('search.hist.mod', 'w') as writer:
    for line in text:
        if (line.startswith('c')):
            writer.write(line)

df_hits = pandas.read_csv("search.hist.mod", header=None, sep='\t')
df_hits.columns = ['k', 'hit', 'score', 'alignment', 'match']
clean_hits = df_hits.drop_duplicates(subset=['k'], keep='first')
clean_hits['percentage'] = (clean_hits['match'].T / clean_hits['alignment']).T

a = pandas.merge(y, clean_hits, on=['k'])
results = pandas.merge(a, df_assignments, on=['k'])

results.to_csv("results.tsv", sep='\t')

EOF
    )

    python3 -c "$PYCMD"
}
OTU_cleanup

rm assignments.txt search.hits
mv results.paf ./$resultsdir/results.paf

    mv results.tsv ./$resultsdir/"${dir##*/}"_OTU_table.tsv
done

for dir in $fastq_dir
    do
    mv $wkdir/"${dir##*/}"_clusters.fasta.crest4/ $resultsdir/
    mv $wkdir/"${dir##*/}"_clusters.fasta $resultsdir/
    mv $wkdir/"${dir##*/}"_concatenated.csv $resultsdir/
done
```

# Finalize

``` bash
cp ./*.cfg ./$resultsdir/*.cfg # Copy the config file into the results folder for documentation.
end=$SECONDS
duration=$(( end - start ))
echo "All good! It took GRIFO $(($duration/3600)) hours, $(($duration/60)) minutes and $(($duration%60)) seconds to complete the job."
```


