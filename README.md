### To do list

- Write a bash or python code for what is now written in R (conversions between csv and fasta formats). - Done
- Write code to make an OTU table. - see "OTU table"
- Generalize basecalling and demultiplexing for most popular library preparation protocols.
- Make a configuration file. - Done

# BRICO
BRICO is an analysis pipeline for nanopore metabarcoding sequence data.
The name of the pipeline is derived of the French verb <em>bricoler</em> (to cobble sth) and describes how it integrates existing tools into a new pipeline.

## Dependencies

Brico uses the clustering module of <em>ashure</em>, which itself relies on the following tools:

``` bash
pip install pandas          # for organizing underlying data
pip install scikit-learn    # for clustering
pip install hdbscan         # for clustering
pip install spoa            # for clustering
```

## Configuration file

Parameters:<br>

``` bash
flow_cell="FLO-MIN106"
library_kit="SQK-PSK004"
locus="COI" # "COI", "18S", "16S", "ITS", "rbcl"
quality_filtering="8" # minimum quality score for Nanofilt
minimum_length="600" # minimum sequence length
maximum_length="800" # maximum sequence length
nit="10" # number of iterations for the clustering algorithm
```

## Doing checks on the data and environment

``` bash
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
source ./*.cfg
```

## Check input files
 
``` bash
fastq_count=$(find . -name "*.fastq.gz" | wc -l)
fastq_dir=$(find . -type f -name "*.fastq.gz" | sed -r 's|/[^/]+$||' |sort |uniq)
fast5_count=$(find . -name "*.fast5" | wc -l)
fast5_dir=$(find . -type f -name "*.fast5" | sed -r 's|/[^/]+$||' |sort |uniq)

if [ $fastq_count -gt 0 ]
then
      echo "found $fastq_count g-zipped fastq files in $fastq_dir, initiating filtering"
elif [ $fast5_count -gt 0 ]
then
      echo  "found $fast5_count fast5 files in $fast5_dir, initiating base calling"
      for val in $fast5_dir
      do
      mkdir -p ./results/basecalling/$val
      /opt/ont-guppy-cpu_3.0.3/bin/guppy_basecaller -i .$val/ -s ./results/basecalling/$val --flowcell $flow_cell --kit $library_kit
      done
else
      echo "\nERROR: BRICO needs sequences in .fast5 or .fastq format as input."
fi
```


## Basecalling

Basecalling with <em>guppy</em>. Specify flow cell type and sequencing kit.<br>
Input: fast5 files<br>
Output: fastq files

``` bash
mkdir ./results/basecalling/
/opt/ont-guppy-cpu_3.0.3/bin/guppy_basecaller -i ./data/ -s ./results/basecalling/ --flowcell $flow_cell --kit $library_kit
```

## Demultiplexing

Demultiplexing with <em>qcat</em><br>
Input: fastq files<br>
Output: fastq files split in different folders by sample barcode

``` bash
mkdir ./results/demultiplexing/
cat ./results/basecalling/*.fastq.gz | qcat -b ./results/demultiplexing/
```

## Quality control and filtering

Filtering with <em>nanofilt</em><br>
Input: unfiltered fastq files<br>
Output: filtered fastq files

``` bash
mkdir ./results/qc/
cd ./results/demultiplexing/
find . -type d -exec mkdir -p -- ./results/qc/{}

for bc in barcode* ;
    do
    cd ./barcode"$bc"/
        for fastq in *.fastq.gz ;  
            do
            echo "Filtering data $fastq..."
            gunzip -c "$fastq" | NanoFilt --length 500 --maxlength 800 -q 8 | gzip > ./results/qc/barcode"$bc"/"$fastq"
        done
done
```

## Concatenate

Input: multiple fastq files<br>
Output: a single fasta file

``` bash
cd ./results/qc/
for bc in barcode*
    do
    cd ./"$bc"/
    cat *.fastq | sed -n '1~4s/^@/>/p;2~4p' > "$bc"_concatenated.fasta
done
```

## Conversion to .csv

~~setwd("C:/Users/pascalh/Documents")<br>
input <- readLines("barcode12_concatenated.fasta")<br>
output <- file("barcode12_output.csv","w")<br>
currentSeq <- 0<br>
newLine <- 0<br>
for(i in 1:length(input)) {<br>
  if(strtrim(input[i], 1) == ">") {<br>
    if(currentSeq == 0) {<br>
      writeLines(paste(input[i],""), output, sep=",")<br>
      currentSeq <- currentSeq + 1<br>
    } else {<br>
      writeLines(paste("\n",input[i],"", sep =""), output, sep=",")<br>
    }<br>
  } else {<br>
    writeLines(paste(input[i]), output, sep="")<br>
  }<br>
}<br>
close(output)~~<br>
\
~~for bc in barcode*<br>
    do<br>
    sed -i '1s/^/id,sequence\n/' "$bc"_output.csv #add column names 'id' and 'sequence'<br>
done~~<br>


```python

fasta = 'TestFASTA.fasta'
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
    
 ```

## Clustering

Clustering with <em>ashure</em>

``` bash
for bc in barcode*
    do
    ./ashure.py clst -i "$bc"_output.csv -o "$bc".csv -iter 3 -r
done
```

## Conversion to .fasta

```{r}
setwd("C:/Users/pascalh/Documents")
#BiocManager::install("Biostrings")
library(Biostrings)
csv = read.csv("barcode12.csv")
seq = csv$sequence
names(seq) = csv$id
dna = DNAStringSet(seq)
writeXStringSet(dna, "barcode12.fasta")
```

something like this in Python:

```python
lines = open("input.csv").readlines()

for line in lines:
    cols = line.split("\t")
    print(">"+cols[0]+"\n"+cols[1])
```

Python by Els - updated

```python
csv = 'TestCSV.csv'
output = 'outputFASTA2.fasta'


out_lines = []
temp_line = ''
with open(csv, 'r') as csv:
    for line in csv:
        cols = line.split("\t")
        out_lines.append(temp_line)
        temp_line = ">" + cols[0] + "\n" + cols[1] + "\n"

out_lines.append(temp_line)

with open(output, 'w') as csv_out:
    csv_out.write(''.join(out_lines)[13:])
```

## Taxonomic assignment

Taxonomic assignement with <em>CREST4</em>

``` bash
for bc in barcode*
    do
    crest4 -f "$bc".fasta -d bold
done
```

## OTU table

Use minimap2 to align all sequences with the cluster centroids. Calculate a mapping percentage by dividing column 10 of the .paf file by column 11. If there are different matches for forward and reverse version of the sequence, remove the one with the smaller mapping percentage. Remove all sequences with a mapping percentage below a specific value (e.g. 0.9). Count the occurrence of each cluster in column 6. Merge this information with the best blast hit for each cluster. The resulting table has 4 columns: cluster number, number of reads, alignment percentage with reference database (SILVA, BOLD, etc.), taxonomic assignment.

``` bash
for dir in $fastq_dir
    do
    minimap2 -c --cs "${dir##*/}"_clusters.fasta $wkdir/results/qc/"${dir##*/}"/"${dir##*/}"_concatenated.fastq.gz > "${dir##*/}".paf
done
```
