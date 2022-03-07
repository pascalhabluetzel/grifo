### To do list

- Write a bash or python code for what is now written in R (conversions between csv and fasta formats). - Done
- Write code to make an OTU table.
- Generalize basecalling and demultiplexing for most popular library preparation protocols.
- Make a configuration file

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
flow_cell = "FLO-MIN106"
library_kit = "SQK-PSK004"
locus = "COI" # "COI", "18S", "16S", "ITS", "rbcl"
quality_filtering = "8" # minimum quality score for Nanofilt
minimum_length = "600" # minimum sequence length
maximum_length = "800" # maximum sequence length
nit = "10" # number of iterations for the clustering algorithm
```

## Check input files
 
``` bash
#!/bin/bash

file_count=$(find . -name "*.fastq.gz" | wc -l)
fast5_count=$(find . -name "*.fast5" | wc -l)

if [ $file_count -gt 0 ]
then
                echo "found $file_count fasta files, initiating filtering"
elif [ $fast5_count -gt 0 ]
then
        echo  "found $fast5_count fast5 files, initiating base calling"
else
        echo "no input data"
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
