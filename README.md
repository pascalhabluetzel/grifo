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

## Basecalling

Basecalling with <em>guppy</em>. Specify flow cell type and sequencing kit.

``` bash
mkdir ./results/basecalling/
/opt/ont-guppy-cpu_3.0.3/bin/guppy_basecaller -i ./data/ -s ./results/basecalling/ --flowcell FLO-MIN106 --kit SQK-PSK004
```

## Demultiplexing

Demultiplexing with <em>qcat</em>

``` bash
mkdir ./results/demultiplexing/
cat ./results/basecalling/*.fastq.gz | qcat -b ./results/demultiplexing/
```

## Quality control and filtering

Filtering with <em>nanofilt</em>

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

``` bash
cd ./results/qc/
    for bc in barcode*
    do
    cd ./"$bc"/
    cat *.fastq | sed -n '1~4s/^@/>/p;2~4p' > "$bc"_concatenated.fasta
    done
```

## Conversion to .csv

```{r}
setwd("C:/Users/pascalh/Documents")
input <- readLines("barcode12_concatenated.fasta")
output <- file("barcode12_output.csv","w")

currentSeq <- 0
newLine <- 0

for(i in 1:length(input)) {
  if(strtrim(input[i], 1) == ">") {
    if(currentSeq == 0) {
      writeLines(paste(input[i],""), output, sep=",")
      currentSeq <- currentSeq + 1
    } else {
      writeLines(paste("\n",input[i],"", sep =""), output, sep=",")
    }
  } else {
    writeLines(paste(input[i]), output, sep="")
  }
}

close(output)
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

## Taxonomic assignment

Taxonomic assignement with <em>CREST4</em>

``` bash
for bc in barcode*
do
crest4 -f $bc".fasta -d bold
done
```
