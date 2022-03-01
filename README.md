# BRICO
BRICO is an analysis pipeline for nanopore metabarcoding sequence data

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

Clustering with {em}ashure{/em}

``` bash
for bc in barcode*
do
./ashure.py clst -i "$bc"_output.csv -o "$bc".csv -iter 3 -r
done
```
