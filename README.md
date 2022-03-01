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
cd ./results/demultiplexing/
for bc in barcode* ;
    cd ./barcode"$bc"/
        for fastq in *.fastq.gz ;  
        do
        echo "Filtering data $fastq..."
        gunzip -c "$fastq" | NanoFilt --length 500 --maxlength 800 -q 8 | gzip > ~/Nanopore_metabarcoding/results/qc/"$fastq"
        done
```
