#!/usr/bin/env python

# Loading libraries
import configparser
import gzip
import os
import glob
import subprocess
import multiprocessing
import shutil
import pandas as pd

'''
Function to obtain paths to the directories where the fastq files are located.
Currently requires demultiplexed data with gzipped files.
'''
def fastq_directory_paths():
    fastq_dir = set()
    file_paths = glob.glob(path_data + '/**/*.fastq.gz', recursive=True)
    directory_paths = map(os.path.dirname, file_paths)
    fastq_dir.update(directory_paths)
    return sorted(fastq_dir)

'''
Function to obtain sample names that are used as base to track the samples through the analysis
'''
basename = lambda paths: [path.split('/')[-2] if path.endswith('/') else path.split('/')[-1] for path in paths]

'''
Functions to count the number of raw sequences per sample
'''
def count_sequences(file_path):
    count = 0
    with gzip.open(file_path, 'rt') as gz_file:
        for line in gz_file:
            if line.startswith('@'):
                count += 1
    return count

def accumulate_counts(folder_path):
    counts = {}
    for folder in folder_path:
        counts[folder] = 0  # Initialize the folder key with count 0
        for file_name in os.listdir(folder):
            file_path = os.path.join(folder, file_name)
            if file_name.endswith('.fastq.gz'):
                counts[folder] += count_sequences(file_path)
    return counts

'''
Wrapper function for NanoFilt
'''
def filtering(folder_list, minimum_length, maximum_length, qscore, base_list):
    for folder, base in zip(folder_list, base_list):
        if os.path.exists(os.path.join('results/qc', base)) and os.path.isdir(os.path.join('results/qc', base)):
            pass
        else:
            os.makedirs(os.path.join('results/qc', base))
        for file_name in os.listdir(folder):
            if file_name.endswith('.fastq.gz'):
                file_path = os.path.join(folder, file_name)
                file_name_unzipped = file_name[:-3]
                command = f'gunzip -c {file_path} | NanoFilt \
                          --length {minimum_length} \
                          --maxlength {maximum_length} \
                          -q {qscore} | gzip > ./results/qc/{base}/{file_name}'
                subprocess.run(command, shell=True)
                
'''
Function to concatenate filtered sequences to one file per sample
'''
def concatenate(bases):
    for base in bases:
        command = f'cat ./results/qc/{base}/*.fastq.gz > ./results/qc/{base}/{base}_concatenated.fastq.gz'
        subprocess.run(command, shell=True)
        
        command = f'gunzip ./results/qc/{base}/{base}_concatenated.fastq.gz'
        subprocess.run(command, shell=True)
        
        command = f'sed -n "1~4s/^@/>/p;2~4p" ./results/qc/{base}/{base}_concatenated.fastq > ./results/qc/{base}/{base}_concatenated.fasta'
        subprocess.run(command, shell=True)

        #command = f'gzip ./results/qc/{base}/{base}_concatenated.fasta'
        #subprocess.run(command, shell=True)
        
'''
Function to count the number of sequences per sample after filtering
'''
def count_sequences_concat(base_name):
    counts = {}
    for base in base_name:
        counts[base] = 0  # Initialize the count for each base
        file_path = './results/qc/' + base + '/' + base + '_concatenated.fasta'
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    counts[base] += 1
    return counts

'''
Function to convert sequence files in fasta format to csv
'''
def fasta2csv(base_name):
    for base in base_name:
        fasta = './results/qc/' + base + '/' + base + '_concatenated.fasta'
        output = './results/qc/' + base + '/' + base + '_concatenated.csv'

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
            
'''
Wrapper function to run the ashure clustering algorithm
'''
def cluster(base):
    os.chdir(wdir)
    shutil.copy('./ashure.py', './results/qc/' + base)
    shutil.copy('./bilge_pype.py', './results/qc/' + base)
    os.chdir('./results/qc/' + base)
    script_path = "./ashure.py"
    input_file = base + "_concatenated.csv"
    output_file = base + "_clusters.csv"    
    
    command = [
        script_path,
        "clst",
        "-i", input_file,
        "-o", output_file,
        "-iter", config.get('ashure', 'niter'),
        "-r"
    ]   
    
    subprocess.run(command)
    os.chdir(wdir)

'''
Function to convert sequence files in csv format to fasta
'''
def csv2fasta(base_name):
    for base in base_name:
        csv_path = './results/qc/' + base + '/' + base + '_clusters.csv'
        output_path = './results/qc/' + base + '/' + base + '_clusters.fasta'

        if os.path.exists(csv_path):
            out_lines = []
            temp_line = ''
            with open(csv_path, 'r') as csv_file:
                for line in csv_file:
                    cols = line.split(",")
                    out_lines.append(temp_line)
                    temp_line = ">" + cols[0] + "\n" + cols[1] + "\n"

            out_lines.append(temp_line)

            with open(output_path, 'w') as csv_out:
                csv_out.write(''.join(out_lines)[13:])
                
'''
Wrapper function to run cutadapt
'''
def remove_primers(base_name):
    for base in base_name:
        file_path = './results/qc/' + base + '/' + base + '_clusters.fasta'
        out_path = './results/qc/' + base + '/' + base + '_clusters_cut.fasta'
        if os.path.exists(file_path):
            command = [
                'cutadapt',
                '-a', 'CAGCAGCCGCGGTAATTCC;max_error_rate=0.20',
                '-g', 'CCCGTGTTGAGTCAAATTAAGC;max_error_rate=0.20',
                '--revcomp',
                '-o', out_path,
                file_path
            ]
            subprocess.run(command)

'''
Wrapper function to run blastn
'''                        
def blast(base_name):
    for base in base_name:
        file_path = './results/qc/' + base + '/' + base + '_clusters_cut.fasta'
        db = config.get('BLAST', 'db')
        if os.path.exists(file_path):
            print("Running blastn on", base)
            output_csv = './results/qc/' + base + '/' + base + '_' + db + '_blastn.csv'

            command = [
                "blastn",
                "-db", config.get('paths', 'path_to_blastdb'),
                "-query", file_path,
                "-task", "blastn",
                "-dust", "no",
                "-num_threads", str(config.get('BLAST', 'numthreads')),
                "-outfmt", "7 delim=, sseqid stitle qacc sacc evalue bitscore length pident",
                "-max_target_seqs", str(config.get('BLAST', 'mts')),
                "-perc_identity", str(config.get('BLAST', 'pct_ident')),
                "-out", output_csv
            ]

            subprocess.run(command)

'''
Function to handle the blastn output files and generate a concatenated table with the taxonomic annotations
'''
def make_output_file(base_name):
    db = config.get('BLAST', 'db')
    for base in base_name:
        input_csv = './results/qc/' + base + '/' + base + '_' + db + '_blastn.csv'
        output_csv = './results/qc/' + base + '/' + base + '_' + db + '_blastn2.csv'
        if os.path.exists(input_csv):
            with open(input_csv, 'r') as infile, open(output_csv, 'w') as outfile:
                for line in infile:
                    if not line.startswith('#'):
                        outfile.write(line)
                    
    for base in base_name:
        input_csv = './results/qc/' + base + '/' + base + '_' + db + '_blastn2.csv'
        print(input_csv)
        output_csv = './results/qc/' + base + '/' + base + '_' + db + '_ASV.csv'
        if os.path.exists(input_csv) and os.path.getsize(input_csv) > 0:
            # load file
            df = pd.read_csv(input_csv, sep=',')

            # add column names
            df.columns=['accession', 'taxonomic_annotation', 'cluster', 'accession', 'evalue', 'bitscore', 'alignment_length', 'percentage_identity']

            # select only rows with alignment length >= 500 bp
            df2 = df[df['alignment_length'] >= 500]

            # arrange rows by match percentage
            df3 = df2.sort_values(by=['percentage_identity'], ascending=False)

            # keep only first row of each ASV
            df4 = df3.drop_duplicates(subset=['cluster'], keep='first', inplace=False, ignore_index=False)

            # add sample name information
            df4['#sample_name'] = base

            df4['taxonomy'] = df4['taxonomic_annotation'].replace('"', '')

            df5 = df4[['#sample_name', 'cluster', 'accession', 'evalue', 'bitscore', 'alignment_length', 'percentage_identity', 'taxonomic_annotation']]

            df5.to_csv(output_csv, sep=';', index=False, header=False)

    intermediate = '_' + db + '_eDNA.csv'
    final = db + '_eDNA.csv'
    if os.path.exists(intermediate):
        os.remove(intermediate)

    for base in base_name:
        file_path = './results/qc/' + base + '/' + base + '_' + db + '_ASV.csv'
        if os.path.exists(file_path):
            with open(file_path, "r") as input_file, open(intermediate, "a") as output_file:
                output_file.write(input_file.read())

    with open(intermediate, "r") as input_file, open(final, "w") as output_file:
        output_file.write("counts,cluster,accession,accession,evalue,bitscore,alignment_length,percentage_identity,taxonomic_annotation\n")
        for line in input_file:
            if not line.startswith("#"):
                output_file.write(line.replace(";", ",").replace("|", ","))




config = configparser.ConfigParser()
config.read('config.ini')
path_data = config.get('paths', 'path_data')





# Load paths to directories where fastq files are located
fastq_dir = fastq_directory_paths()
print(fastq_dir)

# Extract sample names from data
base_name = basename(fastq_dir)
print(base_name)

# Count number of raw sequences per sample
print(accumulate_counts(fastq_dir))

# Filter sequences using NanoFilt
filtering(fastq_dir, config.get('NanoFilt', 'minlength'), config.get('NanoFilt', 'maxlength'), config.get('NanoFilt', 'qscore'), base_name)

# Concatenate filtered sequences to one file per sample
concatenate(base_name)
        
# Count number of sequences per sample after filtering
print(count_sequences_concat(base_name))

# Convert fasta files to csv
fasta2csv(base_name)

# Make sure you save the original working directory before using moving around directories
wdir = os.getcwd()
print(wdir)

# Run the actual clustering
if __name__ == '__main__':
    num_processes = 8  # Number of available CPU cores 
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(cluster, base_name)

# Move back to original working directory after clustering
os.chdir(wdir)

# Convert csv files to fasta
csv2fasta(base_name)

# Remove primers using a wrapper function for cutadapt
remove_primers(base_name)

# Taxonomic annotation using blastn
blast(base_name)

# Handle the blastn output files and generate a concatenated table with the taxonomic annotations
make_output_file(base_name) 
    
                    
