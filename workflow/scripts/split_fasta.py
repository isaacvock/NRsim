#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import math
import os


# Parse commandline arguments
parser = argparse.ArgumentParser(description='This is a script to convert a fasta file to a fastq file')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-f', '--fasta', type=str, required=True, metavar = 'in_file.fasta',
                    help='FASTA file to split up')
requiredNamed.add_argument('-d', '--dir', type=str, required=True, metavar = 'path/to/output_dir/',
                    help='Directory in which to save files')
requiredNamed.add_argument('-i', '--filename', type=str, required=True, metavar = 'filename',
                    help='Fasta filename without rest of directory')
requiredNamed.add_argument('-r', '--reads', required=True, type=int,
                    help='Number of files to split the FASTA into')
parser.add_argument('-n', '--nsplit', default=32, type=int,
                    help='Number of files to split the FASTA into')


args = parser.parse_args()

# Calculate the number of reads that should be in each file
reads_per_file = math.ceil(((args.reads)/(args.nsplit))/1.5)

# name without .fasta suffix
inputName = args.filename

def split_fasta(input_file, output_dir, num_split_files, reads_per_file):
    """
    Splits a FASTA file into multiple smaller files.

    :param input_file: Path to the input FASTA file.
    :param output_dir: Directory where the output files will be saved.
    :param num_split_files: Number of files to split the FASTA file into
    :param reads_per_file: Rough estimate for number of reads that should go in each split FASTA.
        The last file may have more or less than this number of reads depending on how bad the
        estimated total number of reads is. 
    """

    file_count = 1
    current_entry_count = 0
    current_file = open(os.path.join(output_dir, f'{inputName}.{file_count}.fasta'), 'w')

    for record in SeqIO.parse(input_file, "fasta"):

        if file_count < num_split_files and current_entry_count >= reads_per_file:
            current_file.close()
            file_count += 1
            print(current_entry_count)
            current_entry_count = 0
            current_file = open(os.path.join(output_dir, f'{inputName}.{file_count}.fasta'), 'w')

        SeqIO.write(record, current_file, "fasta")
        current_entry_count += 1

    print(file_count)
    current_file.close()

split_fasta(args.fasta, args.dir, args.nsplit, reads_per_file)