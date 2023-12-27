import pandas as pd
import random
from Bio import SeqIO
import argparse
from Bio.Seq import Seq

# Parse commandline arguments
parser = argparse.ArgumentParser(description='This is a script to convert a fasta file to a fastq file')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-f', '--fastq', type=str, required=True, metavar = 'in_file.fastq',
                    help='FASTQ file to introduce mutations in')
requiredNamed.add_argument('-k', '--kinetics', type=str, required=True, metavar = 'kinetics.csv',
                    help='csv file with kinetic parameters for each transcript to derive fraction new from')
parser.add_argument('-r', "--read", default=1, type=int, metavar = '<int>', choices=[1, 2],
                    help='Read 1 or 2?')
parser.add_argument('-m', '--mutrate', default=0.05, type=float,
                    help='s4U labeled read mutation rate')

args = parser.parse_args()

# name without .fasta suffix
inputName = args.fastq.split('.fastq')[0] 

# Output name
outputName = inputName + '.nrseq.fastq'


def modify_nucleotides(sequence, read, mutrate):
    """ Modify nucleotides in the sequence based on the given probability. """

    if read == 1:

        modified_sequence = ""
        for nucleotide in sequence:
            if nucleotide == 'T':
                modified_sequence += 'C' if random.random() < mutrate else 'T'
            else:
                modified_sequence += nucleotide
        return modified_sequence
    
    else:

        modified_sequence = ""
        for nucleotide in sequence:
            if nucleotide == 'G':
                modified_sequence += 'A' if random.random() < mutrate else 'G'
            else:
                modified_sequence += nucleotide
        return modified_sequence

def process_fastq(csv_file, fastq_file, output_fastq, read, mutrate):
    # Load the CSV file
    df = pd.read_csv(csv_file)
    transcript_to_fn = dict(zip(df['transcript_id'], df['fn']))

    # Process the FASTQ file
    with open(output_fastq, "w") as output_handle:
        for record in SeqIO.parse(fastq_file, "fastq"):
            # Find the transcript_id in the read name
            transcript_id = next((tid for tid in transcript_to_fn if tid in record.id), None)
            if transcript_id:
                
                # Generate Bernoulli random variable
                newness = random.random() < transcript_to_fn[transcript_id]
                if newness:
                    # Modify the sequence
                    record.seq = Seq(modify_nucleotides(str(record.seq), read = read, mutrate = mutrate))

            # Write the record to the output file
            SeqIO.write(record, output_handle, "fastq")

# Make NR-seq fastqs
process_fastq(args.kinetics, args.fastq, str(outputName), args.read, args.mutrate)
