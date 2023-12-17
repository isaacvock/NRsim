from Bio import SeqIO
import sys
import argparse

# Parse commandline arguments
parser = argparse.ArgumentParser(description='This is a script to convert a fasta file to a fastq file')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-f', '--fasta', type=str, required=True, metavar = 'in_file.fasta',
                    help='FASTA file to process')
parser.add_argument('--qscore', default='I', type=str,
                    help='qscore ASCII value to impute')

args = parser.parse_args()

# name without .fasta suffix
inputName = args.fasta.split('.fasta')[0] 

# fastq file path
fastq = inputName + '.fastq'


# Function that will process fasta file and produce fastq file
def convert_fasta_to_fastq(fasta_file, q_char, fastq_file):
    try:
        q_score = ord(q_char)  # Convert ASCII character to its ASCII value
        if q_score < 33 or q_score > 126:
            raise ValueError("Q-score ASCII value must be between 33 and 126.")

        with open(fastq_file, 'w') as fq:
            for record in SeqIO.parse(fasta_file, "fasta"):
                fq.write(f"@{record.id}\n")
                fq.write(str(record.seq) + "\n")
                fq.write("+\n")
                fq.write(q_char * len(record.seq) + "\n")

        print(f"FASTQ file generated: {fastq_file}")

    except Exception as e:
        print(f"Error: {e}")

# Convert fasta to fastq
convert_fasta_to_fastq(args.fasta, args.qscore, fastq)
