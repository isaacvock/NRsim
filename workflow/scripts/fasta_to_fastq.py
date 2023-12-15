from Bio import SeqIO
import sys

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

# User input for the file names and Q-score character
fasta_file = input("Enter the path to the FASTA file: ")
q_char = input("Enter the ASCII character for the Q-score: ")
fastq_file = input("Enter the path for the output FASTQ file: ")

convert_fasta_to_fastq(fasta_file, q_char, fastq_file)
