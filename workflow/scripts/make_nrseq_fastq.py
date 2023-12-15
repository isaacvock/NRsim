import pandas as pd
import random
from Bio import SeqIO

def modify_nucleotides(sequence):
    """ Modify nucleotides in the sequence based on the given probability. """
    modified_sequence = ""
    for nucleotide in sequence:
        if nucleotide == 'T':
            modified_sequence += 'C' if random.random() < 0.05 else 'T'
        else:
            modified_sequence += nucleotide
    return modified_sequence

def process_fastq(csv_file, fastq_file, output_fastq):
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
                    record.seq = modify_nucleotides(str(record.seq))

            # Write the record to the output file
            SeqIO.write(record, output_handle, "fastq")

# Example usage
process_fastq("transcripts.csv", "reads.fastq", "modified_reads.fastq")
