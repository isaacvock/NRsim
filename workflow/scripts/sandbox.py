from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os
import random

dict = {"a": {"ab": 1, "ac": 2}, "c": {"cb": 1, "c": 2}}

dict["a"].keys()


key = "a"

dict[key]

def modify_nucleotides(sequence):
    """ Modify nucleotides in the sequence based on the given probability. """
    modified_sequence = ""
    for nucleotide in sequence:
        if nucleotide == 'T':
            modified_sequence += 'C' if random.random() < 0.05 else 'T'
        else:
            modified_sequence += nucleotide
    return modified_sequence



with open("sandbox_nrseq.fastq", "w") as output_handle:

    for record in SeqIO.parse("sandbox.fastq", "fastq"):
        newness = random.random() < 0.5
        if newness:
                record.seq = Seq(modify_nucleotides(str(record.seq)))
        
        SeqIO.write(record, output_handle, "fastq")

