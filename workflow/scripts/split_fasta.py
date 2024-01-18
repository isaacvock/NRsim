from Bio import SeqIO
import math


# Process parameters
reads = snakemake.params.get("reads")
nsplit = snakemake.params.get("nsplit")
outdir = snakemake.params.get("outdir")



# Calculate the number of reads that should be in each file
reads_per_file = math.ceil(reads/nsplit)

# Input
fasta = snakemake.input.get("fasta")


# name without .fasta suffix
inputName = fasta.split('.fasta')[0] 


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
            current_entry_count = 0
            current_file = open(os.path.join(output_dir, f'{inputName}.{file_count}.fasta'), 'w')

        SeqIO.write(record, current_file, "fasta")
        current_entry_count += 1


    current_file.close()

split_fasta(fasta, outdir, nsplit, reads_per_file)