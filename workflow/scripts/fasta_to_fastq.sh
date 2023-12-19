#!/bin/bash

echo "Parse args"

# Source the paths and variables:
cpus=$1
sample=$2
qscore=$3
pyscript=$4
pyscript2=$5
kinetics=$6
output=$7
PE=$8

echo "Set safety mode"

# Exit immediately if any command returns a non-zero status
set -e


if [ "$PE" = "True" ]; then

    read=$9

    # Use seqtk to convert fasta file to fastq
    parallel -j $cpus "seqtk seq -F '$qscore' {1} > ./results/convert_to_fastq/sample_"$sample"_"$read".{#}.fastq" ::: ./results/simulate_fastas/sample_"$sample"_"$read".*.fasta



    # Introduce T-to-C mutations in fastq file to simulate NR-seq data
    parallel -j $cpus "python $pyscript2 -f {1} \
                                                -k $kinetics" ::: ./results/convert_to_fastq/sample_"$sample"_"$read".*.fastq

    # Clean up temp files
    rm -f ./results/convert_to_fastq/sample_"$sample"_"$read".*.fastq


    # Combine NR-seq fragment fastqs and gzip
    cat ./results/split_fasta/sample_"$sample"_"$read".*.nrseq.fastq | pigz -p $cpus > "$output" 

    # Clean up temp files
    rm -f ./results/split_fasta/sample_"$sample"_"$read".*.nrseq.fastq


else

    # Use seqtk to convert fasta to fastq
    parallel -j $cpus "seqtk seq -F '$qscore' {1} > ./results/convert_to_fastq/sample_"$sample".{#}.fastq" ::: ./results/simulate_fastas/sample_"$sample".*.fasta



    # Introduce T-to-C mutations in fastq file to simulate NR-seq data
    parallel -j $cpus "python $pyscript2 -f {1} \
                                                -k $kinetics" ::: ./results/convert_to_fastq/sample_"$sample".*.fastq


    # Clean up temp files
    rm -f ./results/convert_to_fastq/sample_"$sample".*.fastq

    # Combine NR-seq fragment fastqs and gzip
    cat ./results/convert_to_fastq/sample_"$sample".*.nrseq.fastq | pigz -p $cpus > "$output" 

    # Clean up temp files
    rm -f ./results/convert_to_fastq/sample_"$sample".*.nrseq.fastq



fi


