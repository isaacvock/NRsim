#!/bin/bash

echo "Parse args"

# Source the paths and variables:
cpus=$1
sample=$2
qscore=$3
pyscript=$4
kinetics=$5
output=$6
PE=$7
dir=$8
mutrate=$9
seed=${10}

echo "Set safety mode"

# Exit immediately if any command returns a non-zero status
set -e


if [ "$PE" = "True" ]; then

    read=${10}

    # Use seqtk to convert fasta file to fastq
    parallel -j $cpus "seqtk seq -F '$qscore' {1} > "$dir"/sample_"$sample"_"$read".{#}.fastq" ::: "$dir"/sample_"$sample"_"$read".*.fasta



    # Introduce T-to-C mutations in fastq file to simulate NR-seq data
    parallel -j $cpus "python $pyscript -f {1} \
                                                -k $kinetics -r $read -m $mutrate -s $seed" ::: "$dir"/sample_"$sample"_"$read".*.fastq


    # Combine NR-seq fragment fastqs and gzip
    cat "$dir"/sample_"$sample"_"$read".*.nrseq.fastq | pigz -c -p $cpus > "$output" 

    ## Clean up temp files
    rm -f "$dir"/sample_"$sample"_"$read".*.fastq


else

    # Use seqtk to convert fasta to fastq
    parallel -j $cpus "seqtk seq -F '$qscore' {1} > "$dir"/sample_"$sample".{#}.fastq" ::: "$dir"/sample_"$sample".*.fasta



    # Introduce T-to-C mutations in fastq file to simulate NR-seq data
    parallel -j $cpus "python $pyscript -f {1} \
                                                -k $kinetics -m $mutrate" ::: "$dir"/sample_"$sample".*.fastq


    # Combine NR-seq fragment fastqs and gzip
    cat "$dir"/sample_"$sample".*.nrseq.fastq | pigz -c -p $cpus > "$output" 

    # Clean up temp files
    rm -f "$dir"/sample_"$sample".*.fastq



fi


