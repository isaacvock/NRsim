#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
qscore=$3
pyscript=$4
output=$5
read=$6



parallel -j $cpus "python $pyscript -f {1} \
                                            --qscore $qscore" ::: ./results/split_fasta/sample_"$sample"_"$read".*.fasta


cat ./results/split_fasta/sample_"$sample"_"$read".*.fastq | pigz -p $cpus > "$output" 

