#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
read=$3
qscore=$4
pyscript=$5
output=$6



parallel -j $cpus "python $pyscript -f {1} \
                                            --qscore $qscore" ::: ./results/split_fasta/sample_"$sample"_"$read".*.fasta


cat ./results/split_fasta/sample_"$sample"_"$read".*.fastq | pigz -p $cpus > "$output" 

