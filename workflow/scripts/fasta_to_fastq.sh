#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
qscore=$3
pyscript=$4
pyscript2=$5
kinetics=$6
output=$7
PE=$8

# Exit immediately if any command returns a non-zero status
set -e


if [ "$PE" = "True" ]; then

    read=$9

    # parallel -j $cpus "python $pyscript -f {1} \
    #                                             --qscore $qscore" ::: ./results/split_fasta/sample_"$sample"_"$read".*.fasta
    parallel -j $cpus "seqtk seq -F '\"$qscore\"' \
                        {1} > ./results/split_fasta/sample_\"$sample\"_\"$read\".{#}.fastq" ::: ./results/split_fasta/sample_"$sample"_"$read".*.fasta



    parallel -j $cpus "python $pyscript2 -f {1} \
                                                -k $kinetics" ::: ./results/split_fasta/sample_"$sample"_"$read".*.fastq


    cat ./results/split_fasta/sample_"$sample"_"$read".*.nrseq.fastq | pigz -p $cpus > "$output" 


else

    # parallel -j $cpus "python $pyscript -f {1} \
    #                                             --qscore $qscore" ::: ./results/split_fasta/sample_"$sample".*.fasta
    parallel -j $cpus "seqtk seq -F '\"$qscore\"' \
                        {1} > ./results/split_fasta/sample_\"$sample\".{#}.fastq" ::: ./results/split_fasta/sample_"$sample".*.fasta




    parallel -j $cpus "python $pyscript2 -f {1} \
                                                -k $kinetics" ::: ./results/split_fasta/sample_"$sample".*.fastq


    cat ./results/split_fasta/sample_"$sample".*.nrseq.fastq | pigz -p $cpus > "$output" 


fi


