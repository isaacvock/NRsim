#!/bin/sh

input1=$1
input2=$2
temp1=$3
temp2=$4
output1=$5
output2=$6

parallel gzip -c -d "$input1" > "$temp1"
gzip -c -d "$input1" > "$temp1"

parallel -j $cpus --compress --plus "pigz -d -k -c -p 1 {1} " ::: $directory/*_counts*

paste -d '\n' $1 $2 | \
    awk '{

    lines[1] = $0;
    for (i = 2; i <= 8; ++i)
        getline lines[i];

    for (i = 1; i <= 8; i++)
        printf("%s%s", lines[i], i == 8 ? "'"\n"'" : "'"\t"'")

    }' | \

    shuf | \

    tr '\t' '\n' | \

    awk 'NR % 2 == 1 { print >> "shuffled_r1.fastq" }
         NR % 2 == 0 { print >> "shuffled_r2.fastq" }'
