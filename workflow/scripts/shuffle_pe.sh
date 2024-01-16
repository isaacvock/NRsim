#!/bin/sh

# Adapted from: https://gist.github.com/iam28th/418dc7d5048067af194a76ffb5840c90
# Discussed in: https://www.biostars.org/p/9764/

input1=$1
input2=$2
temp1=$3
temp2=$4
output1=$5
output2=$6


parallel --link -j 2 "gzip -d -c {1} > {2} " ::: "$input1" "$input2" \
                                                    ::: "$temp1" "$temp2" 

paste -d '\n' "$temp1" "$temp2" | \
    awk '{

    lines[1] = $0;
    for (i = 2; i <= 8; ++i)
        getline lines[i];

    for (i = 1; i <= 8; i++)
        printf("%s%s", lines[i], i == 8 ? "'"\n"'" : "'"\t"'")

    }' | \

    shuf | \

    tr '\t' '\n' | \

    awk -v file1=$output1 -v file2=$output2 \
        'NR % 2 == 1 { print | "pigz > " file1 }
         NR % 2 == 0 { print | "pigz > " file2 }
