#!/bin/sh

# Adapted from: https://gist.github.com/iam28th/49a245427ea2b8ed5f1f9889c13468bf
# Discussed in: https://www.biostars.org/p/9764/

input=$1
tempf=$2
output=$3
cpus=$4

gzip -d -c "$input" > "$tempf"


awk '{

# read 4 lines
lines[1] = $0;
for (i = 2; i <= 4; ++i)
    getline lines[i];

# and print them tab-separated on a single line
for (i = 1; i <= 4; ++i)
    printf("%s%s", lines[i], i == 4 ? "'"\n"'" : "'"\t"'")
}' "$tempf" | \

# shuffle
shuf | \

# replace all tabs back to newlines
tr '\t' '\n' | pigz -c -p "$cpus" > "$output" 