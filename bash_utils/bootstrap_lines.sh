#!/bin/bash

# parse command-line arguments
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <input_file> <output_file> <bootstrap_lines>"
    exit 1
fi

input_file=$1
output_file=$2
bootstrap_lines=$3

# get the number of lines in the input file
num_lines=$(wc -l < $input_file)

# generate a random sample of line numbers to bootstrap
echo "1" >> line_numbers.txt
shuf -i 2-$num_lines -n $bootstrap_lines > line_numbers.txt

# bootstrap the selected lines from the input file and output to a new file
awk 'NR==FNR{array[$0];next} FNR in array' line_numbers.txt $input_file > $output_file

# clean up
rm line_numbers.txt
