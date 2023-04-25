#!/bin/bash

# Usage: ./random_sample.sh input_file output_file num_samples

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input_file output_file num_samples"
    exit 1
fi

input_file="$1"
output_file="$2"
num_samples="$3"

# Check if input_file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' does not exist."
    exit 1
fi

# Check if num_samples is a positive integer
if ! [[ "$num_samples" =~ ^[0-9]+$ ]]; then
    echo "Error: Number of samples must be a positive integer."
    exit 1
fi

# Get the total number of lines in the input_file
total_lines=$(wc -l < "$input_file")

# Check if num_samples is not greater than total_lines - 1 (excluding header)
if [ "$num_samples" -ge "$total_lines" ]; then
    echo "Error: Number of samples cannot be greater than or equal to the total number of lines in the input file (excluding header)."
    exit 1
fi

# Save the header line to the output_file
head -n 1 "$input_file" > "$output_file"

# Sample random lines from input_file (excluding header) and append them to output_file
tail -n +2 "$input_file" | shuf -n "$num_samples" >> "$output_file"

echo "Sampled $num_samples lines from '$input_file' (excluding header) and saved them to '$output_file' (including header)."

