#!/bin/bash

# Path to the directory containing genome files
genome_dir="/home/anna/storage/arabidopsis/pacbio/pb_27/"
# Path to the input file
input_file="/home/anna/storage/arabidopsis/split_read/scp.fasta"
# Path to the output directory
output_dir="/home/anna/storage/arabidopsis/split_read"

# Iterate over all .fasta files in the genome_dir
for genome_file in "$genome_dir"/*.fasta; do
  # Get the base name of the file without path and extension
  basename=$(basename "$genome_file" .fasta)
  # Form the output path for the current file
  output_path="$output_dir/$basename"
  # Execute the command
  ./sim_in_genome.sh -keepblast -sim 90 -in "$input_file" -genome "$genome_file" -out "$output_path"
done
