#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Making an individual file for each virus contig
# Date: 20240321

import argparse

# Create argument parser
parser = argparse.ArgumentParser(description='Process input file and create output files.')

# Add argument for input file
parser.add_argument('-i', '--input', type=str, required=True, help='Input file name')

# Parse the arguments
args = parser.parse_args()

# Input file provided through command-line argument
input_file = args.input

# Open the input file
with open(input_file, "r") as f:
    # Iterate over each line in the input file
    for line in f:
        # Split the line into columns
        columns = line.strip().split()
        # Extract the unique name (first column)
        unique_name = columns[0]
        # Define the output file name
        output_file = unique_name + ".txt"
        # Write the line to the output file
        with open(output_file, "a") as output:
            output.write(line)
