#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Counting homolog genes at genus or family or phylum level, run the script with the -i or --input argument followed by the desired input type (genus, family, or phylum)
# Date: 20240321

import os
import argparse
from collections import defaultdict

# Function to process each text file based on the argument
def process_file(file, field_index):
    counts = defaultdict(int)
    with open(file, 'r') as f:
        for line in f:
            fields = line.strip().split(';')
            if len(fields) > field_index:
                counts[fields[field_index]] += 1
    
    # Sort the counts
    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    
    # Return sorted counts for further processing
    return sorted_counts

# Main function
def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Process input files and create output files.')
    
    # Add argument for input type
    parser.add_argument('-i', '--input', type=str, required=True, choices=['family', 'genus', 'phylum'], help='Input type (family, genus, or phylum)')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Get the index corresponding to the input type
    if args.input == 'family':
        field_index = 4
        output_filename = 'family-count'
    elif args.input == 'genus':
        field_index = 5
        output_filename = 'genus-count'
    elif args.input == 'phylum':
        field_index = 1
        output_filename = 'phylum-count'
    
    # Get a list of all .txt files in the current directory
    files = [file for file in os.listdir() if file.endswith('.txt')]
    
    # Dictionary to store counts for each file
    file_counts = {}
    
    # Process each file
    for file in files:
        sorted_counts = process_file(file, field_index)
        file_counts[file] = sorted_counts
    
    # Write to output file
    with open(output_filename, 'w') as output_file:
        for file, counts in file_counts.items():
            # Write counts
            for count, value in counts:
                output_file.write(f"sort-{file}\t{count}\t{value}\n")

if __name__ == "__main__":
    main()
