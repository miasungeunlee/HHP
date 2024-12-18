#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Getting the taxID from the accession number
# Date: 20240321

import argparse

def main(accession_file_path, input_file_path, output_file_path):
    # Read accession file and create a mapping of accessions to taxids
    with open(accession_file_path, 'r') as accession_file:
        taxid_map = {}
        for line in accession_file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                taxid_map[parts[0]] = parts[2]

    # Process input file and write output to the specified file
    with open(input_file_path, 'r') as input_file, \
         open(output_file_path, 'w') as output_file:
        for line in input_file:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                accession = parts[1].split('.')[0]
                taxid = taxid_map.get(accession)
                if taxid:
                    output_file.write(f"{parts[0]}\t{taxid}\t{parts[2]}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process accession and taxid mapping files.')
    parser.add_argument('-t', '--accession_file', type=str, help='Path to the accession file')
    parser.add_argument('-i', '--input_file', type=str, help='Path to the input file')
    parser.add_argument('-o', '--output_file', type=str, help='Path to the output file')
    args = parser.parse_args()

    main(args.accession_file, args.input_file, args.output_file)
