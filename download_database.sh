#!/bin/bash
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Download the Refseq Prokaryotic protein databases
# Date: 20241221

# Default values
output_directory=""

# Parse arguments
while getopts "o:h" opt; do
  case $opt in
    o) output_directory="$OPTARG" ;; # Capture the output-directory path
    h)
      echo "Usage: sh $0 -o /path/to/output-directory"
      exit 0
      ;;
    *)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# Check if output_directory was provided
if [[ -z "$output_directory" ]]; then
  echo "Error: Output directory not provided. Use -o to specify it."
  exit 1
fi

# Convert to absolute path
output_directory=$(realpath "$output_directory")

# Ensure the output directory exists
mkdir -p "$output_directory"

# Change to the output directory
cd "$output_directory" || exit 1

echo "Downloading and processing files in $output_directory"

# Download the protein files
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/*protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/*protein.faa.gz

# Unzip the downloaded files
gunzip *.gz

# Concatenate all protein files into one
cat *.faa > Refseq_nr_prokaryotes.faa

# Remove the individual protein files
rm *protein.faa

# Download additional database files
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
gunzip prot.accession2taxid.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# Extract the taxonomy dump
tar -xvf taxdump.tar.gz
diamond makedb --in Refseq_nr_prokaryotes.faa --db Refseq_prokaryotes_all_proteins.dmnd --taxonmap prot.accession2taxid --taxonnodes nodes.dmp --taxonnames names.dmp

echo "Download and processing complete. Files are saved in $output_directory"
