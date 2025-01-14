#!/bin/bash
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Download the Refseq Prokaryotic protein databases
# Date: 20241221

# Default values
database_directory=""

# Parse arguments
while getopts "d:h" opt; do
  case $opt in
    d) database_directory="$OPTARG" ;; # Capture the database_directory path
    h)
      echo "Usage: sh $0 -d /path/to/database_directory"
      exit 0
      ;;
    *)
      echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# Check if database_directory was provided
if [[ -z "$database_directory" ]]; then
  echo "Error: Output directory not provided. Use -d to specify it."
  exit 1
fi

# Convert to absolute path
database_directory=$(realpath "$database_directory")

# Ensure the output directory exists
mkdir -p "$database_directory"

# Change to the output directory
cd "$database_directory" || exit 1

echo "Downloading and processing files in $database_directory"

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
diamond makedb --in Refseq_nr_prokaryotes.faa --db Refseq_prokaryotes_all_proteins.dmnd --taxonmap prot.accession2taxid --taxonnodes nodes.dmp

echo "Download and processing complete. Files are saved in $database_directory"
