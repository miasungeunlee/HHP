#!/bin/bash
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Homolog-based Host Prediction tool bash script for analysing virus genomes
# Date: 20240321

# Set default values
input_fasta=""
data_directory=""
working_directory=""
function usage {
  echo "
  ##################################
  ###    _   _   _   _   _____   ###
  ###   | | | | | | | | |  _  |  ###
  ###   | |_| | | |_| | | |_| |  ###
  ###   |  _  | |  _  | |  ___|  ###
  ###   | | | | | | | | | |      ###
  ###   |_| |_| |_| |_| |_|      ###
  ##################################

  ";
echo "Homolog-based Host Prediction tool";
echo "Predicting host information from virus gene annotation";
echo " ";
echo "Usage: $0 [-i INPUT_FASTA] [-d DATA_DIRECTORY] [-w WORKING_DIRECTORY] [-h]"
}

# Parse command-line options
while getopts "i:d:w:h" flag; do
case "${flag}" in
i) input_fasta="${OPTARG}" ;;
d) data_directory="${OPTARG}" ;;
w) working_directory="${OPTARG}" ;;
h) usage; exit 0 ;;
*) usage; exit 1 ;;
esac
done

while getopts "i:d:w:" flag; do
    case "${flag}" in
        i) input_fasta="${OPTARG}" ;;
        d) data_directory="${OPTARG}" ;;
        w) working_directory="${OPTARG}" ;;
    esac
done


for opt in "${required_opts[@]}"; do
if [[ -z "${!opt}" ]]; then
echo "ERROR: ${opt} is required."
usage
exit 1
fi
done

echo "============================================================================================";
echo "#Host prediction of $input_fasta file";
echo "#Data_directory: $data_directory";
data_directory=$(pwd)
echo "#workng_directory: $workng_directory";
workng_directory=$(pwd)

echo "============================================================================================";

# Directory where the data is to be downloaded
# data_directory="/path/to/data_directory"
# Filenames to check for existence
protein_file="Refseq_nr_prokaryotes.faa"
accession_file="prot.accession2taxid"
taxdump_files=("nodes.dmp" "names.dmp")
diamond_db_file="Refseq_nr_prokaryotes.dmnd"

# Change to the data directory
cd "$data_directory"

# Check if the necessary files exist
if [[ -f "$protein_file" && -f "$accession_file" && -f "nodes.dmp" && -f "names.dmp" && -f "$diamond_db_file" ]]; then
    echo "Required database files and Diamond database already present in the Data_directory: $data_directory. Skipping download and database creation steps."
else
    if [[ -f "$protein_file" && -f "$accession_file" && -f "nodes.dmp" && -f "names.dmp" ]]; then
        echo "Required database files already present in the Data_directory: $data_directory. Skipping download steps."
    else
        echo "Downloading the database in the Data_directory: $data_directory"

        # Download the protein files
        wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/*protein.faa.gz
        wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/*protein.faa.gz

        # Unzip the downloaded files
        gunzip *.gz

        # Concatenate all protein files into one
        cat *.faa > "$protein_file"

        # Remove the individual protein files
        rm *protein.faa

        # Download additional database files
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

        # Extract the taxonomy dump
        tar -xvf taxdump.tar.gz

        echo "Databases successfully downloaded in the Data_directory: $data_directory"
    fi

    # Make the diamond database if it does not exist
    if [[ ! -f "$diamond_db_file" ]]; then
        echo "Making the nr database: $data_directory"
        diamond makedb --in "$data_directory/$protein_file" --db "$data_directory/$diamond_db_file" --taxonmap "$data_directory/$accession_file" --taxonnodes "$data_directory/nodes.dmp" --taxonnames "$data_directory/names.dmp"
        echo "Diamond nr database created: $data_directory"
    else
        echo "Diamond database already present: $data_directory. Skipping database creation step."
    fi
fi
echo "============================================================================================";


### HOMOLOGUE-HOST-PREDICTION ###
cd $workng_directory
base_name="${input_fasta%.fasta}"
# Run Prodigal with the modified output filenames
prodigal -i "$working_directory/$input_fasta" \
         -a "$working_directory/gene_aa_${base_name}.faa" \
         -d "$working_directory/gene_nuc_${base_name}.fna" \
         -o "$working_directory/Prodigal" \
         -p meta

diamond blastp -d "$data_directory/$diamond_db_file" -q $input_fasta -o nr.diamond.${base_name}.tsv --evalue 0.00001  --outfmt 6 qseqid sseqid sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen qcovhsp --threads 16

awk '!x[$1]++' nr.diamond.${base_name}.tsv > best-hit-nr.diamond.${base_name}.tsv # best-hit #
### Find name of nr database using prot_accession2taxid and Kaiju ###
cut -f1,2,12 best-hit-$diamond_file > prot_ids.tsv
python $workng_directory/prot.accession2taxid.py -t $data_directory/prot.accession2taxid -i prot_ids.tsv -o taxids.tsv
#### Add the taxonomy ranks using Kaiju tool ###
cat taxids.tsv > kaiju.out
cut -d ' ' -f1,2 kaiju.out > tmp && mv tmp kaiju.out
awk '{print "C "$0}' kaiju.out > tmp && mv tmp kaiju.out
sed 's/ /\t/g' kaiju.out > tmp && mv tmp kaiju.out
kaiju-addTaxonNames -t $data_directory/nodes.dmp -n $data_directory/names.dmp -i kaiju.out -o kaiju.names.out -r superkingdom,phylum,order,class,family,genus,species

### Add virus-ID ###
awk 'BEGIN {FS="\t"}; {print $2}' kaiju.names.out > ID.txt
sed 's/_[0-9]*$//' ID.txt > tmp && mv tmp ID.txt
paste ID.txt kaiju.names.out > Virus-ID-kaiju.names.out

### Separating file into individual mVC matching file ###
mkdir tmps
mv Virus-ID-kaiju.names.out ./tmps
cd tmps
### Make a individual file for each virus contig ###
python indiv.py -i Virus-ID-kaiju.names.out
python process_files.py -i phylum
python process_files.py -i family
python process_files.py -i genus

mkdir tmps-phylum
mkdir tmps-family
mkdir tmps-genus

mv phylum-count ./tmps-phylum
mv family-count ./tmps-family
mv genus-count ./tmps-genus

### extract two best hits and check three times or not (Banfield's lab threshold) ###
cd tmps-phylum
awk '{ if (++count[$1] <= 2) print $0 }' phylum-count > best-two-phylum-count
python indiv.py -i best-two-phylum-count
ls sort*.txt > Name-list.txt
python selection.py
cat New_* > Final-homologs-based-host-prediction.txt
sed 's/sort-//g' Homologs-based-host-prediction-phylum.txt | sed 's/.txt//g' > tmp && mv tmp Homologs-based-host-prediction-phylum.txt
rm New_* sort-*
cd ..

cd tmps-family
awk '{ if (++count[$1] <= 2) print $0 }' family-count > best-two-family-count
python indiv.py -i best-two-family-count
ls sort*.txt > Name-list.txt

cat New_* > Homologs-based-host-prediction-family.txt
rm New_* sort-*
cd ..

cd tmps-genus
### extract two best hits and check three times or not (Banfield's lab threshold) ###
awk '{ if (++count[$1] <= 2) print $0 }' genus-count > best-two-genus-count
python indiv.py -i best-two-family-count
ls sort*.txt > Name-list.txt
cat New_* > Homologs-based-host-prediction-genus.txt
rm New_* sort-*
cd ..
