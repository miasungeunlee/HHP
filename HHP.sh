#!/bin/bash
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Homolog-based Host Prediction tool bash script for analysing virus genomes
# Date: 20240321

# Set default values
input_fasta=""
data_directory=""
working_directory=""
threads=""
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
echo "Usage: $0 [-i INPUT_FASTA] [-d DATA_DIRECTORY] [-w WORKING_DIRECTORY] [-t THREADS] [-h]"
}

# Parse command-line options
while getopts "i:d:w:t:h" flag; do
case "${flag}" in
i) input_fasta="${OPTARG}" ;;
d) data_directory="${OPTARG}" ;;
w) working_directory="${OPTARG}" ;;
t) threads="${OPTARG}" ;;
h) usage; exit 0 ;;
*) usage; exit 1 ;;
esac
done

while getopts "i:d:w:t" flag; do
    case "${flag}" in
        i) input_fasta="${OPTARG}" ;;
        d) data_directory="${OPTARG}" ;;
        w) working_directory="${OPTARG}" ;;
        t) threads="${OPTARG}" ;;
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
echo "#Database directory: $data_directory";
echo "#working & output directory: $working_directory";
echo "============================================================================================";
protein_file="Refseq_nr_prokaryotes.faa"
accession_file="prot.accession2taxid"
diamond_db_file="Refseq_prokaryotes_all_proteins.dmnd"

#####################################################
echo "============================================================================================";
echo "### STEP 0. Moving files into $working_directory working directory ###"
echo "============================================================================================";
mkdir $working_directory
# Convert to full path
working_directory=$(realpath "$working_directory")
# Ensure the output directory exists
mkdir -p "$working_directory"
cp *.py $working_directory
cp $input_fasta $working_directory
cd $working_directory
echo "### STEP 0. Done ###"
echo "============================================================================================";
#####################################################


#####################################################
echo "============================================================================================";
echo "### STEP 1. Gene prediction ###"


### HOMOLOGUE-HOST-PREDICTION ###
base_name=$(basename "$input_fasta" .fasta)
# Run Prodigal with the modified output filenames
prodigal -i "${base_name}.fasta" \
         -a "gene_aa_${base_name}.faa" \
         -d "gene_nuc_${base_name}.fna" \
         -o "Prodigal" \
         -p meta
rm ${base_name}.fasta
echo "### STEP 1. Gene prediction done ###"
echo "============================================================================================";
#####################################################


#####################################################
echo "============================================================================================";
echo "### STEP 2. Gene annotation ###"

diamond blastp --db $data_directory/$diamond_db_file --query gene_aa_${base_name}.faa --out nr.diamond.${base_name}.tsv --evalue 0.00001  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen qcovhsp --threads $threads
awk '!x[$1]++' nr.diamond.${base_name}.tsv > best-hit-nr.diamond.${base_name}.tsv # best-hit #

mkdir GENE_PREDICTION
mv gene_aa_${base_name}.faa ./GENE_PREDICTION
mv gene_nuc_${base_name}.fna ./GENE_PREDICTION
#mv Prodigal ./GENE_PREDICTION

echo "### STEP 2. Gene annotation done ###"
echo "============================================================================================";
#####################################################

#####################################################
echo "============================================================================================";
echo "### STEP 3. Adding the taxonomic rank using the prot_accession2taxid and Kaiju ###"
cut -f1,2,12 best-hit-nr.diamond.${base_name}.tsv > prot_ids.tsv
python $working_directory/prot.accession2taxid.py -t $data_directory/prot.accession2taxid -i prot_ids.tsv -o taxids.tsv
mkdir GENE_ANNOTATION
mv nr.diamond.${base_name}.tsv ./GENE_ANNOTATION
mv best-hit-nr.diamond.${base_name}.tsv ./GENE_ANNOTATION

cat taxids.tsv > kaiju.out
rm taxids.tsv
cut -d ' ' -f1,2 kaiju.out > tmp && mv tmp kaiju.out
awk '{print "C "$0}' kaiju.out > tmp && mv tmp kaiju.out
sed 's/ /\t/g' kaiju.out > tmp && mv tmp kaiju.out
kaiju-addTaxonNames -t $data_directory/nodes.dmp -n $data_directory/names.dmp -i kaiju.out -o kaiju.names.out -r superkingdom,phylum,order,class,family,genus,species
mv kaiju.out ./GENE_ANNOTATION
### Add virus-ID ###
awk 'BEGIN {FS="\t"}; {print $2}' kaiju.names.out > ID.txt
sed 's/_[0-9]*$//' ID.txt > tmp && mv tmp ID.txt
paste ID.txt kaiju.names.out > Virus-ID-kaiju.names.out
mv kaiju.names.out ./GENE_ANNOTATION
rm ID.txt
echo "### STEP 3. Done ###"
echo "============================================================================================";

echo "============================================================================================";
echo "### STEP 4. Host prediction ###"
mkdir tmps
mv Virus-ID-kaiju.names.out ./tmps
cd tmps
### Make a individual file for each virus contig ###
python $working_directory/indiv.py -i Virus-ID-kaiju.names.out
python $working_directory/process_files.py -i phylum
python $working_directory/process_files.py -i family
python $working_directory/process_files.py -i genus
awk -F '\t' '{print $1 "\t" $3 "\t" $2}' phylum-count > tmp && mv tmp phylum-count
awk -F '\t' '{print $1 "\t" $3 "\t" $2}' genus-count > tmp && mv tmp genus-count
awk -F '\t' '{print $1 "\t" $3 "\t" $2}' family-count > tmp && mv tmp family-count

mkdir tmps-phylum
mkdir tmps-family
mkdir tmps-genus

mv phylum-count ./tmps-phylum
mv family-count ./tmps-family
mv genus-count ./tmps-genus

### extract two best hits and check three times or not (Banfield's lab threshold) ###
cd tmps-phylum
awk '{ if (++count[$1] <= 2) print $0 }' phylum-count > best-two-phylum-count
python $working_directory/indiv.py -i best-two-phylum-count
find . -maxdepth 1 -type f -name 'sort*.txt' | sort > Name-list.txt
sed 's/.\///g' Name-list.txt > tmp && mv tmp Name-list.txt
python $working_directory/selection.py
cat New_* > Homologs-based-host-prediction-phylum.txt
sed 's/sort-//g' Homologs-based-host-prediction-phylum.txt > tmp && mv tmp Homologs-based-host-prediction-phylum.txt
sed 's/\.txt//g' Homologs-based-host-prediction-phylum.txt > tmp && mv tmp Homologs-based-host-prediction-phylum.txt
rm New_* sort-*
cp Homologs-based-host-prediction-phylum.txt $working_directory
cd ..

cd tmps-family
awk '{ if (++count[$1] <= 2) print $0 }' family-count > best-two-family-count
python $working_directory/indiv.py -i best-two-family-count
find . -maxdepth 1 -type f -name 'sort*.txt' | sort > Name-list.txt
sed 's/.\///g' Name-list.txt > tmp && mv tmp Name-list.txt
python $working_directory/selection.py
cat New_* > Homologs-based-host-prediction-family.txt
sed 's/sort-//g' Homologs-based-host-prediction-family.txt > tmp && mv tmp Homologs-based-host-prediction-family.txt
sed 's/\.txt//g' Homologs-based-host-prediction-family.txt > tmp && mv tmp Homologs-based-host-prediction-family.txt
rm New_* sort-*
cp Homologs-based-host-prediction-family.txt $working_directory
cd ..

cd tmps-genus
### extract two best hits and check three times or not (Banfield's lab threshold) ###
awk '{ if (++count[$1] <= 2) print $0 }' genus-count > best-two-genus-count
python $working_directory/indiv.py -i best-two-genus-count
find . -maxdepth 1 -type f -name 'sort*.txt' | sort > Name-list.txt
sed 's/.\///g' Name-list.txt > tmp && mv tmp Name-list.txt
python $working_directory/selection.py
cat New_* > Homologs-based-host-prediction-genus.txt
sed 's/sort-//g' Homologs-based-host-prediction-genus.txt > tmp && mv tmp Homologs-based-host-prediction-genus.txt
sed 's/\.txt//g' Homologs-based-host-prediction-genus.txt > tmp && mv tmp Homologs-based-host-prediction-genus.txt
rm New_* sort-*
cp Homologs-based-host-prediction-genus.txt $working_directory
cp Homologs-based-host-prediction-genus.txt $working_directory
cd $working_directory
python $working_directory/concatenating_files.py
mkdir HOST_PREDICTION
mv Homologs-based-host-prediction* ./HOST_PREDICTION
mv HHP_host_prediction.tsv ./HOST_PREDICTION

echo "### STEP 4. Host prediction done ###"
echo "### HHP Pipeline successfully finished :) ###"
echo "============================================================================================";
