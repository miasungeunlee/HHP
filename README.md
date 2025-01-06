````
  ##################################
  ###    _   _   _   _   _____   ###  
  ###   | | | | | | | | |  _  |  ###   
  ###   | |_| | | |_| | | |_| |  ###    
  ###   |  _  | |  _  | |  ___|  ###  
  ###   | | | | | | | | | |      ###   
  ###   |_| |_| |_| |_| |_|      ###  
  ##################################
````
**HPP** Predicting host information from virus gene annotation

## Installation with Conda
````
conda create -y -n HHP -c bioconda kaiju prodigal diamond
conda activate HHP
pip install pandas
````

## Download HHP directory
````
git clone https://github.com/miasungeunlee/HHP.git
cd HHP
chmod u+x *.sh
````

## Download the Refseq prokaryotic protein sequence databases
````
conda activate HHP
sh download_database.sh [-d database_directory]
````
### Example of the test run sh download_database.sh:
````
conda activate HHP
sh download_database.sh -d /store/ampere/slee/HHP_database 
````

•	````Refseq_prokaryotes_all_proteins.dmnd````: Refseq bacterial and archaeal proteins were downloaded from NCBI (https://www.ncbi.nlm.nih.gov/) and Diamond blastp database was created

•	````names.dmp & nodes.dmp````:

•	````prot.accession2taxid````:


## Quick run
````
conda activate HHP
sh HHP.sh [-h help] [-i fasta_file] [-d database_directory] [-o output_directory] [-t threads] 

### option variable explanation ###
-h: help
-i: virus genome fasta file path (e.g. /store/ampere/slee/HPP/test-genomes/virus-contig.fasta   # Full pathway
-d: database directory path (e.g. /store/ampere/slee/HHP_database   # Full pathway 
-o: output_directory name (e.g. HHP_output)
-t: number of CPUs
````

### Example of the test run sh HHP.sh:
````
conda activate HHP
cd HHP # where the git clone is stored
sh HHP.sh -i /store/ampere/slee/HPP/test-genomes/virus-contig.fasta -d /store/ampere/slee/HHP_database -o HHP_output
````

### What does HHP.sh script do?

**Step 1)** Gene prediction using Prodigal with -p meta option (Hyatt et al. 2010) in ````$output_directory/GENE_PREDICTION````.

**Step 2)** Gene annotation using Diamond blastp against NCBI Refseq bacterial & archaeal proteins database in ````$output_directory/GENE_ANNOTATION```` .

**Step 3)** Getting taxid from prot_accession number using prot_accession2taxid file and adding the taxonomic rank using Kaiju (Menzel et al. 2016) in ````$output_directory/GENE_ANNOTATION````. 

**Step 4)** HHP pipeline for predicting the host from annotated genes. A homolog-based approach was performed as previously described (Al-Shayeb et al. 2020; Lee et al. 2022) with a host predicted as the phylum (family or genus) with ≥3x shared homologs compared to the second
most dominant phylum.

### Output directory and files
````
$output_directory/GENE_PREDICTION
````
•	````gene_aa_${base_name}.faa````: Predicted amino acid sequences of ORFs (Open Reading Frames)
•	````gene_nuc_${base_name}.fna````: Predicted nucleotide sequences of ORFs  
•	````Prodigal````: Prodigal output file (GenBank-like format with detailed annotations for coding sequences (CDS))

````
$output_directory/GENE_ANNOTATION
````
•	````nr.diamond.${base_name}.tsv````: Gene annotation against Refseq protein databases 
•	````best-hit-nr.diamond.${base_name}.tsv````: Top alignment for annotated genes
•	````kaiju.out````: Kaiju output format
•	````kaiju.names.out````: Taxonomic ranks are added (e.g.,superkingdom, phylum, order, class, family, genus, species)

$output_directory/HOST_PREDICTION
````
•	````Homologs-based-host-prediction-phylum.txt````: The host of viral contigs at the phylum level was predicted. Corresponding homolog numbers have been added to the file.
•	````Homologs-based-host-prediction-family.txt````: The host of viral contigs at the family level was predicted. Corresponding homolog numbers have been added to the file.
•	````Homologs-based-host-prediction-genus.txt````: The host of viral contigs at the genus level was predicted. Corresponding homolog numbers have been added to the file.
•	````HPP_host_prediction.txt````: A concatenated file of the phylum-, family-, and genus-level files.





