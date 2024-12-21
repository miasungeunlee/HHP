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
chmod u+x *.sh # make the script executable
````

## Download the Refseq prokaryotic protein sequence databases
````
conda activate HHP
sh download_database.sh [-d database_directory]
````
### Example of the test run sh download_database.sh:
````
conda activate HHP
sh download_database.sh -d /store/ampere/slee/HHP_database # It takes 
````

•	````Refseq_prokaryotes_all_proteins.dmnd````: Refseq bacterial and archaeal proteins were downloaded from NCBI (https://www.ncbi.nlm.nih.gov/) and Diamond blastp database was created

•	````names.dmp & nodes.dmp````:

•	````prot.accession2taxid````:


## Quick run
````
conda activate HHP
sh HHP.sh [-h help] [-i fasta_file] [-d database_directory] [-o working_directory] [-t threads] 

### option variable explanation ###
-h: help
-i: virus genome fasta file path (e.g. /store/ampere/slee/HPP/test-genomes   # Full pathway
-d: database directory path (e.g. /store/ampere/slee/HHP_database   # Full pathway 
-o: working & output directory name (e.g. HHP_output)
-t: number of CPUs
````

### Example of the test run sh HHP.sh:
````
conda activate HHP
cd HHP # where the git clone is stored
sh HHP.sh -i /store/ampere/slee/HPP/test-genomes/virus-contig.fasta -d /store/ampere/slee/HHP_database -o HHP_output
````

### What does HHP.sh script do?

**Step 0)** Making AMOA databases ````AMO.dmnd```` using Diamond tool (Buchfink et al. 2021) in ````$exp_name```` working directory.

**Step 1)** Running the DADA2 pipeline ````dada_AMO.R```` to generate the amplicon sequence variant (ASV) sequences ````out.DADA2.$organism.ASVs.fa```` and ASV count table ````out.DADA2.$organism.ASVs.counts.tsv```` across different samples in $exp_name working directory.

**Step 2)** Select the AMOA ASV sequences according to expected amplicon size using Seqkit tool (Shen et al. 2016). Generating correct ASV sequences ````out.DADA2.correct-size.$organism.ASVs.fa```` and ASV count table ````out.DADA2.correct-size.$organism.ASVs.counts.tsv```` across different samples. 

**Step 3)** Correct AmoA sequence curation using the AMOA database (NR & IMG-JGI) with Diamond blastx (Buchfink et al. 2021) and manually curated conserved AMOA protein sequences using python scripts. Prior to correct translation, first nucleotide & two first nucleotides are removed from AOA & AOB ASV sequences. 
Generating annotated ASV sequences (nucleotide sequence, protein sequence) ````annotated.DADA2.{organism}.ASVs.fa````, ````annotated.DADA2.{organism}.ASVs.faa````, and AMOA-SEQ curated ASV sequences (nucleotide sequence, protein sequence) ````AMOA-SEQ-curated.$organism.ASVs.fa````,  ````AMOA-SEQ-curated.$organism.ASVs.faa````, and the ASV count table ````annotated.{organism}.ASVs.counts.tsv````, ````AMOA-SEQ-curated.{organism}.ASVs.counts.tsv```` across different samples. 

**Step 4)** Comparing the AMOA-SEQ curated ASV sequences to curated AMOA databases with Diamond blastx (Buchfink et al. 2021). Generating diamond blastx output files (every hits and best-hit): ````diamond.output.curateddb.AMOA-SEQ-curated.$organism.ASVs.tsv````, ````besthit.diamond.output.curateddb.AMOA-SEQ-curated.$organism.ASVs.tsv````


### Output directory and files
````
{organism}.ASV-analysis # directory
````
•	````out.DADA2.$organism.ASVs.track-summary.tsv````: DADA2 output summary containing quality control, denoising, number of merged sequences, number of chimeras in each sample. 

•	````out.DADA2.$organism.ASVs.fa````: generated amplicon sequence variants (however, this file could contain the ambiguous sequences, thus not recommended to directly use this ASV file).

•	````out.DADA2.$organism.ASVs.counts.tsv````: DADA2 ASV count table from different samples

•	````out.DADA2.correct-size.$organism.ASVs.fa````: selected ASVs according to expected amplicon size (option -t) 




