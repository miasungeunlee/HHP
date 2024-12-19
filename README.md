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
````

## Installation with Conda
````
conda create -y -n HHP -c bioconda kaiju prodigal diamond
conda activate HHP
````

## Download HHP directory
````
git clone https://github.com/miasungeunlee/HHP.git
cd HHP
chmod u+x HHP.sh # make the script executable
````

## Refseq protein sequences databases
•	````Refseq_prokaryotes_all_proteins.dmnd````: Refseq bacterial and archaeal proteins were downloaded from NCBI (https://www.ncbi.nlm.nih.gov/) and Diamond blastp database was created

•	````names.dmp & nodes.dmp````:

•	````prot.accession2taxid````:

## Quick run
````
conda activate HHP
sh HHP.sh [-h help] [-i fasta_file] [-d database_directory] [-o working_directory] [-t threads] 

### option variable explanation ###
-h: help
-i: fastq.gz file path (e.g. /home/ampere/slee/COMICON-Projet-2022/   # Full pathway for running Rscript) 
-d: database directory path (e.g. /home/ampere/slee/COMICON-Projet-2022/   # Full pathway 
-o: working & output directory name (e.g. output)
-t: number of CPUs
````

### Example of the test run sh HHP.sh:
````
conda activate HHP
sh HHP.sh -i /store/ampere/gnicol/HPP/TEST-fasta/virus-contig.fasta -d /store/ampere/gnicol/test-database -o ouput_directory
````

### What does AMOA-SEQ.sh script do?

**Step 0)** Making AMOA databases ````AMO.dmnd```` using Diamond tool (Buchfink et al. 2021) in ````$exp_name```` working directory.

**Step 1)** Running the DADA2 pipeline ````dada_AMO.R```` to generate the amplicon sequence variant (ASV) sequences ````out.DADA2.$organism.ASVs.fa```` and ASV count table ````out.DADA2.$organism.ASVs.counts.tsv```` across different samples in $exp_name working directory.

**Step 2)** Select the AMOA ASV sequences according to expected amplicon size using Seqkit tool (Shen et al. 2016). Generating correct ASV sequences ````out.DADA2.correct-size.$organism.ASVs.fa```` and ASV count table ````out.DADA2.correct-size.$organism.ASVs.counts.tsv```` across different samples. 

**Step 3)** Correct AmoA sequence curation using the AMOA database (NR & IMG-JGI) with Diamond blastx (Buchfink et al. 2021) and manually curated conserved AMOA protein sequences using python scripts. Prior to correct translation, first nucleotide & two first nucleotides are removed from AOA & AOB ASV sequences. 
Generating annotated ASV sequences (nucleotide sequence, protein sequence) ````annotated.DADA2.{organism}.ASVs.fa````, ````annotated.DADA2.{organism}.ASVs.faa````, and AMOA-SEQ curated ASV sequences (nucleotide sequence, protein sequence) ````AMOA-SEQ-curated.$organism.ASVs.fa````,  ````AMOA-SEQ-curated.$organism.ASVs.faa````, and the ASV count table ````annotated.{organism}.ASVs.counts.tsv````, ````AMOA-SEQ-curated.{organism}.ASVs.counts.tsv```` across different samples. 

**Step 4)** Comparing the AMOA-SEQ curated ASV sequences to curated AMOA databases with Diamond blastx (Buchfink et al. 2021). Generating diamond blastx output files (every hits and best-hit): ````diamond.output.curateddb.AMOA-SEQ-curated.$organism.ASVs.tsv````, ````besthit.diamond.output.curateddb.AMOA-SEQ-curated.$organism.ASVs.tsv````

**Step 5)** Clustering the AMOA ASV sequences into OTUs ````AMOA-SEQ.$organism.OTUs.fa```` with 97% of sequence identity using CDHIT tool (Li et al. 2006) and generating OTU count table across different samples ````AMOA-SEQ.$organism.OTUs.counts.tsv```` and annotating OTUs with curated AMOA databases with Diamond blastx (Buchfink et al. 2021). Generating diamond blastx output files (every hits and best-hit): ````diamond.output.curateddb.AMOA-SEQ.$organism.OTUs.tsv````, ````besthit.diamond.output.curateddb.AMOA-SEQ.$organism.OTUs.tsv````

**Step 6)** Translating the AMOA-SEQ curated ASV sequences to protein sequence variant (PSV) sequences ````AMOA-SEQ.$organism.PSVs.faa```` using Seqkit tool (Shen et al. 2016). Dereplicating the PSV sequences ````{organism}.PSV.faa````using CDHIT tool (Li et al. 2006)

**Step 7)** Annotating the PSV sequences against curated AMOA database using BLASTp. Generating diamond blastx output files (every hits and best-hit): ````blastp.output.AMOA-SEQ.$organism.PSVs.tsv````, ````besthit.blastp.output.AMOA-SEQ.$organism.PSVs.tsv````

**Step 8)** Aligning of the PSV sequences ````AMOA-SEQ.{organism}.PSVs.faa```` and curated AMOA sequences ````ref.{organism}.amoA.faa```` using MUSCLE (Edgar et al. 2004) and spurious sequences or poorly aligned regions were removed using trimAI (Capella-Gutiérrez · 2009). ````tree.$organism.trim.afa```` is used for generating phylogenetic tree ````tree.{organism}.nwk```` using FastTree (Price et al. 2009) and IQTree (Nguyen et al. 2015).



### Output directory and files
````
{organism}.ASV-analysis # directory
````
•	````out.DADA2.$organism.ASVs.track-summary.tsv````: DADA2 output summary containing quality control, denoising, number of merged sequences, number of chimeras in each sample. 

•	````out.DADA2.$organism.ASVs.fa````: generated amplicon sequence variants (however, this file could contain the ambiguous sequences, thus not recommended to directly use this ASV file).

•	````out.DADA2.$organism.ASVs.counts.tsv````: DADA2 ASV count table from different samples

•	````out.DADA2.correct-size.$organism.ASVs.fa````: selected ASVs according to expected amplicon size (option -t) 

•	````out.DADA2.correct-size.$organism.ASVs.counts.tsv````: correct size ASV count table from different samples

•	````annotated.DADA2.$organism.ASVs.fa````: DADA2 ASVs with correct size and matched to AMOA database. 

•	````annotated.DADA2.$organism.ASVs.faa````: translated DADA2 ASVs with correct size and matched to AMOA database. 

•	````AMOA-SEQ-curated.$organism.ASVs.fa````: Ambiguous sequences were removed using using the AMO databases and conserved protein sequence (those are confident and genuine AMOA sequences).

•	````AMOA-SEQ-curated.$organism.ASVs.faa````: translated AMOA-SEQ curated ASVs. 

•	````AMOA-SEQ-curated.$organism.ASVs.counts.tsv````: AMOA-SEQ curated ASV count table from different samples (recommanded to use this file for further bioinformatic analysis). 

•	````diamond.output.DADA2.$organism.ASVs.tsv````: annotation of DADA2 ASVs using total AMOA database

•	````besthit.diamond.output.DADA2.$organism.ASVs.tsv````: besthit of annoated ASVs using total AMOA database

•	````diamond.output.curateddb.AMOA-SEQ-curated.$organism.ASVs.tsv````: annotation of AMOA-SEQ curated ASVs using curated AMOA database  

•	````besthit.diamond.output.curateddb.AMOA-SEQ-curated.$organism.ASVs.tsv````: besthit of AMOA-SEQ curated ASVs using curated AMOA database  



````
{organism}.OTU-analysis # directory
````
•	````AMOA-SEQ.{organism}.OTUs.fa.clstr````: cd-hit clustering output file

•	````AMOA-SEQ.{organism}.OTUs.fa````: generated OTU sequences

•	````AMOA-SEQ.{organism}.OTUs.counts.tsv````: OTU count table from different samples

•	````diamond.output.curateddb.AMOA-SEQ.{organism}.OTUs.tsv````: annotation of OTUs using curated AMOA database  

•	````besthit.diamond.output.curateddb.AMOA-SEQ.{organism}.OTUs.tsv````: besthit of annoated OTUs using curated AMOA database  

•	````AMOA-SEQ.{organism}.OTUs.taxa.tsv````: OTU taxa information 


````
{organism}.PSV-analysis # directory
````

•	````AMOA-SEQ.{organism}.PSVs.faa.clstr````: clustering of translated ASV sequences with 100% identity

•	````AMOA-SEQ.{organism}.PSVs.faa````: Unique protein sequence variants

•	````blastp.output.AMOA-SEQ.{organism}.PSVs.tsv````: (recommended for beta-diversity analysis (e.g. phylogenetic tree)): annotation of PSVs using curated AMOA database

•	````besthit.blastp.output.AMOA-SEQ.{organism}.PSVs.tsv````: besthit of annotated PSVs using curated AMOA database 

````
{organism}.phylogeny-analysis # directory
````

•	````tree.{organism}.faa````: PSV sequences + curated AMOA sequence for phylogeny analysis

•	````tree.{organism}.afa````: AMOA-aligments for each PSV and curated AMOA sequences using MUSCLE with -super5 option

•	````tree.{organism}.trim.afa````: Ambiguous regions and gaps removed using trimal with -nogaps option

•	````FastTree.{organism}.nwk````: AMOA phylogenetic tree generated using FastTREE file in Newick format  

•	````IQTree.{organism}.treefile````: AMOA phylogenetic tree generated using IQTree file in Newick format 



