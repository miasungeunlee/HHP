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
conda create -y -n AMOA-SEQ -c bioconda seqkit fasttree muscle blast diamond trimal diamond=0.9.19 cd-hit cutadapt iqtree
source activate AMOA-SEQ
pip install biopython
pip install pandas
````

## DADA2 tool installation 
````
https://benjjneb.github.io/dada2/dada-installation.html
````

## Download AMOA-SEQ directory
````
git clone https://github.com/miasungeunlee/AMOA-SEQ.git
cd AMOA-SEQ
chmod u+x AMOA-SEQ.sh # make the script executable
````

## AMOA sequence databases
•	````AMO_database.faa````: all AMOA sequences downloaded from JGI IMG (https://img.jgi.doe.gov/) and NCBI (https://www.ncbi.nlm.nih.gov/). List of accessions and detail of the gene set were shown in AMO_database.tsv 

•	````ref.AOA.amoA.faa````: curated archaeal amoA sequences with defined lineage from Alves et al. 2018

•	````ref.COM.amoA.faa````: curated archaeal amoA sequences from Palomo et al. 2022

•	````ref.AOB.amoA.faa````: curated Nitrosospira amoA sequences from Aigle et al. 2019 and Nitrosomonas amoA sequences were downloaded from JGI IMG site and curated (Lee et al. 2023) 

## Quick run
````
source activate AMOA-SEQ
sh AMOA-SEQ.sh [-h help] [-e output_directory_name] [-i fastq_directory] [-f forward_primer] [-r reverse_primer] [-m minimum_read_length] [-l truncation_read_length] [-c just_concatenating_option] [-t expected merged sequence_length] [-t expected_merged_sequence_length] [-n number_nucleotide] [-o AO_type]

### option variable explanation ###
-h: help
-e: output directory name (e.g. output)
-i: fastq.gz file path (e.g. /home/ampere/slee/COMICON-Projet-2022/   # Full pathway for running Rscript) 
# fastq file must end with either "_R1_001.fastq.gz" or "_R2_001.fastq.gz" pattern (directly from MiSeq sequencing) #
-f: forward primer (e.g. ATGGTCTGGCTWAGACG for AOA forward primer, GGGGTTTCTACTGGTGGT for AOB forward primer and AGGNGAYTGGGAYTTCTGG for COMAMMOX forword primer)
-r: reverse primer (e.g. GCCATCCATCTGTATGTCCA for AOA reverse primer, CCCCTCKGSAAAGCCTTCTTC for AOB reverse primer and CGGACAWABRTGAABCCCAT for COMAMMOX reverse primer)
-m: minimum read length (e.g. 200 bp for AOA, 232 bp for AOB, 204 bp for COMMAMMOX in order to overlap between R1 and R2)
-l: truncation read length (e.g. 200 bp for AOA, 232 bp for AOB, 204 bp for COMMAMMOX)
-c: TRUE for AOA AMOA, just concatenating forward and reverse reads, for other AMOA, merging forward and reverse reads are possible, use FALSE option
-t: expected merged sequence length (410 bp for AOA, 452 bp for AOB and 396 bp for COMAMMOX; these lengths are after removing the primer lengths)
-n: number of nucleotides to be removed prior to correct translation (first nucleotide & two first nucleotides were removed from AOA & AOB ASV sequences)
-o: organism; It can be either AOA, AOB, or COM, depending on your dataset.
````

### Example of the test run AMOA-SEQ.sh for AOA, AOB and Comammox AMOA amplicon sequencing:
````
source activate AMOA-SEQ
# For Miseq reagent kit V3, 2 x 300 bp
sh AMOA-SEQ.sh -e AOA-output -i /home/ampere/slee/AMOA-SEQ/TEST-AOA-Fastq -f ATGGTCTGGCTWAGACG -r GCCATCCATCTGTATGTCCA -m 200 -l 200 -c TRUE -t 410 -n 2 -o AOA
sh AMOA-SEQ.sh -e AOB-output -i /home/ampere/slee/AMOA-SEQ/TEST-AOB-Fastq -f GGGGTTTCTACTGGTGGT -r CCCCTCKGSAAAGCCTTCTTC -m 232 -l 232 -c FALSE -t 452 -n 3 -o AOB
sh AMOA-SEQ.sh -e COM-output -i /home/ampere/slee/AMOA-SEQ/TEST-COM-Fastq -f AGGNGAYTGGGAYTTCTGG -r CGGACAWABRTGAABCCCAT -m 204 -l 204 -c FALSE -t 396 -n 1 -o COM

# For Miseq reagent kit V2, 2 x 250 bp
sh AMOA-SEQ.sh -e AOA-output -i /home/ampere/slee/AMOA-SEQ/TEST-AOA-Fastq -f ATGGTCTGGCTWAGACG -r GCCATCCATCTGTATGTCCA -m 200 -l 200 -c TRUE -t 410 -n 2 -o AOA
sh AMOA-SEQ.sh -e AOB-output -i /home/ampere/slee/AMOA-SEQ/TEST-AOB-Fastq -f GGGGTTTCTACTGGTGGT -r CCCCTCKGSAAAGCCTTCTTC -m 200 -l 200 -c TRUE -t 410 -n 3 -o AOB
sh AMOA-SEQ.sh -e COM-output -i /home/ampere/slee/AMOA-SEQ/TEST-COM-Fastq -f AGGNGAYTGGGAYTTCTGG -r CGGACAWABRTGAABCCCAT -m 204 -l 204 -c FALSE -t 396 -n 1 -o COM
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



