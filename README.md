# Phylodyanmic-Subsampling


# Repository description:
This code and data were used for the analysis presented in “Using multiple sampling strategies to estimate SARS-CoV-2 epidemiological parameters from genomic sequencing data" - https://www.medrxiv.org/content/10.1101/2022.02.04.22270165v1.

The repository contains the following elements:

1. Data

Contains some of the data (genomic data has been excluded due to permissions) used for analysis that is:

* Case data for Amazonas State, Brazil and Hong Kong
* Example XML files for Hong Kong and Amazonas State
* Accession ID's for genomic data

2. Code

Contains the scripts needed for processing and plotting using R. 

# Running the code

In this example we show how to run the process for the Hong Kong Dataset with proportional sampling.

Step 1: Obtainting and naming genomic data

* Download sequences and metadata from GISAID for Hong Kong up to 7th May 2020
* Run script Code/creation_of_EPI_summary.R to generate table for naming of sequences
* Run script Code/proccess_covid-19_sequences_for_beast_from_GISAID.R to correctly name sequences

Step 2: Process genomic data

* Align sequences using MAFFT 
* Remove first 130 and last 50 base pairs (Mega X was used to do this)

Step 3: Quality control using NexClade pipeline

nextclade \
--in-order \
--input-fasta data/sars-cov-2/sequences.aln.fasta \
--output-tsv output/nextclade.tsv \
--input-dataset data/sars-cov-2 \
--output-tree output/nextclade.auspice.json \
--output-dir output/ \
--output-basename nextclade

Run the R script Code/Quality_Control 

Step 4: Sub-sample dataset

