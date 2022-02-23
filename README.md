# Phylodyanmic-Subsampling


# Repository description:
This code and data were used for the analysis presented in “Using multiple sampling strategies to estimate SARS-CoV-2 epidemiological parameters from genomic sequencing data" - https://www.medrxiv.org/content/10.1101/2022.02.04.22270165v1.

The repository contains the following elements:

1. Data

Contains some of the data (genomic data has been excluded due to permissions) used for analysis that is:

* Case data for Amazonas State, Brazil and Hong Kong
* XML files
* Accession ID's for genomic data

2. Data Processing Script

Contains the scripts needed for processing sequence metadata and sequence data using R. 

3. Plotting Script

Contains modified script for Treemmer which allows for pairs to be protected using Python.

# Running the code

In this example we show how to run the process for Hong Kong:

Run the R script Code/sampling.R 

* The output of this is a migration-informed sequenced data

Align sequences to reference sequence using minimap2 embedded within Pangolin

conda activate pangolin
pangolin --alignment gisaid_hcov-19_2022_01_20_09ta_processedNames.fas

Quality control through the NextClade pipeline

nextclade \
--in-order \
--input-fasta data/sars-cov-2/sequences.aln.fasta \
--output-tsv output/nextclade.tsv \
--input-dataset data/sars-cov-2 \
--output-tree output/nextclade.auspice.json \
--output-dir output/ \
--output-basename nextclade

Run the R script Code/Quality_Control 

* Output of this is QC migration-informed sequence data.

Create ML-tree using IQtree

iqtree -s -m GTR+I+G -alrt 1000

Run the R script Code/treemmer_metadata.R

* Output of this is metadata needed to run treemmer

Run the Python Script Code/Treemmer_v0.3_mod.py

python3 Code/Treemmer_v0.3_mod.py  -lm Data/treemmer_metadata_519.txt -lmc Data/metadata_proportions.txt  Data/Mex_B_1_1_519_12_07_2021.fasta.treefile -pp -X 4000

Run R script Code/prune_treemmer.R

* Output of this is 519 sequences after pruning by treemmer

Re-run IQtree to get final ML-tree
