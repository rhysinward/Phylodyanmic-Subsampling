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

Step 3: Quality control

* See here for Next clade pipeline - https://docs.nextstrain.org/projects/nextclade/en/stable/
* Run dataset through IQTREE to generate ML-tree
* Process using tempest and remove erroneous sequences using Code/remove_sequences_from_tempest.R

Step 4: Sub-sample dataset

* Run script Code/sampling_code_Hong_Kong.R to generate sub-sampled Hong Kong dataset and plots

Step 5: Analysis  

* Run dataset through BEAST/BEAST2 (see XML's) for TMRCA + BDSKY
* Run dataset through Code/Skygrowth_ML+BEAST.R for Skygrowth

Step 6: Analysis and plotting of results

* Run BDSKY model through Code/Skyline_BDSKY.r
* Subsequently run through Code/final_merge_HPD_BDSKY.R or final_merge_HPD_Skygrowth.R to plot EpiFilter against BDSKY or Skygrowth
* Run Code/Jensen_Shannon_Divergence_BDSKY_Hong_Kong.R to get JSD
* Run Code/plot_TMRCA_and_R0.R to get TMRCA and R0 plots 
* Run Code/substitution_rate_plot.R to get substitution rate plot
