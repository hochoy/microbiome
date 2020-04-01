# Analysis of microbiome 16S rRNA data (fastq format)

Following the experimental design of 

> Microbiota changes 		= ~ diet + phase + diet:phase

There are a number of approaches that take into account interactions between group variables, whilst managing multiple testing corrections and other errors in the underlying data creation and processing.
In this example, we use 3 `R` language packages:

1. DADA2
2. phyloseq
3. DESeq2
4. `+` good old tidyverse


## Overview of steps

1. Removal of primers ←- Not performed due to lack of primer knowledge
2. Filter and truncate reads based on base-calling quality
3. Identify and extract unique sequences from reads
4. Merge read-pairs based on overlapping sequences
5. Remove chimeras based on parent sequences
6. Assign taxonomy (Silva dataset)
7. Plot taxonomic distribution (phyloseq)
8. Plot Richness and Diversity (phyloseq)
9. Plot Ordination (phyloseq)
10. Detect metrics (taxonomy, genes, sequences) that have changed in distribution across sample groups via DESeq2

Steps 1 to 6 are common sequence pre-processing steps with a number of popular tools for each purpose, i.e. fastqc, QIIME, mothur, trimmomatic. Here, we use DADA2

Steps 7 to 8 are common plots showing the distribution of taxonomic groups, and microbial richness and diversity across samples.

Step 9 or Ordination allows us to reduce a high-dimensional dataset into a low dimensional dataset (2D graph) for easier human visualization. 

Step 10 I am not an expert on the DESeq2 method, but in simple terms, DESeq2 attempts to fit our data onto a negative binomial curve based on the experimental design (variables + interaction terms) and determine if any of the model's coefficients are statistically different from zero; i.e. have a statistically significant effect on the overall model.

