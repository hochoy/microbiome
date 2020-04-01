# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.10")
# BiocManager::install("phyloseq")
# BiocManager::install("Biostrings")
# BiocManager::install("DESeq2")
# install.packages("tidyverse")

# packageVersion("tidyverse") # 1.3.0
# packageVersion("dada2") # 1.14.1
# packageVersion("phyloseq") # 1.30.0
# packageVersion("Biostrings") # 2.54.0
# packageVersion("DESeq2") # 1.26.0


library(dada2); 
library(tidyverse);
library(R.utils)
library(phyloseq)
library(Biostrings)
library(DESeq2)

# DADA2: https://benjjneb.github.io/dada2/tutorial.html

# helper function to get files in this dir
get_this_path <- function(){
  path <- rstudioapi::getSourceEditorContext()$path
  dir <- dirname(path)
  list(path = path, dir=dir)
}

# get files in this dir
(curr_dir <- get_this_path()$dir)
(dir_files <- list.files(curr_dir))

# target folders for read/write
(fastq_gz_folder <- paste0(curr_dir,"/fastq_gz/"))
(fastq_folder <- paste0(curr_dir,"/fastq/"))


# gunzip files
(fastq_gz_paths <- list.files(fastq_gz_folder,full.names = TRUE))

# fastq_paths <- sapply(fastq_gz_paths, function(gz_path){
#   out_name <- basename(gz_path) %>% str_replace("\\.gz","")
#   out_path <- paste0(fastq_folder,out_name)
#   
#   gunzip(gz_path,out_path,TRUE,FALSE)
#   out_path
# })
(fastq_paths <- list.files(fastq_folder,full.names=TRUE))


# Place fastq files into a df
(fast_df <- data.frame(
  raw_fq = fastq_paths
) %>% 
  mutate(
    phase = str_extract(raw_fq,"Pre|Post"),
    diet = str_match(raw_fq,"Trans\\.(EEN|HF|LF)")[,2],
    mouse_id = paste(diet, str_match(raw_fq,"\\.(M[0-9]+)_R")[,2],sep="_"),
    sample_id = paste(phase,mouse_id,sep="_"),
    grp_id = paste(phase,diet,sep="_"),
    read_pair = str_match(raw_fq, "_(R[1,2])\\.fastq")[,2]
  )
) %>% head()


# view a sample fastq file
(forward_qualities <- plotQualityProfile(
  fast_df %>% filter(read_pair == 'R1') %>% .$raw_fq,
  aggregate = TRUE))

(reverse_qualities <- plotQualityProfile(
  fast_df %>% filter(read_pair == 'R2') %>% .$raw_fq,
  aggregate = TRUE))

forward_qualities + expand_limits(x= 260)

# split the rows by forward/reverse-read, making sure to sort the rows to match the sample
(forward <- fast_df %>% arrange(sample_id) %>% filter(read_pair == "R1"))
(reverse <- fast_df %>% arrange(sample_id) %>% filter(read_pair == "R2"))
(forward_paths <- forward$raw_fq %>% sapply(toString))
(reverse_paths <- reverse$raw_fq %>% sapply(toString))
sample.names <- forward$sample_id


# assign filename for filtered files
(filtFs <- forward_paths %>% str_replace("\\.fastq","_Filt\\.fastq.gz"))
(filtRs <- reverse_paths %>% str_replace("\\.fastq","_Filt\\.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# filter reads and trim 
# trunLen: truncate the ends of the reads by a fixed length
# maxN:
# maxEE: maximum error rate allowed in a read before discarded
# truncQ: sequences with values below this threshold are truncated
# rm.phix: PhiX is a control spike of DNA for calibration and QC of sequencing runs
out <- filterAndTrim(forward_paths,filtFs,reverse_paths,filtRs,
                     truncLen = c(240, 160),
                     maxN = 0,
                     maxEE = c(2,2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress=TRUE, multithread = TRUE)

# Fit potential errors in base-calling to a loess curve
# and see if it matches expected transition rates with increasing quality
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


# sample inference: identifying unique sequence from reads
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)

dadaFs[[1]]
dadaFs[[2]]


# merging forward/reverse reads using the unique sequences 
# and some alignment rules
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

# construct an Amplicon Sequence Variant (ASV) table
seqtab <- makeSequenceTable(mergers)

# sanity-check sequence lengths
table(nchar(getSequences(seqtab)))

nchar(getSequences(seqtab)) %>% hist(breaks=20)
#hmm, looks like there might be primers, if expected length is ~250

# remove chimeras using parent sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# final check of processing pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track_w_perc <- track %>% as.data.frame() %>%  mutate(perc_final = round(nonchim/input * 100,1))
row.names(track_w_perc) <- row.names(track)
track_w_perc
# assign taxonomy using Silva
taxa <- assignTaxonomy(seqtab.nochim, paste0(curr_dir,"/silva/silva_nr_v138_train_set.fa.gz"), multithread=TRUE)

# optionally, add species assignment by exact matching to ASV sequences
taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v138.fa.gz")

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# some phyloseq analysis
theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)

sample_table <- forward
rownames(sample_table) <- samples.out

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(sample_table),
  tax_table(taxa)
  )

# Rename ASV to shorter strings
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# plot richness
plot_richness(ps, x="phase", measures=c("Observed","Shannon","Simpson"), color="diet")
estimate_richness(ps)
plot_richness(measures)
# plot ordination NMDS
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="grp_id", title="Bray NMDS")


# plot top 20 taxa
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="grp_id", fill="Family") + facet_wrap(~grp_id, scales="free_x",ncol = 6)

# Deseq2
# http://joey711.github.io/phyloseq-extensions/DESeq2.html
# deseq2_data <- phyloseq_to_deseq2(ps, ~ diet + phase + diet:phase)
deseq2_data <- phyloseq_to_deseq2(ps, ~ diet + phase + diet:phase)
deseq2_data <- DESeq(deseq2_data, test="Wald", fitType = "parametric")

# view results
resultsNames(deseq2_data)
# res = results(deseq2_data,
#               contrast = c("diet","EEN","HF"),
#               cooksCutoff = FALSE)
res <- results(deseq2_data,contrast = c("diet","HF","LF"))
# res <- results(deseq2_data,contrast = c("diet","EEN","LF"))
# res <- results(deseq2_data,contrast = c("diet","HF","LF"))
# res <- results(deseq2_data,contrast = c("diet","EEN","HF"))
# res <- results(deseq2_data,contrast = c("phase","Post","Pre"))
res
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

# plot 
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
