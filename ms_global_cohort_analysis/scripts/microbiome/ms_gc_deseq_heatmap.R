####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Microbiome: Deseq2 for differential abundance: LRT test
####################

# In this script, we will determine differentially abundant OTUs across cohorts using DESeq2.  The most discriminatory taxa will be selected using sPLS-DA, and visualized in a heatmap.

# load packages
library(phyloseq)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(DESeq2)
library(mixOmics)
library(pheatmap)
# helper functions
source("ms_global_cohort_analysis/scripts/functions/functions_otu_transformations.R")
source("ms_global_cohort_analysis/scripts/functions/matrix_missing_values_cleanup.R")

# load data
ps <- readRDS("ms_global_cohort_analysis/Rdata/human_raw_data/gc_physeq.rds")

# load colour attributes for plotting
gc.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_mb_cols.csv")
phylum.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_phylum_cols.csv")
family.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_network_fam_cols.csv")


# Filter OTU table to only include abundant OTUs ####
# Create matrix of OTU data
otu.tab <- t(data.frame(otu_table(ps)))

# retain only OTUs present in at least 6% of samples
otu.f <- remove.low.counts(as.data.frame(otu.tab), 1, 6)
otu.f <- otu.f$df 
dim(otu.f) # 2075 OTUs

# subset physeq to these taxa
taxtab <- data.frame(tax_table(ps))

ps.f <- subset_taxa(ps, OTU.name %in% colnames(otu.f))

# run DESeq2 ####
design <- ~Cohort

dds <- phyloseq_to_deseq2(ps.f, design)

ds.lrt <- DESeq(dds, test = "LRT", reduced = ~1)

# Get deseq results ####
res <- data.frame(results(ds.lrt))
res$OTU.name <- rownames(res)
# Add taxonomy to results
res <- join(res, taxtab, by = "OTU.name")

# get significant results only
sig <- filter(res, padj < 0.01)
dim(sig) # 462 OTUs at padj = 0.01
sig <- sig[order(sig$padj),]


# PLS-DA for significant OTUs ####
# prepare normalized OTU data
# get transformed counts
size.fac <- estimateSizeFactors(dds)
otu.dds <- counts(size.fac, normalized = TRUE)

# get variance stabalized counts
otu.vst <- varianceStabilizingTransformation(dds, fitType = "local")
otu.mat <- assay(otu.vst)
otu.mat <- t(otu.mat)

otu.sig <- otu.mat[,which(colnames(otu.mat) %in% sig$OTU.name)]

# run PLS-DA
X <- otu.sig
Y <- sample_data(ps)$Cohort

# GC PLSDA ####
gc.plsda <- splsda(X = X,
                   Y = Y,
                   ncomp = 4)

# Heatmap of top discriminatory taxa ####

# get top taxa
# function to retrieve top features from a PLS object
mc.get.features <- function(mixomics.obj, comp){
  
  # get loadings from mixomics object
  loadings <- data.frame(mixomics.obj$loadings$X)
  loadings$feature <- rownames(loadings)
  
  # modify data frame to only include comp of interest
  loadings <- loadings[,c("feature", comp)]
  loadings$raw.loadings <- loadings[,2]
  
  # add absolute loadings and order features by strength
  loadings$abs.loadings <- abs(loadings$raw.loadings)
  loadings <- loadings[order(loadings$abs.loadings, decreasing = TRUE),]
  
  # remove redundant column
  loadings <- loadings[,-2]
  return(loadings)
  
}

# prepare data
length(unique(sig$Family)) # 37 families in total

c1.load <- mc.get.features(gc.plsda, "comp.1")
features <- c1.load$feature[1:50]

# format data matrix
mat <- X[,which(colnames(X) %in% features)]

# Create plotting taxonomy table ####
# create labels of taxonomic rank with otu name
tax.sig <- filter(taxtab, OTU.name %in% sig$OTU.name)
tax.sig$lab <- ifelse(as.character(tax.sig$Genus) == as.character(tax.sig$Species),
                      paste(as.character(tax.sig$OTU.name), as.character(tax.sig$Genus)),
                      paste(as.character(tax.sig$OTU.name), as.character(tax.sig$Genus), as.character(tax.sig$Species)))

# truncated version, with species names where applicable
tax.sig$lab.short <- ifelse(as.character(tax.sig$Genus) == as.character(tax.sig$Species),
                            paste(as.character(tax.sig$Genus)),
                            paste(as.character(tax.sig$Genus), as.character(tax.sig$Species)))

# prepare row annotations (OTU class)
row.tax <- filter(tax.sig, OTU.name %in% features)
row.annot <- row.tax[order(colnames(mat)),]

rownames(row.annot) <- row.annot$lab
row.annot <- subset(row.annot, select = c("Phylum", "Family"))

# prepare column annotations (Cohort)
col.annot <- data.frame(sample_data(ps))
col.annot <- col.annot[order(rownames(mat)),]
rownames(col.annot) <- col.annot$X.SampleID
col.annot <- subset(col.annot, select = c("Cohort"))

colnames(mat) <- row.tax$lab
mat <- t(mat)

# set colour attributes
cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')


# specific cols: phylum
phyla <- unique(as.character(row.annot$Phylum))
phy.select <- filter(phylum.cols, Phylum %in% phyla)
phy.cols <- as.character(phy.select$phylum.cols)
names(phy.cols) <- phy.select$Phylum

# family
fams <- unique(as.character(row.annot$Family))
fam.select <- filter(family.cols, Family %in% fams)

# diffs
length(which(row.annot$Family %in% family.cols$Family))
diffs <- setdiff(row.annot$Family, family.cols$Family)
diffs # "F_Rickettsiaceae" 

cols.extra <- family.cols[-which(family.cols$col %in% fam.select$col),]
cols.extra <- unique(cols.extra$col)

fam.cols <- as.character(fam.select$col)
names(fam.cols) <- fam.select$Family
fam.cols <- c(fam.cols, "F_Rickettsiaceae" = as.character(cols.extra[1]))

# cohort
coh.cols <- as.character(gc.cols$Color)
names(coh.cols) <- gc.cols$Site

ann_colors = list(Phylum = phy.cols,
                  Family = fam.cols,
                  Cohort = coh.cols)

pheatmap(mat,
         show_colnames = FALSE,
         annotation_row = row.annot,
         annotation_col = col.annot,
         annotation_colors = ann_colors,
         labels_row = as.character(row.tax$lab.short),
         fontsize_row = 8)

# save as ms_gc_heatmap

# END ####