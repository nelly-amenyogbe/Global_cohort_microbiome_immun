####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# sPLS: OTU data preparation
####################

# In this script, we prepare microbiome data for use with sPLS.  This includes filtering the OTU tables to only select OTUs repesented across at least 30% of subjects within each cohort, and normalizing using the CLR transformation.

# load packages
library(phyloseq)
library(mixOmics)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
# helper functions
source("ms_global_cohort_analysis/scripts/functions/matrix_missing_values_cleanup.R")
source("ms_global_cohort_analysis/scripts/functions/functions_otu_transformations.R")

# load data
# OTU data (Phyloseq object)
ps <- readRDS("ms_global_cohort_analysis/Rdata/human_raw_data/gc_physeq.rds")

# luminex data
lmx.dfs <- readRDS("ms_global_cohort_analysis/Rdata/R_export/ms_gc_lmx_filtered_sPLS.rds")
lmx.all <- lmx.dfs$all

#### Prepare OTU data for filtering ####
gc <- subset_samples(ps, InnateID %in% unique(as.character(lmx.all$Sample)))

# create separate phyloseq objects for every cohort
cohorts <- c("CAD", "SAF", "ECD", "BLG")

cad <- subset_samples(gc, Cohort == "CAD")
saf <- subset_samples(gc, Cohort == "SAF")
ecd <- subset_samples(gc, Cohort == "ECD")
blg <- subset_samples(gc, Cohort == "BLG")

ps <- list(cad, saf, ecd, blg)
names(ps) <- cohorts

# prune each cohort to only include OTUs present in that cohort 

otu.tabs <- llply(ps, function(i){
  
  tab <- data.frame(t(otu_table(i)))
  rownames(tab) <- sample_data(i)$InnateID
  tab
  
})

otu.tss.tabs <- llply(otu.tabs, function(i){
  
  otu.tss <- TSS.transform(i)
  otu.tss <- as.data.frame(otu.tss)
  otu.tss
  
})

x <- otu.tss.tabs$CAD

otu.rlc <- llply(otu.tss.tabs, function(i){
  
  tab <- remove.low.counts(i, 1e-5, 50)
  return(tab$df)
  
})

otus <- llply(otu.rlc, function(i){
  
  otus <- as.character(colnames(i))
  otus
  
})

llply(otu.rlc, function(i){
  
  dim(i)
  
}) # check how many OTUs were retained

overlap.otus <- Reduce(intersect, otus)
length(overlap.otus) # 135 OTUs selected across cohorts

# Prepare OTU table from all samples together
gc.otu <- as.data.frame(otu_table(gc))
gc.meta <- data.frame(sample_data(gc))
colnames(gc.otu) <- gc.meta$InnateID

gc.otu.s <- gc.otu[which(rownames(gc.otu) %in% overlap.otus),]
gc.otu.s <- as.data.frame(t(gc.otu.s))

otu.rlc <- c(otu.rlc, list(gc.otu.s))
names(otu.rlc)[5] <- "ALL"

# Normalized OTU datasets ####

# CLR
otu.clr <- llply(otu.rlc, function(i){
  mat <- as.matrix(i)
  mat.clr <- clr(mat)
  mat.clr
})

# save relevant files ####
#saveRDS(otu.clr, "ms_global_cohort_analysis/Rdata/R_export/otu_clr.rds")

# END ####
