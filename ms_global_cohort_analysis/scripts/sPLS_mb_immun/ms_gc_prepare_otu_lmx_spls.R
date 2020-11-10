####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# sPLS: OTU_LMX matrix preparation
####################

# In this script, we prepare microbiome and OTU data for sPLS.  This involves creating data matrices from luminex data, and ensuring both data types have matching row names. These matrices are saved for ease of use for sPLS models.

# load packages
library(plyr)
library(dplyr)
library(reshape2)
library(missForest)
# helper functions
source("ms_global_cohort_analysis/scripts/functions/matrix_missing_values_cleanup.R")

# load data 
# microbiome
clr.dat <- readRDS("ms_global_cohort_analysis/Rdata/R_export/otu_clr.rds")

# luminex
lmx <- readRDS("ms_global_cohort_analysis/Rdata/R_export/ms_gc_lmx_filtered_sPLS.rds")

# taxonomy
taxonomy <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_taxonomy_table.csv")
taxonomy$OTU.num <- taxonomy$X

# Prepare luminex matrices ####
# make matrices
lmx.mats <- llply(lmx, function(i){
  mat <- dcast(i, Sample ~ cyto.stim, value.var = "final.concentration")
  rownames(mat) <- mat$Sample
  mat <- mat[,-1]
})

# remove cytokines with >15% missing values
lmx.beads <- llply(lmx.mats, function(i){
  
  mat <- beads.rm(i, 15)
  
})

lmx.beads # none to remove

# remove subjects with >15% missing data
lmx.sub <- llply(lmx.mats, function(i){
  mat <- which.sub.rm(i, 15)
})

lmx.sub #ECD1010, SAF1031, BLG1033, SAF1026 removed

lmx.mats <- llply(lmx.mats, function(i){
  mat <- subj.rm(i, 15)
})

# prepare OTU tables to include taxonomic names
otu.clr.tax <- llply(clr.dat, function(i){
  tax.f <- filter(taxonomy, OTU.num %in% colnames(i))
  tax.f <- tax.f[order(tax.f$OTU.num),]
  tax.f$OTU.names <- paste(tax.f$OTU.num, tax.f$Genus)
  colnames(i) <- tax.f$OTU.names
  return(i)
})

# Prepare sPLS matrices ####
names(lmx.mats)
names(otu.clr.tax)
lmx.mats.reordered <- lmx.mats[c(3, 4, 2, 1, 5)]
names(lmx.mats.reordered) # ensure all are in same order as microbiome data

spls.data <- llply(as.list(1:length(lmx.mats)), function(i){
  
  X <- otu.clr.tax[[i]]
  Y <- lmx.mats.reordered[[i]]
  
  # get overlapping subjects
  subs <- Reduce(intersect, list(rownames(X), rownames(Y)))
  
  # Arrange matrices to have the same subjects and in the same order
  X <- X[rownames(X) %in% subs,]
  X <- X[order(rownames(X)),]
  
  Y <- Y[rownames(Y) %in% subs,]
  Y <- Y[order(rownames(Y)),]
  Y <- missForest(Y)
  Y <- Y$ximp
  
  to.return <- list(X, Y)
  names(to.return) <- c("Xotu", "Ylmx")
  return(to.return)
  
})
names(spls.data) <- names(otu.clr.tax)

# check matching rownames
llply(as.list(1:length(spls.data)), function(i){
  dat <- spls.data[[i]]
  X <- dat$Xotu
  Y <- dat$Ylmx
  
  test <- length(which(rownames(X) == rownames(Y))) / nrow(X)
  test
}) # if all equal 1 then you are good to go

# save data ###
#saveRDS(spls.data, "ms_global_cohort_analysis/Rdata/R_export/spls_data.rds")

# END ####