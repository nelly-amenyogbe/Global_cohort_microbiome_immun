####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# sPLS: hyperheometric test for selected OTUs
####################

# In this script, we perform the hypergeometric test to determine whether any cytokines were enriched for in the sPLS results for each cohort

# load packages
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)

# load data 
all.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/blg_cad_ecd_spls_res.rds")
blg.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/blg_spls_res.rds")
cad.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/cad_spls_res.rds")
ecd.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/ecd_spls_res.rds")
saf.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/saf_spls_res.rds")

stim.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_stim_cols.csv")

# taxonomy table
taxonomy <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_taxonomy_table.csv")
taxonomy$OTU.num <- taxonomy$X
taxonomy$names <- paste(taxonomy$OTU.num, taxonomy$Genus) 

# run hg test: all ####

# function for hyper test ####

otu.phyper.test <- function(spls.res, features, cohort, comp){
  
  # list of all input features
  otu.vars <- data.frame(colnames(spls.res$spls.obj$X))
  colnames(otu.vars)[1] <- "names"
  # add families
  otu.meta <- join(otu.vars, taxonomy[,c("names", "Family")], by = "names")
  
  
  # list of selected features
  select <- data.frame(features)
  colnames(select)[1] <- "names"
  select.meta <- join(select, taxonomy[,c("names", "Family")], by = "names")
  
  
  # filter to only include families that are in the dataset at least 4 times
  freq <- data.frame(table(as.character(otu.meta$Family)))
  no.fams <- length(otu.meta$Family)
  freq$occurence <- (freq$Freq/no.fams)*100
  top <- filter(freq, Freq > 4)
  families <- unique(as.character(top$Var1))
  
  
  all.res <- ldply(as.list(families), function(i){
    
    qt <- length(which(select.meta$Family == i)) # success
    mt <- length(which(otu.meta$Family == i)) # number of feats available
    nt <- length(otu.meta$Family) - mt # number of non-stim to choose from
    kt <- length(select.meta$Family) # number of features selected
    
    # test for over-representation
    res <- phyper(q = qt, m = mt, n = nt, k = kt, lower.tail = FALSE)
    
    Family <- i
    to.return <- c(res, Family)
    
  })
  
  colnames(all.res) <- c("p.value", "Family")
  all.res$cohort <- cohort
  all.res$component <- comp
  all.res$p.adj <- p.adjust(all.res$p.value, method = "BH")
  all.res
  
}

# do all tests
# all samples ####
all.res <- otu.phyper.test(all.spls, all.spls$otus.select, "all", "comp.1")

# BLG ####
# comp 1 and 2

blg1.res <- otu.phyper.test(blg.spls, blg.spls$otu1.select, "BLG", "comp.1")

blg2.res <- otu.phyper.test(blg.spls, blg.spls$otu2.select, "BLG", "comp.2")

# CAD ####
# comp1 and 2
cad1.res <- otu.phyper.test(cad.spls, cad.spls$otu1.select, "CAD", "comp.1")

cad2.res <- otu.phyper.test(cad.spls, cad.spls$otu2.select, "CAD", "comp.2")

# ECD ####
# comp 1 and 3
ecd1.res <- otu.phyper.test(ecd.spls, ecd.spls$otu1.select, "ECD", "comp.1")

ecd3.res <- otu.phyper.test(ecd.spls, ecd.spls$otu3.select, "ECD", "comp.3")

# SAF ####
# comps 1, 2, and 3
saf1.res <- otu.phyper.test(saf.spls, saf.spls$otu1.select, "SAF", "comp.1")

saf2.res <- otu.phyper.test(saf.spls, saf.spls$otu2.select, "SAF", "comp.2")

saf3.res <- otu.phyper.test(saf.spls, saf.spls$otu3.select, "SAF", "comp.3")

# combine and plot #### 
combined.res <- rbind(all.res, blg1.res, blg2.res, cad1.res, cad2.res, ecd1.res, ecd3.res, saf1.res, saf2.res, saf3.res)

# sig findings
sig <- filter(combined.res, p.adj < 0.05)
sig # prevotella for all/ECD.
# END ####
