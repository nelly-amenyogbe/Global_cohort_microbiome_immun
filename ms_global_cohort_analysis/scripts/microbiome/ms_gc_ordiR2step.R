####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Microbiome: OrdiR2step
####################

# In this script, we will determine contribution of host factors to microbiome community composition among all children, and within each cohort separately.

# load packages
library(phyloseq)
library(ggplot2)
library(vegan)
library(mixOmics)
library(missForest)
library(plyr)
library(dplyr)
source("ms_global_cohort_analysis/scripts/functions/functions_otu_transformations.R")
source("ms_global_cohort_analysis/scripts/functions/matrix_missing_values_cleanup.R")

# load data
ps <- readRDS("ms_global_cohort_analysis/Rdata/human_raw_data/gc_physeq.rds")
meta <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_metadata.csv")
bfeed <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_breastfeeding.csv") # breastfeeding duration for Canadian and Ecuadorean children

# add metadata to physeq
ps.meta <- data.frame(sample_data(ps))

meta <- filter(meta, Microbiome == "Y")
rownames(meta) <- meta$Microbiome_ID
colnames(meta)

# create age in months
meta$age <- ifelse(is.na(meta$age.stool.months), meta$age.months, meta$age.stool.months)

# add breastfeeding data to metadata ####
meta <- join(meta, bfeed[,c("SUBJECT_ID", "ever.breastfed", "breastnow", "agesuspendbrfeed.mo")], by = "SUBJECT_ID")

# calculate time since breastfeeding
meta$time.since.bfeed <- ifelse(is.na(meta$breastnow) | meta$breastnow == "N", meta$age - meta$agesuspendbrfeed.mo, 
                                ifelse(meta$breastnow == "Y", 0, NA))

# change negative values to 0
meta$time.since.bfeed <- ifelse(meta$time.since.bfeed < 0, 0, meta$time.since.bfeed)

meta$bfeed.duration <- ifelse(is.na(meta$breastnow) | meta$breastnow == "N", meta$agesuspendbrfeed.mo,
                              ifelse(meta$breastnow == "Y", meta$age, NA))


meta.select <- meta[,c("Microbiome_ID","Cohort", "ETHNICITY", "COUNTRY", "SEX", "GEST_AGE", "BIRTH_WT", "DELIVERY", "MOM_AGE", "WAZ", "WLZ", "HAZ", "age", "time.since.bfeed", "bfeed.duration")]

# order to be same as metadata
meta.select <- meta.select[match(ps.meta$X.SampleID, meta.select$Microbiome_ID),]
length(which(meta.select$Microbiome_ID == rownames(ps.meta)))
rownames(meta.select) <- meta.select$Microbiome_ID

sampledat <- sample_data(meta.select)
ps <- merge_phyloseq(ps, sampledat)

# Adonis test for effect of Cohort ####
ord <- ordinate(ps, distance = "bray", method = "NMDS")

dist <- phyloseq::distance(ps, "bray")
ad.test <- adonis(dist ~ Cohort, data = data.frame(sample_data(ps)))
ad.test # R2 = 0.15; p = 0.001

# ordiR2step analysis for all host factors ####
# 1. Prepared data for redundancy analysis by applying the CCS transformation
# 2. Run rda on OTU data using intercept only, followed by all host factors of interest
# 3. Run ordiR2step on both models

# Normaize OTU data
otu.tab <- data.frame(t(otu_table(ps)))

# prune rare taxa
otu.f <- remove.low.counts(otu.tab, 1, 5)
otu.f <- otu.f$df 
dim(otu.f) # 2391 OTUs present in at least 5% of individuals

otu.nzv <- nearZeroVar(otu.f, freqCut = 95/5, uniqueCut = 10)

x <- otu.nzv$Metrics 
otu.nzv.f <- otu.f[,-otu.nzv$Position]
dim(otu.nzv.f) # 1856 OTUs fit the bill

# prepare metadata for capscale ####
# no missing values are allowed
meta <- data.frame(sample_data(ps))
colnames(meta)

vars <- c("Cohort","SEX", "DELIVERY", "age", "GEST_AGE", "BIRTH_WT", "MOM_AGE", "WLZ", "WAZ", "HAZ")

meta.s <- meta[,which(colnames(meta) %in% vars)]
meta.mf <- missForest(meta.s)
meta.mf$OOBerror
meta.rda <- meta.mf$ximp

# capscale with distance-based data ###
# run RDA-Cohort
cap.1 <- capscale(otu.nzv.f ~ 1, data = meta.rda, dist = "bray")
cap.all <- capscale(otu.nzv.f ~ ., data = meta.rda, dist = "bray")

# ordiR2step ####

step.cap <- ordiR2step(cap.1, scope = formula(cap.all))
step.cap$anova #Cohort only has significant finding with R2 = 0.12655 (p = 0.002).  R2 total is 0.12655.  Remainder of variables contribute an R2 of 0.011 (0.11 variance) for a total R2 of 0.13777, or 13% total variance explained.

# OrdiR2 step by cohort ####
cad <- subset_samples(ps, Cohort == "CAD")
ecd <- subset_samples(ps, Cohort == "ECD")
saf <- subset_samples(ps, Cohort == "SAF")
saf <- prune_taxa(taxa_sums(saf) > 3, saf)
blg <- subset_samples(ps, Cohort == "BLG")

physeqs <- list(cad, ecd, saf, blg)
names(physeqs) <- c("CAD", "ECD", "SAF", "BLG")

# create OTU tables for each cohort separately
otu.tabs <- llply(physeqs, function(i){
  
  otu.tab <- data.frame(t(otu_table(i)))
  otu.nzv <- remove.nzv(otu.tab)
  otu.nzv <- otu.nzv$df.no.nzv
  
  otu.f <- remove.low.counts(otu.nzv, 1, 5)
  otu.f <- otu.f$df
  
  otu.f
  
})

names(physeqs)
# create cohort-specific metadata
meta.dfs.sb <- llply(physeqs[c(3:4)], function(i){
  meta <- data.frame(sample_data(i))
  
  vars <- c("Cohort","SEX", "DELIVERY", "age", "GEST_AGE", "BIRTH_WT", "MOM_AGE", "WLZ", "WAZ", "HAZ")
  
  meta.s <- meta[,which(colnames(meta) %in% vars)]
  meta.mf <- missForest(meta.s)
  meta.rda <- meta.mf$ximp
  meta.rda
})

# Canada and Ecuador, with breastfeeding 
meta.dfs.ec <- llply(physeqs[c(1:2)], function(i){
  meta <- data.frame(sample_data(i))
  
  vars <- c("Cohort","SEX", "DELIVERY", "age", "GEST_AGE", "BIRTH_WT", "MOM_AGE", "WLZ", "WAZ", "HAZ", "bfeed.duration", "time.since.bfeed")
  
  meta.s <- meta[,which(colnames(meta) %in% vars)]
  meta.mf <- missForest(meta.s)
  meta.rda <- meta.mf$ximp
  meta.rda
})

# CAD ordiR2 ####
cad.1 <- capscale(otu.tabs$CAD ~1, data = meta.dfs.ec$CAD, dist = "bray")
cad.all <- capscale(otu.tabs$CAD ~., data = meta.dfs.ec$CAD, dist = "bray")
cad.step <- ordiR2step(cad.1, scope = formula(cad.all))
cad.step$anova # WLZ R2 = 0.030; p = 0.012

# ECD ordiR2 ####

# run ECD capscale
ecd.1 <- capscale(otu.tabs$ECD ~1, data = meta.dfs.ec$ECD, dist = "bray")
ecd.all <- capscale(otu.tabs$ECD ~ SEX + DELIVERY + MOM_AGE + WLZ + WAZ + HAZ + age + GEST_AGE + time.since.bfeed + bfeed.duration, data = meta.dfs.ec$ECD, dist = "bray")

ecd.step <- ordiR2step(ecd.1, scope = formula(ecd.all)) 
ecd.step$anova # None

# SAF ordiR2 ####
saf.1 <- capscale(otu.tabs$SAF ~1, data = meta.dfs.sb$SAF, dist = "bray")

saf.all <- capscale(otu.tabs$SAF ~ BIRTH_WT + WLZ + WAZ + HAZ + GEST_AGE, data = meta.dfs.sb$SAF, dist = "bray")

saf.all2 <- capscale(otu.tabs$SAF ~ SEX + MOM_AGE, data = meta.dfs.sb$SAF, dist = "bray")

saf.step <- ordiR2step(saf.1, scope = formula(saf.all)) 
saf.step$anova # WLZ R2 = 0.174, p = 0.016

saf.step2 <- ordiR2step(saf.1, scope = formula(saf.all2)) 
saf.step2$anova # none

# BLG ordiR2 ####

blg.1 <- capscale(otu.tabs$BLG ~1, data = meta.dfs.sb$BLG, dist = "bray")

blg.all <- capscale(otu.tabs$BLG ~ BIRTH_WT + WLZ + WAZ + HAZ, data = meta.dfs.sb$BLG, dist = "bray")

blg.all2 <- capscale(otu.tabs$BLG ~ GEST_AGE + BIRTH_WT + MOM_AGE, data = meta.dfs.sb$BLG, dist = "bray")

blg.step <- ordiR2step(blg.1, scope = formula(blg.all)) # no significant effects
blg.step$anova # none

blg.step2 <- ordiR2step(blg.1, scope = formula(blg.all2)) # no significant effects
blg.step2$anova # none

# END ####