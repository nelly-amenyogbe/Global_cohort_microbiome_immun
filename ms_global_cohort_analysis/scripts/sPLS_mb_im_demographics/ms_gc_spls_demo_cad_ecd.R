#####################
# Nelly Amenyogbe
# Global Cohort Manuscript Analysis
# sPLS: microbiome-immune demographics: CANADA and ECUADOR
#####################

# In this script, we determine whether host demographic factors associated with microbiome-immune correlations among both Canadian and Ecuadorean children.


# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(reshape2)
library(mixOmics)
library(cowplot)
# helper function
source("ms_global_cohort_analysis/scripts/functions/function_pairwise_lm.R")
source("ms_global_cohort_analysis/scripts/functions/function_block_spls_loadings.R")
source("ms_global_cohort_analysis/scripts/functions/function_pairwise_lm.R")

# load data
meta.all <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_metadata.csv")
bfeed <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_breastfeeding.csv")
spls.dat <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_data.rds")

# colours that will be used for plotting
cohort.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_mb_cols.csv")

# prepare data for sPLS ####
# spls data
spls.data <- spls.dat$ALL

# add breastfeeding data to metadata ####
# create age in months
meta.all$age <- ifelse(is.na(meta.all$age.stool.months), meta.all$age.months, meta.all$age.stool.months)

meta.all <- join(meta.all, bfeed[,c("SUBJECT_ID", "ever.breastfed", "breastnow", "agesuspendbrfeed.mo")], by = "SUBJECT_ID")

# calculate time since breastfeeding
meta.all$time.since.bfeed <- ifelse(is.na(meta.all$breastnow) | meta.all$breastnow == "N", meta.all$age - meta.all$agesuspendbrfeed.mo, 
                                ifelse(meta.all$breastnow == "Y", 0, NA))

# change negative values to 0
meta.all$time.since.bfeed <- ifelse(meta.all$time.since.bfeed < 0, 0, meta.all$time.since.bfeed)

meta.all$bfeed.duration <- ifelse(is.na(meta.all$breastnow) | meta.all$breastnow == "N", meta.all$agesuspendbrfeed.mo,
                              ifelse(meta.all$breastnow == "Y", meta.all$age, NA))


# prepare spls metadata ####
# metadata
spls.subjects <- c(rownames(spls.dat$CAD$Xotu), rownames(spls.dat$ECD$Xotu))

meta <- filter(meta.all, SUBJECT_ID %in% spls.subjects)

rownames(meta) <- meta$SUBJECT_ID
meta <- meta[order(meta$SUBJECT_ID),]

# set factor levels to metadata
Z <- meta[,c("SEX", "GEST_AGE", "DELIVERY", "MOM_AGE", "WAZ", "WLZ", "HAZ", "time.since.bfeed", "bfeed.duration")]

rownames(Z) <- meta$SUBJECT_ID

Z$SEX <- ifelse(Z$SEX == "Male", 0,
                ifelse(Z$SEX == "Female", 1, Z$SEX))

Z$DELIVERY <- ifelse(Z$DELIVERY == "Caesarean", 0,
                     ifelse(Z$DELIVERY == "Vaginal", 1, Z$DELIVERY))

Z <- as.matrix(Z)

# CAD_ECD sPLS ####
X <- spls.dat$ALL$Xotu
Xs <- X[which(rownames(X) %in% spls.subjects),]
Xs <- Xs[order(rownames(Xs)),]

Y <- log10(spls.dat$ALL$Ylmx)
Ys <- Y[which(rownames(Y) %in% spls.subjects),]
Ys <- Ys[order(rownames(Ys)),]

length(which(rownames(Xs) == rownames(Ys)))
length(which(rownames(Xs) == rownames(Z)))

# set parameters
dat <- list(mb = Xs, lmx = Ys, meta = Z)

design = matrix(c(1,1,1,1,1,1,1,1,1), byrow = TRUE, nrow = 3)
diag(design) <- 0

res <- block.spls(X = dat,
                  indY = 3,
                  ncomp = 3,
                  keepX = list(mb = c(20, 20, 20), lmx = c(20, 20, 20), meta = c(7, 7, 7)),
                  design = design,
                  mode = "regression")

# plot sPLS results ####

# Correlation circle
cor.circle <- plotVar(res,
                      col = c("blue", "red", "#006d2c"),
                      cex = c(2,2,2),
                      comp = c(1,2),
                      overlap = TRUE,
                      title = "ECD_CAD CLR",
                      plot = TRUE)


plotLoadings(res, comp = 1, title = "ECD_CAD comp1")

plotLoadings(res, comp = 2, title = "ECD_CAD comp2")

##### Regression: Z:Z ####
# OTU:meta
# get lm data
lm.dat.xz <- lm.df.dataprep(Xs, Z)
lm.dat.yz <- lm.df.dataprep(Ys, Z)
lm.dat <- rbind(lm.dat.xz, lm.dat.yz)

lm.dat$Cohort <- substr(lm.dat$sample.id, start = 1, stop = 1)
lm.dat$Cohort <- ifelse(lm.dat$Cohort == "E", "ECD",
                        ifelse(lm.dat$Cohort == "F", "CAD", lm.dat$Cohort))

# run pairwise lm
lm.res.xz <- x.y.lm.scaled(lm.dat.xz)
lm.res.yz <- x.y.lm.scaled(lm.dat.yz)

lm.res <- rbind(lm.res.xz, lm.res.yz)

# cohort-specific
cad.dat <- filter(lm.dat, Cohort == "CAD")
lm.res.cad <- x.y.lm.scaled(cad.dat)

ecd.dat <- filter(lm.dat, Cohort == "ECD")
lm.res.ecd <- x.y.lm.scaled(ecd.dat)

# get loadings
loadings.x <- data.frame(res$loadings$mb)
loadings.y <- data.frame(res$loadings$lmx)
loadings.z <- data.frame(res$loadings$meta)

# select vars to plot
x.ld <- block.spls.get.features(loadings.x, "comp.1")
y.ld <- block.spls.get.features(loadings.y, "comp.1")

# get data to plot
# get features
x.feats <- x.ld$feature[1:15]
y.feats <- y.ld$feature[1:8]
feats <- c(x.feats, y.feats)

# CAD_ECD
lm.res.select <- filter(lm.res, X.var %in% feats, Y.var %in% c("MOM_AGE", "HAZ", "WAZ", "MOM_AGE", "WLZ", "SEX", "GEST_AGE", "time.since.bfeed", "bfeed.duration"))
lm.res.select$Estimate <- ifelse(lm.res.select$p.value < 0.05, lm.res.select$Estimate, 0)

# CAD
cad.lm.res.select <- filter(lm.res.cad, X.var %in% feats, Y.var %in% c("MOM_AGE", "HAZ", "WAZ", "MOM_AGE", "WLZ", "SEX", "GEST_AGE", "time.since.bfeed", "bfeed.duration"))
cad.lm.res.select$Estimate <- ifelse(cad.lm.res.select$p.value < 0.05, cad.lm.res.select$Estimate, 0)

# ECD
ecd.lm.res.select <- filter(lm.res.ecd, X.var %in% feats, Y.var %in% c("MOM_AGE", "HAZ", "WAZ", "MOM_AGE", "WLZ", "SEX", "GEST_AGE", "time.since.bfeed", "bfeed.duration"))
ecd.lm.res.select$Estimate <- ifelse(ecd.lm.res.select$p.value < 0.05, ecd.lm.res.select$Estimate, 0)

# prepare hm data
x.hmdat <- dcast(lm.res.select, X.var ~ Y.var, value.var = "Estimate")
rownames(x.hmdat) <- x.hmdat$X.var
x.hmdat <- x.hmdat[,-1]

# set internal colors
color.edge <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))
breaks <- seq(from = -1, to = 1.01, by = 0.01)
cell.cols <- color.edge(length(breaks)-1)

# plot heatmap
pheatmap(x.hmdat,
         color = cell.cols,
         breaks = breaks,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "CAD_ECD: COMP1")

# exported as figures/sPLS_mb_im_demographics/cad_ecd/cad_ecd_hm_comp1.pdf

# CAD HM ####
cad.x.hmdat <- dcast(cad.lm.res.select, X.var ~ Y.var, value.var = "Estimate")
rownames(cad.x.hmdat) <- cad.x.hmdat$X.var
cad.x.hmdat <- cad.x.hmdat[,-1]

# set internal colors
color.edge <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))
breaks <- seq(from = -1, to = 1.01, by = 0.01)
cell.cols <- color.edge(length(breaks)-1)

# plot heatmap
pheatmap(cad.x.hmdat,
         color = cell.cols,
         breaks = breaks,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "CAD")

# exported as figures/sPLS_mb_im_demographics/cad_ecd/can_only_comp1.pdf

# ECD HM ####
ecd.x.hmdat <- dcast(ecd.lm.res.select, X.var ~ Y.var, value.var = "Estimate")
rownames(ecd.x.hmdat) <- ecd.x.hmdat$X.var
ecd.x.hmdat <- ecd.x.hmdat[,-1]

# set internal colors
color.edge <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))
breaks <- seq(from = -1, to = 1.01, by = 0.01)
cell.cols <- color.edge(length(breaks)-1)

# plot heatmap
pheatmap(ecd.x.hmdat,
         color = cell.cols,
         breaks = breaks,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "ECD")

# exported as figures/sPLS_mb_im_demographics/cad_ecd/ecd_only_comp1.pdf

# Plot selected features ####

# explore top comp1 loadings ####
plotLoadings(res, comp = 1, title = "ECD_CAD comp1") # HAZ and MIP-1b_PAM, OTU_1 G_Bacteroides

# comp1 loadings: OTU:WAZ ####
gc.cols <- as.character(cohort.cols$Color)
names(gc.cols) <- cohort.cols$Site

# combine lm.data for otu and lmx
lm.dat <- rbind(lm.dat.xz, lm.dat.yz)
lm.dat$Cohort <- substr(lm.dat$sample.id, start = 1, stop = 1)
lm.dat$Cohort <- ifelse(lm.dat$Cohort == "E", "ECD",
                           ifelse(lm.dat$Cohort == "F", "CAD", lm.dat$Cohort))


comp1.dat <- filter(lm.dat, X.var %in% c("OTU_1 G_Bacteroides", "OTU_2 G_Prevotella"), Y.var %in% c("MOM_AGE", "HAZ"))

# CAD/ECD together
ggplot(comp1.dat, aes(x = Y.value, y = X.value)) + 
  facet_wrap(~Y.var + X.var, scales = "free", nrow = 1) +
  geom_point(shape = 21, size = 3, aes(fill = Cohort)) +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_manual(values = gc.cols) +
  labs(x = "Host Factor", y = "OTU Abundance") +
  theme(legend.position = "bottom")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad_ecd/cad_ecd_comp1_plot.pdf", device = "pdf", dpi = 300, width = 12, height = 3)

# CAD/ECD separate
ggplot(comp1.dat, aes(x = Y.value, y = X.value)) + 
  facet_wrap(Cohort~Y.var + X.var, scales = "free", nrow = 2) +
  geom_point(shape = 21, size = 3, aes(fill = Cohort)) +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_manual(values = gc.cols) +
  labs(x = "Host Factor", y = "OTU Abundance") +
  theme(legend.position = "bottom")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad_ecd/cad_ecd_comp1_plot_separate.pdf", device = "pdf", dpi = 300, width = 12, height = 6)

# write results
otu.res.all <- x.y.lm.scaled(comp1.dat)
otu.res.all <- otu.res.all[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(otu.res.all, "ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad_ecd/graph_annotations/cad_ecd_otu_plot_annotations.csv") 

# cad only
cad.dat.otu <- filter(comp1.dat, Cohort == "CAD")
cad.res.otu <- x.y.lm.scaled(cad.dat.otu)
cad.res.otu <- cad.res.otu[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(cad.res.otu, "ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad_ecd/graph_annotations/cad_otuplot_annotations.csv") 

# ecd only
ecd.dat.otu <- filter(comp1.dat, Cohort == "ECD")
ecd.res.otu <- x.y.lm.scaled(ecd.dat.otu)
ecd.res.otu <- ecd.res.otu[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(ecd.res.otu, "ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad_ecd/graph_annotations/ecd_otuplot_annotations.csv") 

# END ####