#####################
# Nelly Amenyogbe
# Global Cohort Manuscript Analysis
# sPLS: microbiome, immune, and host demographic data for Canadian children
#####################

# In this script, we determine whether host demographic factors associated with microbiome-immune correlations among Canadian children.

# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(mixOmics)
library(reshape2)
library(pheatmap)
library(cowplot)
# helper function
source("ms_global_cohort_analysis/scripts/functions/function_plot_taxa.R")
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
spls.data <- spls.dat$CAD

# metadata
spls.subjects <- rownames(spls.data$Xotu)
meta <- filter(meta.all, SUBJECT_ID %in% spls.subjects)

rownames(meta) <- meta$SUBJECT_ID
meta <- meta[order(meta$SUBJECT_ID),]

length(which(rownames(meta) ==  rownames(spls.data$Xotu))) # 19 samples. all match.

# add breastfeeding data to metadata ####
# create age in months
meta$age <- ifelse(is.na(meta$age.stool.months), meta$age.months, meta$age.stool.months)

meta <- join(meta, bfeed[,c("SUBJECT_ID", "ever.breastfed", "breastnow", "agesuspendbrfeed.mo")], by = "SUBJECT_ID")

# calculate time since breastfeeding
meta$time.since.bfeed <- ifelse(is.na(meta$breastnow) | meta$breastnow == "N", meta$age - meta$agesuspendbrfeed.mo, 
                                ifelse(meta$breastnow == "Y", 0, NA))

# change negative values to 0
meta$time.since.bfeed <- ifelse(meta$time.since.bfeed < 0, 0, meta$time.since.bfeed)

meta$bfeed.duration <- ifelse(is.na(meta$breastnow) | meta$breastnow == "N", meta$agesuspendbrfeed.mo,
                              ifelse(meta$breastnow == "Y", meta$age, NA))


# set factor levels to metadata
Z <- meta[,c("SEX", "GEST_AGE", "DELIVERY", "MOM_AGE", "WAZ", "WLZ", "HAZ", "bfeed.duration", "time.since.bfeed")]

Z$SEX <- ifelse(Z$SEX == "Male", 0,
                ifelse(Z$SEX == "Female", 1, Z$SEX))

Z$DELIVERY <- ifelse(Z$DELIVERY == "Caesarean", 0,
                     ifelse(Z$DELIVERY == "Vaginal", 1, Z$DELIVERY))

Z <- as.matrix(Z)
rownames(Z) <- meta$SUBJECT_ID

# CAD sPLS ####
X <- spls.data$Xotu
Y <- log10(spls.data$Ylmx)

# set parameters
dat <- list(mb = X, lmx = Y, meta = Z)

design = matrix(c(1,1,1,1,1,1,1,1,1), byrow = TRUE, nrow = 3)
diag(design) <- 0

res <- block.spls(X = dat,
                  indY = 3,
                  ncomp = 3,
                  keepX = list(mb = c(19, 19, 19), lmx = c(19, 19, 19), meta = c(7, 7, 7)),
                  design = design,
                  mode = "canonical")

# Correlation circle
cor.circle <- plotVar(res,
                      col = c("blue", "red", "#006d2c"),
                      cex = c(2,2,2),
                      comp = c(1,2),
                      overlap = TRUE,
                      title = "CAD CLR",
                      plot = TRUE)


plotLoadings(res, comp = 1, title = "CAD") # SEX, DELIVERY with IFNa2_R848, OTU_2217 F_Lachnospiraceae

plotLoadings(res, comp = 2, title = "CAD") # GEST_AGE, LAZ with IFNg_pIC, OTU_25 F_Rikenellaceae

# plot hm sig findings ####
##### Regression: Z:Z ####
# OTU:meta
# get lm data
lm.dat.xz <- lm.df.dataprep(X, Z)
lm.dat.yz <- lm.df.dataprep(Y, Z)
lm.dat <- rbind(lm.dat.xz, lm.dat.yz)

# run pairwise lm
lm.res.xz <- x.y.lm.scaled(lm.dat.xz)
lm.res.yz <- x.y.lm.scaled(lm.dat.yz)

lm.res <- rbind(lm.res.xz, lm.res.yz)

# get loadings
loadings.x <- data.frame(res$loadings$mb)
loadings.y <- data.frame(res$loadings$lmx)
loadings.z <- data.frame(res$loadings$meta)

# select vars to plot
x.ld <- block.spls.get.features(loadings.x, "comp.1")
y.ld <- block.spls.get.features(loadings.y, "comp.1")

# get data to plot
# get features
x.feats <- x.ld$feature[1:10]
y.feats <- y.ld$feature[1:10]
feats <- c(x.feats, y.feats)

# prepare hm res
lm.res.select <- filter(lm.res, X.var %in% feats, Y.var %in% c("MOM_AGE", "HAZ", "WAZ", "SEX", "DELIVERY", "WLZ", "GEST_AGE", "time.since.bfeed", "bfeed.duration"))
lm.res.select$Estimate <- ifelse(lm.res.select$p.value < 0.1, lm.res.select$Estimate, 0)

# prepare HM data
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
         main = "CAD sPLS_DM: COMP 1")

# exported figures/sPLS_mb_im_demographics/cad/cad_dm_hm_comp1.pdf

# COMP 2 HM ####
# select vars to plot
x.ld <- block.spls.get.features(loadings.x, "comp.2")
y.ld <- block.spls.get.features(loadings.y, "comp.2")

# get data to plot
# get features
x.feats <- x.ld$feature[1:10]
y.feats <- y.ld$feature[1:10]
feats <- c(x.feats, y.feats)

# prepare hm res
lm.res.select <- filter(lm.res, X.var %in% feats, Y.var %in% c("MOM_AGE", "HAZ", "WAZ", "SEX", "DELIVERY", "WLZ", "GEST_AGE", "time.since.bfeed", "bfeed.duration"))
lm.res.select$Estimate <- ifelse(lm.res.select$p.value < 0.05, lm.res.select$Estimate, 0)

# prepare HM data
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
         main = "CAD sPLS_DM: COMP2")

# exported figures/sPLS_mb_im_demographics/cad/cad_dm_hm_comp2.pdf


# COMP 3 HM ####
# select vars to plot
x.ld <- block.spls.get.features(loadings.x, "comp.3")
y.ld <- block.spls.get.features(loadings.y, "comp.3")

# get data to plot
# get features
x.feats <- x.ld$feature[1:10]
y.feats <- y.ld$feature[1:10]
feats <- c(x.feats, y.feats)

# prepare hm res
lm.res.select <- filter(lm.res, X.var %in% feats, Y.var %in% c("MOM_AGE", "HAZ", "WAZ", "SEX", "DELIVERY", "WLZ", "GEST_AGE", "time.since.bfeed", "bfeed.duration"))
lm.res.select$Estimate <- ifelse(lm.res.select$p.value < 0.05, lm.res.select$Estimate, 0)

# prepare HM data
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
         main = "CAD sPLS_DM: COMP3")

# exported figures/sPLS_mb_im_demographics/cad/cad_dm_hm_comp3.pdf

# plot selected features ####

# prepare plotting data ####
gc.cols <- as.character(cohort.cols$Color)
names(gc.cols) <- cohort.cols$Site

# combine lm.data for otu and lmx
lm.dat <- rbind(lm.dat.xz, lm.dat.yz)
lm.dat$Cohort <- substr(lm.dat$sample.id, start = 1, stop = 1)
lm.dat$Cohort <- ifelse(lm.dat$Cohort == "E", "ECD",
                        ifelse(lm.dat$Cohort == "F", "CAD", lm.dat$Cohort))

# COMP1: DELIVERY ####

# explore top comp1 loadings ####
# DELIVERY, OTU_444 G_Bacteroides, OTU_71 G_Bifidobacterium, TNFa_PGN, TNFa_LPS

comp1.dat <- filter(lm.dat, X.var %in% c("OTU_444 G_Bacteroides", "OTU_878 F_Lachnospiraceae", "TNFa_PGN", "TNFa_LPS"), Y.var %in% c("DELIVERY"))

# plot DELIVERY
# 0 = Caesarean, 1 = Vaginal
comp1.dat$DeliveryMode <- ifelse(comp1.dat$Y.value == 0, "Caesarean", "Vaginal")

ggplot(comp1.dat, aes(x = DeliveryMode, y = X.value)) + 
  facet_wrap(~Y.var + X.var, scales = "free", nrow = 1) +
  #geom_smooth(method = "lm", color = "black") +
  geom_boxplot(outlier.size = NA) +
  geom_point(shape = 21, size = 3, aes(fill = Cohort)) +
  scale_fill_manual(values = gc.cols) +
  labs(x = "Delivery Mode", y = "Normalized Feature") +
  theme(legend.position = "none")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad/cad_delivery_plots.pdf", device = "pdf", dpi = 300, width = 12, height = 3)


# plot annotations
deliv.res <- x.y.lm.scaled(comp1.dat)
deliv.res <- deliv.res[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(deliv.res, "ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad/graph_annotations/cad_dm_plot_annotations.csv")

# get annotations

# COMP 1: SEX ####
comp1.dat.s <- filter(lm.dat, X.var %in% c("OTU_878 F_Lachnospiraceae", "OTU_1705 O_Clostridiales", "IFNa2_R848", "IL-12p70_R848", "IL-23_LPS"), Y.var %in% c("SEX"))

comp1.dat.s$Sex <- ifelse(comp1.dat.s$Y.value == 0, "Male", "Female")

# plot SEX
ggplot(comp1.dat.s, aes(x = Sex, y = X.value)) + 
  facet_wrap(~Y.var + X.var, scales = "free", nrow = 1) +
  #geom_smooth(method = "lm", color = "black") +
  geom_boxplot(outlier.size = NA) +
  geom_point(shape = 21, size = 3, aes(fill = Cohort)) +
  scale_fill_manual(values = gc.cols) +
  labs(x = "Sex", y = "Normalized Feature") +
  theme(legend.position = "none")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad/cad_sex_plots.pdf", device = "pdf", dpi = 300, width = 12, height = 3)

# plot annotations
sex.res <- x.y.lm.scaled(comp1.dat.s)
sex.res <- sex.res[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(sex.res, "ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad/graph_annotations/cad_sex_plot_annotations.csv")

# COMP 2: breastfeeding ####
comp2.dat <- filter(lm.dat, X.var %in% c("OTU_948 G_Lachnospira", "OTU_30 G_Roseburia", "OTU_58 G_Faecalibacterium"), Y.var %in% c("time.since.bfeed"))

# plot bf
ggplot(comp2.dat, aes(x = Y.value, y = X.value)) + 
  facet_wrap(~Y.var + X.var, scales = "free", nrow = 1) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(shape = 21, size = 3, aes(fill = Cohort)) +
  scale_fill_manual(values = gc.cols) +
  labs(x = "Months since weaning", y = "OTU Abundance") +
  theme(legend.position = "bottom")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad/cad_comp2_bfeed_plots.pdf", device = "pdf", dpi = 300, width = 12, height = 3)

# get plot annotations
bf.res <- x.y.lm.scaled(comp2.dat)
bf.res <- bf.res[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(bf.res, "ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad/graph_annotations/comp2_bf_plot_annotations.csv")

# COMP3: WLZ/WAZ ####

comp3.dat <- filter(lm.dat, X.var %in% c("OTU_16 G_Lachnospira", "OTU_1690 F_Rikenellaceae"), Y.var %in% c("WAZ", "WLZ"))

ggplot(comp3.dat, aes(x = Y.value, y = X.value)) + 
  facet_wrap(~Y.var + X.var, scales = "free", nrow = 1) +
  geom_smooth(method = "lm", color = "black") +
  #geom_boxplot(outlier.size = NA) +
  geom_point(shape = 21, size = 3, aes(fill = Cohort)) +
  scale_fill_manual(values = gc.cols) +
  labs(x = "Host Factor", y = "Normalized Feature") +
  theme(legend.position = "none")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad/cad_waz_plots.pdf", device = "pdf", dpi = 300, width = 12, height = 3)

waz.res <- x.y.lm.scaled(comp3.dat)
waz.res <- waz.res[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(waz.res, "ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad/graph_annotations/cad_waz_plot_annotations.csv")

# END ####