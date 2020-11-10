#####################
# Nelly Amenyogbe
# Global Cohort Manuscript Analysis
# sPLS: microbiome-immune demographics: ECUADOR
#####################

#  In this script, we determine whether host demographic factors associated with microbiome-immune correlations among Ecuadorean children.

# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(mixOmics)
library(reshape2)
library(pheatmap)
# helper function
source("ms_global_cohort_analysis/scripts/functions/function_plot_taxa.R")
source("ms_global_cohort_analysis/scripts/functions/function_block_spls_loadings.R")
source("ms_global_cohort_analysis/scripts/functions/function_pairwise_lm.R")

# load data
meta.all <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_metadata.csv")
bfeed <- read.csv("Rdata/metadata/ecd_only/ECD_bfeed_duration.csv")
spls.dat <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_data.rds")

# colours that will be used for plotting
cohort.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_mb_cols.csv")

# prepare data for sPLS ####
# spls data
spls.data <- spls.dat$ECD

# metadata
spls.subjects <- rownames(spls.data$Xotu)
meta <- filter(meta.all, Cohort == "ECD", SUBJECT_ID %in% spls.subjects)

# add breastfeeding duration to metadata
colnames(bfeed)[1] <- "Microbiome_ID"
meta <- join(meta, bfeed, by = "Microbiome_ID")

meta$time.since.bfeed <- meta$age.stool.months - meta$agesuspendbrfeed.mo

# for negative values, change to 0
meta$time.since.bfeed <- ifelse(meta$time.since.bfeed < 0, 0, as.numeric(meta$time.since.bfeed))

# Only one infant was being breastfed at the time of sampling

rownames(meta) <- meta$SUBJECT_ID
meta <- meta[order(meta$SUBJECT_ID),]

length(which(rownames(meta) ==  rownames(spls.data$Xotu))) # 41 samples. all match.

# set factor levels to metadata
Z <- meta[,c("SEX", "GEST_AGE", "DELIVERY", "MOM_AGE", "WAZ", "WLZ", "HAZ", "time.since.bfeed", "agesuspendbrfeed.mo")]

colnames(Z)[which(colnames(Z) == "agesuspendbrfeed.mo")] <- "bfeed.duration"

Z$SEX <- ifelse(Z$SEX == "Male", 0,
                ifelse(Z$SEX == "Female", 1, Z$SEX))

Z$DELIVERY <- ifelse(Z$DELIVERY == "Caesarean", 0,
                         ifelse(Z$DELIVERY == "Vaginal", 1, Z$DELIVERY))

Z <- as.matrix(Z)

# ECD sPLS ####
X <- spls.data$Xotu
Y <- log10(spls.data$Ylmx)

# set parameters
dat <- list(mb = X, lmx = Y, meta = Z)

design = matrix(c(1,1,1,1,1,1,1,1,1), byrow = TRUE, nrow = 3)
diag(design) <- 0

res <- block.spls(X = dat,
                  indY = 3,
                  ncomp = 3,
                  keepX = list(mb = c(20, 20, 20), lmx = c(20, 20, 20), meta = c(8, 8, 8)),
                  design = design,
                  mode = "canonical")

# plot sPLS results ####

# Correlation circle
cor.circle <- plotVar(res,
                      col = c("blue", "red", "#006d2c"),
                      cex = c(2,2,2),
                      comp = c(1,2),
                      overlap = TRUE,
                      title = "ECD CLR",
                      plot = TRUE)


plotLoadings(res, comp = 1, title = "ECD") # WAZ and HAZ negatively associated, age.susp.bf with OTU_1404 G_Prevotella, MIP-1b_PGN

plotLoadings(res, comp = 2, title = "ECD") # DELIVERY, SEX associate with TNFa_R848, OTU_4848.G_Dorea

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
         main = "ECD sPLS_DM: COMP 1")

# exported figures/sPLS_mb_im_demographics/ecd/ecd_dm_hm_comp1.pdf

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
         main = "ECD sPLS_DM: COMP2")

# exported figures/sPLS_mb_im_demographics/ecd/ecd_dm_hm_comp2.pdf

# plot selected features: COMP1 ####

# explore top comp1 loadings ####
# DELIVERY, OTU_444 G_Bacteroides, OTU_71 G_Bifidobacterium, TNFa_PGN, TNFa_LPS

gc.cols <- as.character(cohort.cols$Color)
names(gc.cols) <- cohort.cols$Site

# combine lm.data for otu and lmx
lm.dat <- rbind(lm.dat.xz, lm.dat.yz)
lm.dat$Cohort <- substr(lm.dat$sample.id, start = 1, stop = 1)
lm.dat$Cohort <- ifelse(lm.dat$Cohort == "E", "ECD",
                        ifelse(lm.dat$Cohort == "F", "CAD", lm.dat$Cohort))


comp1.dat <- filter(lm.dat, X.var %in% c("OTU_1969 G_Prevotella", "OTU_1555 G_Prevotella"), Y.var %in% c("MOM_AGE"))

# plot MOM_AGE

ggplot(comp1.dat, aes(x = Y.value, y = X.value)) + 
  facet_wrap(~Y.var + X.var, scales = "free", nrow = 1) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(shape = 21, size = 3, aes(fill = Cohort)) +
  scale_fill_manual(values = gc.cols) +
  labs(x = "Maternal Age", y = "Normalized Feature") +
  theme(legend.position = "none")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/ecd/ecd_comp1_momage_plots.pdf", device = "pdf", dpi = 300, width = 6, height = 3)

# plot annotations
ma.res <- x.y.lm.scaled(comp1.dat)
ma.res <- ma.res[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(ma.res, "ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/ecd/graph_annotations/ecd_ma_plot_annotations.csv")


# COMP 2: DELIVERY ####
comp2.dat.dm <- filter(lm.dat, X.var %in% c("OTU_93 F_Lachnospiraceae", "OTU_4148 G_Dorea", "TNFa_R848", "IL-10_R848"), Y.var %in% c("DELIVERY"))

comp2.dat.dm$DeliveryMode <- ifelse(comp2.dat.dm$Y.value == 0, "Caesarean", "Vaginal")

# plot GA
ggplot(comp2.dat.dm, aes(x = DeliveryMode, y = X.value)) + 
  facet_wrap(~Y.var + X.var, scales = "free", nrow = 1) +
  #geom_smooth(method = "lm", color = "black") +
  geom_boxplot(outlier.size = NA) +
  geom_point(shape = 21, size = 3, aes(fill = Cohort)) +
  scale_fill_manual(values = gc.cols) +
  labs(x = "Delivery Mode", y = "Normalized Feature") +
  theme(legend.position = "none")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/ecd/ecd_comp2_dm_plots.pdf", device = "pdf", dpi = 300, width = 12, height = 3)

# get plot annotations
dm.res <- x.y.lm.scaled(comp2.dat.dm)
dm.res <- dm.res[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(dm.res, "ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/ecd/graph_annotations/ecd_dm_plot_annotations.csv")

# END ####