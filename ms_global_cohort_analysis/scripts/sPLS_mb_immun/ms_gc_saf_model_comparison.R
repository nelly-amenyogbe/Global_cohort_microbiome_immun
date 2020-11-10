####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# sPLS: Compare model results with/without saf
####################

# In this script, we will compare sPLS results with and without SAF included, to illustrate overlap of features selected at the expense of significant findings

# load packages
library(mixOmics)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(reshape2)
# helper function
source("ms_global_cohort_analysis/scripts/functions/function_pairwise_lm.R")

# load data
spls.ws <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/all_subjects_spls_res.rds")
spls.ns <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/blg_cad_ecd_spls_res.rds")

gc.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_mb_cols.csv")

# Comparison of features selected: OTU ####
ws.otu <- spls.ws$loadings$X.sum
ns.otu <- spls.ns$loadings$X.sum

colnames(ws.otu)[2:6] <- paste0("ws.", colnames(ws.otu)[2:6])
colnames(ns.otu)[2:6] <- paste0("ns.", colnames(ns.otu)[2:6])

# how many features overlap?
length(which(ws.otu$X.var %in% ns.otu$X.var)) # 61 OTUs

# set overlapping OTU data
otu.both <- ws.otu[which(ws.otu$X.var %in% ns.otu$X.var),]
otu.both <- join(otu.both, ns.otu, by = "X.var")

# plot selected OTUs ####

# loadings
ggplot(otu.both, aes(x = ws.loading, y = ns.loading)) +
  geom_point()

# plot position
ggplot(otu.both, aes(x = ws.n, y = ns.n)) +
  geom_point()

# plot percent sig
ggplot(otu.both, aes(x = ws.percent.x.sig, y = ns.percent.x.sig)) +
  geom_point() +
  geom_hline(yintercept = 30, color = "red") +
  geom_vline(xintercept = 30, color = "blue") +
  theme(axis.text = element_text(size = 12)) +
  labs(x = "all children: selected OTUs", y = "BLG,CAD,ECD children: selected OTUs")

length(which(otu.both$ns.percent.x.sig > 30)) # 22 OTUs pass the threshold for BLG/CAD/ECD
length(which(otu.both$ws.percent.x.sig > 30)) # 3 OTUs with SAF

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/model_comparison/percent_sig_comparison_OTU.pdf", device = "pdf", width = 5, height = 4)

# Comparison of features selected: LMX ####

ws.lmx <- spls.ws$loadings$Y.sum
ns.lmx <- spls.ns$loadings$Y.sum

colnames(ws.lmx)[2:6] <- paste0("ws.", colnames(ws.lmx)[2:6])
colnames(ns.lmx)[2:6] <- paste0("ns.", colnames(ns.lmx)[2:6])

# how many features overlap?
length(which(ws.lmx$Y.var %in% ns.lmx$Y.var)) # 67 Cytokines

# set overlapping OTU data
lmx.both <- ws.lmx[which(ws.lmx$Y.var %in% ns.lmx$Y.var),]
lmx.both <- join(lmx.both, ns.lmx, by = "Y.var")

# plot selected cytokines ####

# loadings
ggplot(lmx.both, aes(x = ws.loading, y = ns.loading)) +
  geom_point()

# plot position
ggplot(lmx.both, aes(x = ws.n, y = ns.n)) +
  geom_point()

# plot percent sig
ggplot(lmx.both, aes(x = ws.percent.y.sig, y = ns.percent.y.sig)) +
  geom_point() +
  geom_hline(yintercept = 30, color = "red") +
  geom_vline(xintercept = 30, color = "blue") +
  theme(axis.text = element_text(size = 12)) +
  labs(x = "all children: selected cyokines", y = "BLG,CAD,ECD children: selected cytokines")

length(which(lmx.both$ns.percent.y.sig > 30)) # 21 cytokines without saf
length(which(lmx.both$ws.percent.y.sig > 30)) # 11 cytokines with saf

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/model_comparison/percent_sig_comparison_LMX.pdf", device = "pdf", width = 5, height = 4)

# Conclude: number of cytokines co-varying roughly doubles with the exclusion of SAF children

# plot selected correlations ####
# Plot the three top relationships with/without SAF 

# get top loadings
ns.select.x <- filter(otu.both, ns.percent.x.sig > 30)
ns.select.y <- filter(lmx.both, ns.percent.y.sig > 30)

# for cytokines select first to PAM responses: MIP-1a_PAM, IL-1b_PAM
# for otus, select first Prevotella and first Bacteroides OTUs: OTU_18 G_Prevotella, OTU_3760 G_Bacteroides

# get lm results for selected features
otu.select <- c("OTU_1 G_Bacteroides", "OTU_3760 G_Bacteroides")
cyto.select <- c("MIP-1a_PAM", "IL-1b_PAM")

ns.lm.res <- spls.ns$lm.res
ws.lm.res <- spls.ws$lm.res

ns.lm.res.select <- filter(ns.lm.res, X.var %in% otu.select, Y.var %in% cyto.select)
ws.lm.res.select <- filter(ws.lm.res, X.var %in% otu.select, Y.var %in% cyto.select)

# plot selected: with SAF
lm.dat <- spls.ws$lm.dat

lm.dat.select <- filter(lm.dat, X.var %in% otu.select, Y.var %in% cyto.select)

# set colours
#group.cols <- filter(cohort.cols, Site != "SAF")
g.cols <- as.character(gc.cols$Color)
names(g.cols) <- as.character(gc.cols$Site)

g.shapes <- gc.cols$Shape
names(g.shapes) <- gc.cols$Site

# set cohorts
lm.dat.select$cohort <- substr(lm.dat.select$sample.id, start = 1, stop = 1)

lm.dat.select$Site <- ifelse(lm.dat.select$cohort == "B", "BLG",
                             ifelse(lm.dat.select$cohort == "F", "CAD",
                                    ifelse(lm.dat.select$cohort == "E", "ECD",
                                           ifelse(lm.dat.select$cohort == "S", "SAF", "uk"))))


# plot with SAF
ggplot(lm.dat.select, aes(x = Y.value, y = X.value)) +
  geom_point(aes(color = Site, shape = Site)) +
  scale_color_manual(values = g.cols) +
  scale_shape_manual(values = g.shapes) +
  geom_smooth(method = "lm", color = "black") +
  facet_grid(~ Y.var + X.var, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) +
  labs(x = "Cytokine concentration", y = "OTU abundance")

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/model_comparison/all_subjects_comparison_plots.pdf", device = "pdf", dpi = 300, width = 10, height = 2.5)


# plot without SAF
ggplot(filter(lm.dat.select, Site != "SAF"), aes(x = Y.value, y = X.value)) +
  geom_point(aes(color = Site, shape = Site)) +
  scale_color_manual(values = g.cols) +
  scale_shape_manual(values = g.shapes) +
  geom_smooth(method = "lm", color = "black") +
  facet_grid(~ Y.var + X.var, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) +
  labs(x = "Cytokine concentration", y = "OTU abundance")

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/model_comparison/blg_cad_ecd_comparison_plots.pdf", device = "pdf", dpi = 300, width = 10, height = 2.5)

# perform statistics to embed into plots
# with SAF
res.wsaf <- x.y.lm.scaled(lm.dat.select)
res.wsaf <- res.wsaf[,c("Estimate", "p.value", "X.var", "Y.var")]

# save results
#write.csv(res.wsaf, "ms_global_cohort_analysis/figures/sPLS-mb_immun/model_comparison/ms_gc_all_subjects_comparison_annotation.csv")

# Without
dat.nsaf <- filter(lm.dat.select, Site != "SAF")
res.nsaf <- x.y.lm.scaled(dat.nsaf)
res.nsaf <- res.nsaf[,c("Estimate", "p.value", "X.var", "Y.var")]

# save results
#write.csv(res.nsaf, "ms_global_cohort_analysis/figures/sPLS-mb_immun/model_comparison/ms_gc_blg_cad_ecd_comparison_annotation.csv")

# END ####