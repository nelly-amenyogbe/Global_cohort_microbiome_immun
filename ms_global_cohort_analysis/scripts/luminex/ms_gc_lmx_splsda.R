####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Luminex: sPLS-DA analysis
####################

# In this script, we will perform sPLS-DA analysis for the luminex cytokine dataset including all children, excluding South African children, to identify discriminatory cytokine signatures for Belgian, Canadian, and Ecuadorean infants.  Our previous findings illustrated that South African infants under-responded to PRR stimulation compared to all other site.  This dramatic effect failed to highlight more subtle differences among the other three sites. For this reason, we will focus on the other three sites here.

# Load packages
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mixOmics)
library(missForest)
library(tidyr)
# helper functions
source("ms_global_cohort_analysis/scripts/functions/matrix_missing_values_cleanup.R")
source("ms_global_cohort_analysis/scripts/functions/functions_mixomics_objects.R")

# load data
lmx <- readRDS("ms_global_cohort_analysis/Rdata/R_export/ms_gc_lmx_filtered_sPLS.rds")
lmx <- lmx$all

# results from Kruskal test
kw.res <- read.csv("ms_global_cohort_analysis/Rdata/R_export/ms_gc_lmx_cohort_kruskalwallis_res.csv")

# specify graphing attributes ####
# specify cohort colours
gc.cols <- c("BLG" = "#fec44f", "CAD" = "#ef3b2c", "ECD" = "#225ea8", "SAF" = "#41ab5d")
gc.shapes <- c("BLG" = 19, "CAD" = 17, "ECD" = 18, "SAF" = 15)
gc.size <- c(4,4,5,4)

# Prep lmx matrix for sPLS-DA ####
kw.res$padj <- p.adjust(kw.res$kw.p, method = "BH")
kw.res.sig <- filter(kw.res, padj < 0.05) # 62 cytokines

lmx.sig <- filter(lmx, cyto.stim %in% kw.res.sig$cyto.stim)

lmx.w <- dcast(lmx.sig, Sample + Site ~ Bead.Name + Stim, value.var = "final.concentration")

rownames(lmx.w) <- lmx.w$Sample
colnames(lmx.w)

lmx.mat <- lmx.w[,-which(colnames(lmx.w) %in% c("Sample", "Site"))]

# remove subjects with >10% missing data
sub.rm <- which.sub.rm(lmx.mat, 10) 
sub.rm # 4 removed: BLG1033, ECD1010, SAF1026, SAF1031
lmx.m.sr <- subj.rm(lmx.mat, 10)

# remove beads with >10% missing value
bead.rm <- beads.rm(lmx.m.sr, 10)
bead.rm # IL-23

lmx.rm <- remove.beads(lmx.m.sr, 10)

# replace missing values
lmx.mf <- missForest(lmx.rm)
lmx.final <- lmx.mf$ximp
rownames(lmx.final) <- rownames(lmx.rm)

# prepare metadata
meta.final <- lmx.w[which(lmx.w$Sample %in% rownames(lmx.final)),]

length(which(meta.final$Sample == rownames(lmx.final))) # 93. all good.

meta.final$Cohort <- ifelse(meta.final$Site == "Canada", "CAD",
                            ifelse(meta.final$Site == "Ecuador", "ECD",
                                   ifelse(meta.final$Site == "Belgium", "BLG", 
                                          ifelse(meta.final$Site == "South Africa", "SAF", meta.final$Site))))


# sPLSDA BLG, CAD, ECD ####
X <- log10(lmx.final)
Y <- meta.final$Cohort

ns <- filter(meta.final, Site != "South Africa")
X.ns <- X[rownames(X) %in% ns$Sample,]
length(which(rownames(X.ns) == ns$Sample))

Y.ns <- ns$Cohort

# let's try with a bunch more first
lmx.splsda.ns <- splsda(X = X.ns,
                        Y = Y.ns,
                        mode = "canonical",
                        keepX = c(30, 30, 30),
                        ncomp = 3)

err.ns <- perf(lmx.splsda.ns, validation = "Mfold", folds = 6, nrepeat = 100)
err.ns$error.rate.class
err.ns$error.rate.all

# plot error

er.dat <- get.error.dat(err.ns)

p.error <- ggplot(filter(er.dat, variable != "comp.4"), aes(x =  variable, y = value, group = group, color = group)) + 
  geom_line(size = 1.2) +
  theme_classic() +
  scale_color_manual(values = c(gc.cols, "Overall" = "#636363")) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "right") +
  labs(x = "Component", y = "Max error rate") +
  ylim(c(0,1.0))

p.error

# conclude:  Class error rate does not improve after 2nd component.
# plot error

# plot ordination
plotIndiv(lmx.splsda.ns, legend = TRUE, comp = c(1,2)) # CAD from BLG/ECD
plotIndiv(lmx.splsda.ns, legend = TRUE, comp = c(1,3)) # No discriminatory groupings
plotIndiv(lmx.splsda.ns, legend = TRUE, comp = c(2,3)) # BLG vs ECD/CAD

# The final model selects only ~5 features, as that is all that is needed to achieve a good error rate.  However, chances are you can select more features than that without sacrificing error, for the sake of biological exploration.  We will gage how many features to retain by capping it at roughly the number of features that still remained significant in KW results.  A crude, but acceptable strategy. 

# tune comp1-CAD
c1.cad <- plotLoadings(lmx.splsda.ns, comp = 1, contrib = "max", method = "mean")
c1.feat <- mc.get.features(lmx.splsda.ns, comp = "comp.1")
c1.feat <- filter(c1.feat, abs.loadings > 0)
colnames(c1.feat)[1] <- "cyto.stim"

c1.sig <- kw.res[,c("cyto.stim", "Belgium...Canada", "Canada...South.Africa", "Canada...Ecuador")]

c1.feat <- join(c1.feat, c1.sig, by = "cyto.stim")

c1.feat$index <- c(1:length(c1.feat$cyto.stim))

ggplot(c1.feat, aes(x = index, y = Belgium...Canada)) +
  geom_point() +
  geom_point(color = "red", aes(x=index, y = Canada...Ecuador)) +
  coord_flip() +
  geom_hline(yintercept = 0.05, color = "red")

# All remain signifiant until ~12 features in.  Choose 12.

# tune comp2-BLG
c2.blg <- plotLoadings(lmx.splsda.ns, comp = 2, contrib = "max", method = "mean")
c2.feat <- mc.get.features(lmx.splsda.ns, comp = "comp.2")
c2.feat <- filter(c2.feat, abs.loadings > 0)
colnames(c2.feat)[1] <- "cyto.stim"

c2.sig <- kw.res[,c("cyto.stim", "Belgium...Ecuador", "Belgium...South.Africa", "Belgium...Canada")]

c2.feat <- join(c2.feat, c2.sig, by = "cyto.stim")

c2.feat$index <- c(1:length(c2.feat$cyto.stim))

ggplot(c2.feat, aes(x = index, y = Belgium...Canada)) +
  geom_point() +
  geom_point(color = "red", aes(x = index, y = Belgium...Ecuador)) +
  coord_flip() +
  geom_hline(yintercept = 0.05, color = "red") # up to 14 features separate BLG from the rest

# re-run
lmx.splsda.ns <- splsda(X = X.ns,
                        Y = Y.ns,
                        mode = "canonical",
                        keepX = c(12, 14),
                        ncomp = 2)

# Plot sPLSDA-NS ####
ns.cols <- gc.cols[c(1:3)]
ns.shapes <- gc.shapes[c(1:3)]
ns.size <- gc.size[c(1:3)]

# plot indiv
p.ind <- plotIndiv(lmx.splsda.ns, legend = TRUE, comp = c(1,2), ind.names = FALSE)
p.dat <- p.ind$df
p.labs <- p.ind$graph$labels[2:3]

gg.ind.12 <- ggplot(p.dat, aes(x = x, y = y, colour = group, shape = group, size = group)) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  labs(x = p.labs[[1]], y = p.labs[[2]]) +
  scale_color_manual("Cohort", values = c(ns.cols)) +
  scale_shape_manual("Cohort", values = c(ns.shapes)) +
  scale_size_manual("Cohort", values = c(ns.size)) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right") 

gg.ind.12

#ggsave("ms_global_cohort_analysis/figures/luminex/gc_lmx_splsda12_ord.pdf", device = "pdf", dpi = 300, width = 7.46, height = 5.83)

# plot error NS ####
err.ns <- perf(lmx.splsda.ns, validation  = "Mfold", folds = 6)
er.dat.ns <- get.error.dat(err.ns)

p.error.ns <- ggplot(filter(er.dat.ns, variable != "comp.4"), aes(x =  variable, y = value, group = group, color = group)) + 
  geom_line(size = 1.2) +
  theme_classic() +
  scale_color_manual("Group", values = c(gc.cols, "Overall" = "#636363")) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right") +
  labs(x = "Component", y = "Max error rate") +
  ylim(c(0,1.0))

p.error.ns

#ggsave("ms_global_cohort_analysis/figures/luminex/gc_lmx_splsda_error.pdf", device = "pdf", dpi = 300, width = 4.8, height = 4)

# Plot selected features ####

# prepare lmx.ns to have stims and cytos separated

plot.lmx.xcyto <- function(df, Cyto.Stim, cols){
  
  df <- filter(df, cyto.stim %in% Cyto.Stim)
  
  p <- ggplot(df, aes(x = Bead.Name, y = log10(final.concentration), colour = Cohort)) +
    geom_boxplot(outlier.size = NA) +
    geom_point(size = 0.8, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1)) +
    theme_classic() +
    labs(x = "Cytokine", y = "Log10 - Concentration") +
    theme(axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1.0, vjust = 1.0),
          axis.text.y = element_text(size = 12, face = "bold"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)) +
    scale_color_manual(values = c(cols))
  
  p
  
}

plot.lmx.xstim <- function(df, Cyto.Stim, cols){
  
  df <- filter(df, cyto.stim %in% Cyto.Stim)
  
  p <- ggplot(df, aes(x = Stim, y = log10(final.concentration), colour = Cohort)) +
    geom_boxplot(outlier.size = NA) +
    geom_point(size = 0.8, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1)) +
    theme_classic() +
    labs(x = "Stimulus", y = "Log10 - Concentration") +
    theme(axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1.0, vjust = 1.0),
          axis.text.y = element_text(size = 12, face = "bold"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)) +
    scale_color_manual(values = c(cols))
  
  p
  
}

# prepare lmx data
lmx$Cohort <- ifelse(lmx$Site == "Canada", "CAD",
                     ifelse(lmx$Site == "Belgium", "BLG", 
                            ifelse(lmx$Site == "Ecuador", "ECD", 
                                   ifelse(lmx$Site == "South Africa", "SAF", lmx$Site))))

lmx.ns <- filter(lmx, Cohort!= "SAF")
lmx.ns$cyto.stim <- paste0(lmx.ns$Bead.Name, "_", lmx.ns$Stim)

# Plot comp1 CAD dotplots ####
cad.loadings <- mc.get.features(lmx.splsda.ns, "comp.1")
cad.loadings <- filter(cad.loadings, abs.loadings > 0)
cad.ld.up <- filter(cad.loadings, raw.loadings > 0)
cad.ld.down <- filter(cad.loadings, raw.loadings < 0)

# Canada lower ####

lmx.ns$Stim <- factor(lmx.ns$Stim, levels = c("Unstim", "PAM", "PGN", "pIC", "LPS", "R848"))

p.cad.low <- plot.lmx.xcyto(filter(lmx.ns), Cyto.Stim = c(cad.ld.up$feature), cols = gc.cols) + ggtitle("Canada: lower responses") + facet_grid(~Stim, scales = "free", space = "free")

p.cad.low

#ggsave("ms_global_cohort_analysis/figures/luminex/cad_features_down.pdf", device = "pdf", dpi = 300, height= 3.5, width = 7)

# Canada higher
p.cad.high <- plot.lmx.xcyto(filter(lmx.ns), Cyto.Stim = c(cad.ld.down$feature), cols = gc.cols) + ggtitle("Canada: higher responses") + facet_grid(~Stim, scales = "free", space = "free")

p.cad.high

#ggsave("ms_global_cohort_analysis/figures/luminex/cad_features_up.pdf", device = "pdf", dpi = 300, width = 3, height = 3.5)

# Plot BLG ####
blg.loadings <- mc.get.features(lmx.splsda.ns, "comp.2")
blg.loadings <- filter(blg.loadings, abs.loadings > 0)
blg.ld.up <- filter(blg.loadings, raw.loadings > 0)
blg.ld.down <- filter(blg.loadings, raw.loadings < 0)

# BLG higher
p.blg.up <- plot.lmx.xstim(lmx.ns, Cyto.Stim = c(blg.ld.up$feature), cols = gc.cols) + ggtitle("BLG: higher responses") + facet_grid(~Bead.Name, scales = "free", space = "free")

p.blg.up

#ggsave("ms_global_cohort_analysis/figures/luminex/blg_features_up.pdf", device = "pdf", dpi = 300, width = 4.4, height = 3.5)

# BLG lower
p.blg.down <- plot.lmx.xcyto(lmx.ns, Cyto.Stim = c(blg.ld.down$feature), cols = gc.cols) + ggtitle("BLG: lower responses") + facet_grid(~Stim, scales = "free", space = "free")

p.blg.down

#ggsave("ms_global_cohort_analysis/figures/luminex/blg_features_down.pdf", device = "pdf", dpi = 300, width = 8, height = 3.5)

# END ####