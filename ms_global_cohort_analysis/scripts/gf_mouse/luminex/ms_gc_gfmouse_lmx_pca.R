####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Germ-Free Mouse: Luminex cytokine Principal Components Analysis
####################

# In this script, we will perform PCA for luminex cytokine data, to visualize the variance due to stimulus, cohort, and stool donor. We will also perform the wilcoxon test to determine which cytokine responses differed between SAF and CAD mice.

# load packages
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mixOmics)
# helper function
source("ms_global_cohort_analysis/scripts/functions/pca_labs_function.R")

# load data
lmx <- read.csv("ms_global_cohort_analysis/Rdata/mouse_raw_data/gc_gfmouse_lmx.csv")
stim.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_stim_cols.csv")
gc.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_mb_cols.csv")

# prepare data for plotting
# values of "-" denote values removed during QC due to low bead count.  Change this to NA.
lmx.p <- lmx

lmx.p$ConcFinal <- ifelse(lmx.p$ConcFinal == "-", NA, as.character(lmx.p$ConcFinal))

lmx.p$ConcFinal <- as.numeric(lmx.p$ConcFinal)

# annotate the cage ID ####
# mice CAD1-3 are cage1, 4-5 cage2, 7-9 cage 3.  Ditto for SAF mice.

lmx.p$mouse <- substr(lmx.p$MouseID, start = 4, stop = 4)

lmx.p$cage <- ifelse(lmx.p$mouse %in% c("1", "2", "3"), "1",
                     ifelse(lmx.p$mouse %in% c("4", "5", "6"), "2",
                            ifelse(lmx.p$mouse %in% c("7", "8", "9"), "3", lmx.p$mouse)))

lmx.p$cage.id <- paste0(lmx.p$GroupID, lmx.p$cage)

# do cursory plotting ####
ggplot(lmx.p, aes(x = Stim, y = log10(ConcFinal), color = GroupID)) +
  geom_boxplot(outlier.size = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.3)) +
  facet_grid(~Bead.Name, space = "free", scales = "free")

# No IL-23 produced at all. Remove from dataset.

lmx.p <- filter(lmx.p, Bead.Name != "IL-23")

ggplot(lmx.p, aes(x = Stim, y = log10(ConcFinal), color = GroupID)) +
  geom_boxplot(outlier.size = NA) +
  geom_point(size = 1, position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.9)) +
  facet_grid(~Bead.Name, scales = "free") + 
  theme_bw() 

#  Wilcoxon Test ####
# define variables to test: remove unstim values, IL-23 (not produced), IFN-a LPS (not produced).
lmx.p$stim.cyto <- paste(lmx.p$Stim, lmx.p$Bead.Name)

features <- lmx.p[,c("Stim", "Bead.Name", "stim.cyto")]
features <- unique(features)

features.stim <- filter(features, Bead.Name != "IL-23", Stim != "US")
features.stim <- filter(features.stim, stim.cyto != "LPS IFN-a")

stim.cyto <- unique(as.character(features.stim$stim.cyto))

# unstimulated features
#Only IL-10, IL-6, MIP-1b.  Rest were not detected.

unstim.features <- filter(features, Stim == "US", Bead.Name %in% c("IL-10", "IL-6", "MIP-1b"))

unstim.features <- as.character(unstim.features$stim.cyto)

all.features <- c(stim.cyto, unstim.features)

# run wilcox test ####

wilcox.res <- ldply(as.list(all.features), function(i){
  
  df.s <- filter(lmx.p, stim.cyto == i, GroupID == "SAF")
  df.c <-  filter(lmx.p, stim.cyto == i, GroupID == "CAD")
  w.res <- wilcox.test(df.s$ConcFinal, df.c$ConcFinal)
  
  p.v <- w.res$p.value
  stim.cyto <- i
  stim <- unique(as.character(df.s$Stim))
  cyto <- unique(as.character(df.s$Bead.Name))
  
  res <- data.frame(p.v, stim.cyto, stim, cyto)
  res
  
})

wilcox.res <- wilcox.res[order(wilcox.res$p.v),] 
wilcox.res # 6 sig US MIP-1b  LPS IL-10  R848 IL-10 R848 IL-6  LPS IFN-g  US IL-10

wilcox.res$p.adj <- p.adjust(wilcox.res$p.v, method = "BH") # 5 under 0.1, US MIP-1b  LPS IL-10  R848 IL-10 R848 IL-6  LPS IFN-g

# add significance stars
wilcox.res$sig.annot <- ifelse(wilcox.res$p.adj <= 0.0001, "****",
                               ifelse(wilcox.res$p.adj <= 0.001, "***",
                                      ifelse(wilcox.res$p.adj <= 0.01, "**",
                                             ifelse(wilcox.res$p.adj <= 0.05, "*",
                                                    ifelse(wilcox.res$p.adj <= 0.1, "+",
                                                           ifelse(wilcox.res$p.adj > 0.1, "ns", wilcox.res$p.adj))))))

wilcox.res$annot.size <- ifelse(wilcox.res$sig.annot %in% c("+", "ns"), 2, 3)


wilcox.sum <-  ddply(lmx.p,
                     c("Stim", "Bead.Name", "stim.cyto"),
                     summarise,
                     max = max(log10(ConcFinal), na.rm = TRUE))

wilcox.res <- join(wilcox.res, wilcox.sum, by = "stim.cyto")

# get significant features
wilcox.sig <- filter(wilcox.res, p.adj < 0.1)
wilcox.feats <- as.character(wilcox.sig$stim.cyto)

# make barplots
grp.cols <- as.character(gc.cols$Color)
names(grp.cols) <- gc.cols$Site
shapes <- c("SAF" = 22, "CAD" = 24)

# make boxplots: Sig features
ggplot(filter(lmx.p, stim.cyto %in% wilcox.feats), aes(x = Bead.Name, y = log10(ConcFinal), fill = GroupID, shape = GroupID)) +
  geom_boxplot(outlier.size = NA, alpha = 0.7) +
  geom_point(size = 2, color = "black", position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.6), aes(fill = GroupID, shape = GroupID)) +
  facet_grid(~Stim, scales = "free", space = "free") +
  theme_bw() +
  scale_fill_manual(values = grp.cols) +
  scale_shape_manual(values = shapes) +
  labs(x = "", y = "Log10 Concentration") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 1.0, size = 10, face = "bold"),
        strip.background = element_rect(fill = "#d9d9d9", color = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") +
  geom_text(data = filter(wilcox.res, stim.cyto %in% wilcox.feats), aes(x = Bead.Name, y = max + 0.2, label = sig.annot), inherit.aes = FALSE, size = 6)

#ggsave("ms_global_cohort_analysis/figures/gf_mouse/luminex/gf_sig_cytokines.pdf", device = "pdf", dpi = 300, width = 5, height = 4)

# make boxplots: all features

ggplot(lmx.p, aes(x = Stim, y = log10(ConcFinal), fill = GroupID)) +
  geom_boxplot(outlier.size = NA, alpha = 0.7) +
  geom_point(size = 1.2, color = "black", position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.6), aes(fill = GroupID, shape = GroupID)) +
  facet_grid(~Bead.Name, scales = "free", space = "free") +
  scale_shape_manual(values = shapes) +
  theme_bw() +
  scale_fill_manual(values = grp.cols) +
  labs(x = "", y = "Log10 Concentration") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 1.0, size = 10, face = "bold"),
        strip.background = element_rect(fill = "#d9d9d9", color = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") +
  geom_text(data = wilcox.res, aes(x = Stim, y = max + 0.3, label = sig.annot), inherit.aes = FALSE, size = 6)

#ggsave("ms_global_cohort_analysis/figures/gf_mouse/luminex/gf_all_cytokines.pdf", device = "pdf", dpi = 300, width = 9, height = 3.5)

# Conclude: Significant features after adjustment (using q < 0.1) are MIP-1b at baseline, IFN-g and IL-10 in response to LPS, and IL-6 in response to R848. 
#In addition, nominally significant features also include IL-10 and IFN-a in response to R848, and baseline levels of IL-10.

# PCA anaysis ####
# Zero-variance features should be excluded from PCA analysis.  To this end, only features tested under univariate conditions will apply here.

# prepare data matrix
lmx.cast <- dcast(lmx.p, Stim + MouseID + GroupID + FecalSample + cage.id ~ Bead.Name, value.var = "ConcFinal")

lmx.cast$unique <- paste(lmx.cast$MouseID, lmx.cast$Stim)
length(unique(lmx.cast$unique)) # these are all uniqe

# Mouse CAD2 is missing R848 data, and CAD3 is missing LPS data.
mouse.remove <- c("CAD3 LPS", "CAD2 R848")
lmx.cast.f <- lmx.cast[-which(lmx.cast$unique %in% mouse.remove),]

# separate to metadata and matrix
colnames(lmx.cast.f)

# metadata
lmx.meta <- lmx.cast.f[,c("unique", "Stim", "MouseID", "GroupID", "FecalSample", "cage.id")]

# matrix
lmx.mat <- lmx.cast.f[,c("unique", "IFN-g", "IL-10", "IL-6", "MIP-1b", "TNF-a")]
rownames(lmx.mat) <- lmx.mat$unique
lmx.mat <- lmx.mat[,-1]

# prepare lmx.meta


# make PCA ####
lmx.pca <- prcomp(log10(lmx.mat), scale. = TRUE)

pca.dat <- as.data.frame(lmx.pca$x)
pca.dat$unique <- rownames(pca.dat)
pca.dat <- join(pca.dat, lmx.meta, by = "unique")

# get pca axis labels
pca.labs <- get.pca.dat(lmx.pca)

# plot PCA ####
# need cols for CAD/SAF
pca.cols <- c("CAD US" = "#636363", "SAF US" = "#bdbdbd", "CAD LPS" = "#2c7fb8", "SAF LPS" = "#7fcdbb", "CAD R848" = "#d95f0e", "SAF R848" = "#fec44f")

pca.shapes <- c("CAD US" = 22, "SAF US" = 22, "CAD LPS" = 21, "SAF LPS" = 21, "CAD R848" = 24, "SAF R848" = 24)

pca.dat$site.stim <- paste(pca.dat$GroupID, pca.dat$Stim)
pca.dat$site.stim <- factor(pca.dat$site.stim, levels = c("CAD US", "SAF US", "CAD LPS", "SAF LPS", "CAD R848", "SAF R848"))

# comp 1 vs 2
pca.12 <- ggplot(pca.dat, aes(x = PC1, y = PC2, fill = site.stim, shape = site.stim)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = pca.labs$pc1, y = pca.labs$pc2) +
  scale_fill_manual(values = pca.cols) +
  scale_shape_manual(values = pca.shapes) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

pca.12

#ggsave("ms_global_cohort_analysis/figures/gf_mouse/luminex/gfmouse_lmx_pca_12.pdf", device = "pdf", dpi = 300, width = 4.5, height = 4.5)

# comp 1 vs 3
pca.13 <- ggplot(pca.dat, aes(x = PC1, y = PC3, fill = site.stim, shape = site.stim)) +
  geom_point(size = 3) +
  labs(x = pca.labs$pc1, y = pca.labs$pc3) +
  scale_fill_manual(values = pca.cols) +
  theme_classic()+
  scale_shape_manual(values = pca.shapes) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

pca.13
#ggsave("ms_global_cohort_analysis/figures/gf_mouse/luminex/gfmouse_lmx_pca_13.pdf", device = "pdf", dpi = 300, width = 4.5, height = 4.5)

# END ####

