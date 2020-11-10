####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Germ-Free Mouse: Microbiome Ordination
####################

# In this script, we determine the Bray-Curtis distance between all microbiome samples, and visualize via Non-metric multidimensional scaling (NMDS). We then perform this analysis for the Ielum, Jejunum, and Feces separately, and perform the Adonis test to determine the variance explained by fecal donor country on microbial community composition.

# load packages
library(phyloseq)
library(ggplot2)
library(plyr)
library(dplyr)
library(vegan)

# load data
ps <- readRDS("ms_global_cohort_analysis/Rdata/mouse_raw_data/gc_gfmouse_phyloseq.rds")
gc.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_mb_cols.csv")

#### Set aesthetics for plotting ####

meta <- data.frame(sample_data(ps))

meta$Donor <- ifelse(meta$human.donor == "FM1005", "CAD1", 
                     ifelse(meta$human.donor == "FM1008", "CAD2", 
                            ifelse(meta$human.donor == "FM1011", "CAD3",
                                   ifelse(meta$human.donor == "HEU001G", "SAF1",
                                          ifelse(meta$human.donor == "HEU047G", "SAF2",
                                                 ifelse(meta$human.donor == "HEU089G", "SAF3", meta$human.donor))))))

sample_data(ps)$Donor <- meta$Donor

#### Ordinate samples ####
ord <- ordinate(ps, distance = "bray", method = "NMDS")

ord.dat <- plot_ordination(ps, ord, type = "samples")
ord.dat <- ord.dat$data

# ordinate tissues
ord.dat$species <- ifelse(ord.dat$species == "human", "Human",
                          ifelse(ord.dat$species == "mouse", "Mouse", ord.dat$species))

ord.dat$sample.type <- tolower(ord.dat$sample.type) # change to lower-case

# plot all samples
plot.ord.tissue <- ggplot(ord.dat, aes(x = NMDS1, y = NMDS2, color = species.sample, shape = cohort)) +
  geom_point(size = 3) +
  theme_classic() + 
  scale_shape_manual("Cohort", values = c("SAF" = 15, "CAD" = 17)) +
  scale_color_manual("Sample Type", values = c('#1b9e77','#d95f02','#7570b3','#e7298a')) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right")

plot.ord.tissue

#ggsave("ms_global_cohort_analysis/figures/gf_mouse/microbiome/gfmouse_ord_alltissues.pdf", device = "pdf", dpi = 300, width = 6, height = 4)

# each tissue separately
ord.dat$sample.type <- capitalize(ord.dat$sample.type)

plot.ord.cohort <- ggplot(filter(ord.dat, species != "Human"), aes(x = NMDS1, y = NMDS2, color = cohort, shape = cohort)) +
  geom_point(size = 2.5) +
  theme_bw() + 
  scale_shape_manual("Cohort", values = c("SAF" = 15, "CAD" = 17)) +
  scale_color_manual("Cohort", values = c("SAF" = "#41ab5d", "CAD" = "#ef3b2c")) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(size = 1.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right",
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "white", size = 1.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  facet_wrap(~sample.type, scales = "free")

plot.ord.cohort

#ggsave("ms_global_cohort_analysis/figures/gf_mouse/microbiome/gfmouse_ord_cohort.pdf", device = "pdf", dpi = 300, width = 8.67, height = 3)

# Adonis Test ####
# get R2 values for ordination
# get separate tissues
jej <- subset_samples(ps, sample.species == "Jejunum mouse")
il <- subset_samples(ps, sample.species == "Ileum mouse")
fec <- subset_samples(ps, sample.species == "Feces mouse")


# Jejunum
j.otu <- data.frame(t(otu_table(jej)))
j.meta <- data.frame(sample_data(jej))

j.adonis <- adonis(j.otu ~ cohort, data = j.meta, method = "bray")
j.adonis # R2 = 0.17, p = 0.019

# Ileum
i.otu <- data.frame(t(otu_table(il)))
i.meta <- data.frame(sample_data(il))

i.adonis <- adonis(i.otu ~ cohort, data = i.meta, method = "bray")
i.adonis # R2 = 0.21, p = 0.002

# Feces
f.otu <- data.frame(t(otu_table(fec)))
f.meta <- data.frame(sample_data(fec))

f.adonis <- adonis(f.otu ~ cohort, data = f.meta, method = "bray")
f.adonis # R2 = 0.20, p = 0.001

# END ####
