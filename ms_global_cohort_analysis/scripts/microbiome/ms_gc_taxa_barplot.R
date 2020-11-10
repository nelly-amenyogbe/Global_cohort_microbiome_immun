####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Microbiome: Relative Abundance
####################

# In this script, we will prepare figures to visualize taxonomic composition of all samples in barplots, colored by the taxonomic representation of samples

# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
# helper function
source("ms_global_cohort_analysis/scripts/functions/function_plot_taxa.R")

# load data
ps <- readRDS("ms_global_cohort_analysis/Rdata/human_raw_data/gc_physeq.rds")

# trim data
ps.lcr <- prune_taxa(taxa_sums(ps) > 3, ps)

# transform to relative abundance data
ps.rab <- transform_sample_counts(ps.lcr, function(x) x/sum(x))

# Phylum-level data
phylum.df <- get_taxbar_data(ps.rab, num.taxa = 13, rank = "Phylum")

# Family-level data
family.df <- get_taxbar_data(ps.rab, num.taxa = 10, rank = "Family")

# Genus-level data
genus.df <- get_taxbar_data(ps.rab, num.taxa = 10, rank = "Genus")

# plot ####

# assign colours
cols.15 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0', "#525252","#969696")

cols.10 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a', "#525252",  "#969696")

# plot Phylum ####
phylum.df$Cohort <- factor(phylum.df$Cohort, levels = c("CAD", "BLG", "ECD", "SAF"))

# aggregate by Bacteroidetes rab
phylum.ag <- aggregate(Abundance ~ Sample + Cohort + Taxa, data = phylum.df, sum)
bac <- filter(phylum.ag, Taxa == "P_Bacteroidetes")
bac <- bac[order(bac$Abundance, decreasing = TRUE),]

phylum.ag$Sample <- factor(phylum.ag$Sample, levels = c(as.character(bac$Sample)))

p.phylum <- ggplot(phylum.ag, aes(x = Sample, y = Abundance, fill = Taxa)) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.text.y = element_text(face = "bold", size = 10),
        panel.spacing.x = unit(0.1, "lines")) +
  facet_grid(~Cohort, scale = "free", space = "free") + 
  scale_fill_manual(name = 'Phylum', values = c(cols.15)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Phylum Abundance") +
  labs(x = "Samples", y = "Relative Abundance")

p.phylum

# plot Family ####
family.df$Cohort <- factor(family.df$Cohort, levels = c("CAD", "BLG", "ECD", "SAF"))

# aggregate by Bacteroidetes rab
family.ag <- aggregate(Abundance ~ Sample + Cohort + Taxa, data = family.df, sum)

# CAD and BLG: by Bacteroidaceae
bac.f <- filter(family.ag, Cohort %in% c("CAD", "BLG"), Taxa == "F_Bacteroidaceae")
bac.f <- bac.f[order(bac.f$Abundance, decreasing = TRUE),]

# SAF and ECD: by Prevotellaceae
prev.f <- filter(family.ag, Cohort %in% c("SAF", "ECD"), Taxa == "F_Prevotellaceae")
prev.f <- prev.f[order(prev.f$Abundance),]

f.samples <- c(as.character(bac.f$Sample), as.character(prev.f$Sample))

family.ag$Sample <- factor(family.ag$Sample, levels = c(f.samples))

# plot
p.family <- ggplot(family.ag, aes(x = Sample, y = Abundance, fill = Taxa)) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 0.8),
        axis.text.y = element_text(face = "bold", size = 10),
        panel.spacing.x = unit(0.1, "lines")) +
  facet_grid(~Cohort, scale = "free", space = "free") + 
  scale_fill_manual(name = 'Phylum', values = c(cols.10)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Family Abundance") +
  labs(x = "Samples", y = "Relative Abundance")

p.family

# plot Genus ####
genus.df$Cohort <- factor(genus.df$Cohort, levels = c("CAD", "BLG", "ECD", "SAF"))

# aggregate by Bacteroides rab
genus.ag <- aggregate(Abundance ~ Sample + Cohort + Taxa, data = genus.df, sum)
# order by bacteroides abundance for CAD/BLG
bac.g <- filter(genus.ag, Taxa == "G_Bacteroides", Cohort %in% c("CAD", "BLG"))
# order by prevotella in ECD/SAF
prev.g <- filter(genus.ag, Taxa == "G_Prevotella", Cohort %in% c("ECD", "SAF"))

# combine datasets
bac.g <- bac.g[order(bac.g$Abundance, decreasing = TRUE),]
prev.g <- prev.g[order(prev.g$Abundance),]
g.samples <- c(as.character(bac.g$Sample), as.character(prev.g$Sample))

genus.ag$Sample <- factor(genus.ag$Sample, levels = c(g.samples))

p.genus <- ggplot(genus.ag, aes(x = Sample, y = Abundance, fill = Taxa)) + 
  geom_bar(stat = "identity") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.ticks.x = element_blank()) +
  facet_grid(~Cohort, scale = "free", space = "free") + 
  scale_fill_manual(name = 'Genus', values = c(cols.10)) +
  scale_y_continuous(expand = c(0,0)) +
  #ggtitle("Genus Abundance") +
  labs(x = "Samples", y = "Relative Abundance")

p.genus

#ggsave("ms_global_cohort_analysis/figures/microbiome/genus_barplot.pdf", device = "pdf",dpi = 300, height = 4.5, width = 9.53)

# END ####
