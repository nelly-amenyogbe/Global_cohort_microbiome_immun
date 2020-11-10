####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Microbiome: Beta diversity
####################

# In this script, we will perform an ordination of all samples via Bray-Curtis Distance, and visualize using non-metric multidimensional scaling.

# Load packages
library(phyloseq)
library(ggplot2)
library(plyr)
library(dplyr)

# load data
ps <- readRDS("ms_global_cohort_analysis/Rdata/human_raw_data/gc_physeq.rds")
meta <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_metadata.csv")

# add metadata to physeq
ps.meta <- data.frame(sample_data(ps))

meta <- filter(meta, Microbiome == "Y")
rownames(meta) <- meta$Microbiome_ID
meta.select <- meta[,c("Microbiome_ID", "ETHNICITY", "COUNTRY")]

# order to be same as metadata
meta.select <- meta.select[match(ps.meta$X.SampleID, meta.select$Microbiome_ID),]
length(which(meta.select$Microbiome_ID == rownames(ps.meta)))

sampledat <- sample_data(meta.select)
ps <- merge_phyloseq(ps, sampledat)

# create separate objects for each cohort
cad <- subset_samples(ps, Cohort == "CAD")
cad <- prune_taxa(taxa_sums(cad) > 0, cad)

saf <- subset_samples(ps, Cohort == "SAF")
saf <- prune_taxa(taxa_sums(saf) > 0, saf)

ecd <- subset_samples(ps, Cohort == "ECD")
ecd <- prune_taxa(taxa_sums(ecd) > 0, ecd)

blg <- subset_samples(ps, Cohort == "BLG")
blg <- prune_taxa(taxa_sums(blg) > 0, blg)

# Helper Functions ####
# plot.nmds
plot.nmds <- function(physeq, col, shp, sz){
  ord <- ordinate(physeq, method = "NMDS", distance = "bray")
  
  nmds <- plot_ordination(physeq, ord, type = "samples", color = col, shape = shp) + geom_point(aes(size = Cohort, alpha = 0.8)) + theme_classic()
}
# Bray-Curtis Distance NMDS
# physeq = phyloseq object, col = "color of dots", shp = "shape of dots", sz = size of dots (e.g. 4)

# get.nmds.data
get.nmds.data <- function(physeq){
  ord <- ordinate(physeq, method = "NMDS", distance = "bray")
  p <- plot_ordination(physeq, ord)
  dat <- p$data
  return(dat)
}

#### Ordinations ####
gc.cols <- c("BLG" = "#fec44f", "CAD" = "#ef3b2c", "ECD" = "#225ea8", "SAF" = "#41ab5d")
gc.shapes <- c("BLG" = 19, "CAD" = 17, "ECD" = 18, "SAF" = 15)
gc.size <- c(4,4,5,4)

# 1.  Cohort ####
gc.nmds.dat <- get.nmds.data(ps)

p.cohort <- ggplot(gc.nmds.dat, aes(x = NMDS1, y =NMDS2, color = Cohort, shape = Cohort, size = Cohort)) +
  geom_point(alpha = 0.9) +
  theme_classic() +
  scale_shape_manual(values = gc.shapes) +
  scale_colour_manual(values = gc.cols) + 
  scale_size_manual(values = gc.size) +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 14),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16))

p.cohort

#ggsave("ms_global_cohort_analysis/figures/microbiome/all_samples_nmds_ordination.pdf", device = "pdf", dpi = 300, width = 5.3, height = 4)

#### 2.  Belgium Country of Birth ####
blg.cols <- c("Belgium" = "#fec44f", "Congo" = "#2171b5", "Germany" = "#525252", "Mauritania" = "#006d2c", "Morocco" = "#e31a1c", "Uganda" = "#fc9272", "Canada" = "#bdbdbd", "Ecuador" = "#bdbdbd", "South Africa" = "#bdbdbd")
blg.shapes <- c(19,  19, 19, 19, 19, 19, 17, 18, 15)
blg.size <- c(4, 4, 4, 4, 4, 4, 4, 5, 4)


gc.nmds.dat$COUNTRY <- factor(gc.nmds.dat$COUNTRY, levels = c("Belgium", "Germany", "Congo", "Mauritania", "Morocco", "Uganda", "Canada", "Ecuador", "South Africa"))

# Reorder data frame
blg.dat <- filter(gc.nmds.dat, Cohort == "BLG")
nblg.dat <- filter(gc.nmds.dat, Cohort != "BLG")
blg.re <- rbind(nblg.dat, blg.dat)

blg.p <- blg.re[-which(is.na(blg.re$COUNTRY)),]

p.blg.cob <- ggplot(blg.p, aes(x = NMDS1, y = NMDS2, color = COUNTRY, shape = COUNTRY, size = COUNTRY)) +
  geom_point(alpha = 0.9) + 
  theme_classic() +
  scale_shape_manual(name = "Country of Birth", values = blg.shapes) +
  scale_colour_manual(name = "Country of Birth", values = blg.cols) + 
  scale_size_manual(name = "Country of Birth", values = blg.size) +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 14),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16))

p.blg.cob

#ggsave("ms_global_cohort_analysis/figures/microbiome/gc_birthcountry_nmds_ordination.pdf", device = "pdf", dpi = 300, width = 5.7, height = 4)

# END ####