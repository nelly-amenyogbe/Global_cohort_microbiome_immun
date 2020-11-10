####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Germ-Free Mouse: DESeq2 for differentially abundant OTUs
####################

# In this script, we peform the DESeq2 test to identify differentially-abundant OTUs between mice gavaged with either CAD or SAF feces.  We then visualize these specific OTUs using boxplots.

# load packages
library(phyloseq)
library(ggplot2)
library(DESeq2)
library(plyr)
library(dplyr)
library(pheatmap)
library(reshape2)
# helper functions
source("ms_global_cohort_analysis/scripts/functions/functions_otu_transformations.R")

# load data
gf.physeq <- readRDS("ms_global_cohort_analysis/Rdata/mouse_raw_data/gc_gfmouse_phyloseq.rds")
phy.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_phylum_cols.csv")

# create human donor column with non-unique values
sample_data(gf.physeq)$donor <- as.numeric(sample_data(gf.physeq)$human.donor)

levs <- matrix(c(1,2,3,4,5,6,1,2,3,1,2,3), nrow = 6)
levs <- as.data.frame(levs)
colnames(levs)[1:2] <- c("donor", "donor.rep")

x <- sample_data(gf.physeq)
x <- data.frame(x)

x <- join(x, levs, by = "donor")
x$donor.rep <- factor(x$donor.rep)

sample_data(gf.physeq)$donor.rep <- x$donor.rep

# Filter OTUs ####
# 1.  Only OTUs present in all of the 7/9 mice per group
# 2.  Optional:  limit this to OTUs also detected in human samples ???

# prune data
# human
human <- subset_samples(gf.physeq, species == "human") #1848 taxa
human <- prune_taxa(taxa_sums(human) > 3, human)
human # 190 taxa with taxa sums greater than 3.

h.otu.tab <- data.frame(t(otu_table(human)))

# mouse
mouse <- subset_samples(gf.physeq, species == "mouse") #1848 taxa
mouse <- prune_taxa(taxa_sums(mouse) > 3, mouse)
mouse # 306 taxa after pruning

m.otu.tab <- data.frame(t(otu_table(mouse)))

# subset data to intestinal sites

# split by tissue type
feces <- subset_samples(mouse, sample.type == "Feces")
feces <- prune_taxa(taxa_sums(feces) > 3, feces)
feces # 253

ileum <- subset_samples(mouse, sample.type == "Ileum")
ileum <- prune_taxa(taxa_sums(ileum) > 3, ileum)
ileum # 167 taxa

jejunum <- subset_samples(mouse, sample.type == "Jejunum")
jejunum <- prune_taxa(taxa_sums(jejunum) > 3, jejunum)
jejunum # 142 taxa

physeqs <- list(feces, ileum, jejunum)
names(physeqs) <- c("feces", "ileum", "jejunum")

#### Get Deseq sigtabs ####

ds.lrt <- function(physeq){
  dds <- phyloseq_to_deseq2(physeq, ~ donor.rep + cohort)
  ds.lrt <- DESeq(dds, test = "LRT", fitType = "parametric", reduced = ~ donor.rep)
  res <- results(ds.lrt, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha),]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab),], "matrix"))
  sigtab$otu <- rownames(sigtab)
  return(sigtab)
}

sigtabs <- llply(physeqs, function(i){
  
  sigtab <- ds.lrt(i)
  
})

#### Plot OTU summaries ####
# set colours
phy.cols$Phylum <- gsub("P_", "p__", phy.cols$Phylum)
phy.cols$Phylum <- gsub("K_", "k__", phy.cols$Phylum)

phylum.cols <- as.character(phy.cols$phylum.cols)
names(phylum.cols) <- phy.cols$Phylum

# Plot sig OTUs: select only those in human samples
ds.genus.h <- function(sigtab, site){
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels = names(x))
  
  
  p <- ggplot(sigtab, aes( x = Genus, y = log2FoldChange, color = Phylum, shape = found.in.human)) +
    geom_point(size = 6, alpha = 0.8) +
    theme_bw() +
    ggtitle(site) +
    geom_hline(color = "red", yintercept = 0) +
    theme(plot.title = element_text(face = "bold")) +
    scale_color_manual(values = phylum.cols) +
    labs(x = "", y = "Log2FoldChange CAD:SAF") + 
    theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1.0, vjust = 1.0),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.border = element_rect(size = 1.5),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = "right",
          strip.text.x = element_text(size = 12, face = "bold"),
          strip.background = element_rect(fill = "white", size = 1.5))
  
  p
  
}

# feces
f.dat <- sigtabs$feces
f.dat$found.in.human <- ifelse(rownames(f.dat) %in% colnames(h.otu.tab), "TRUE", "FALSE")

f.gen.h <- ds.genus.h(f.dat, "Feces")
f.gen.h
#ggsave("ms_global_cohort_analysis/figures/gf_mouse/microbiome/mouse_feces_deseq.pdf", device = "pdf", dpi = 300, width = 6, height = 4.8)

# ileum
i.dat <- sigtabs$ileum
i.dat$found.in.human <- ifelse(rownames(i.dat) %in% colnames(h.otu.tab), "TRUE", "FALSE")

i.gen.h <- ds.genus.h(i.dat, "Ileum")
i.gen.h
#ggsave("ms_global_cohort_analysis/figures/gf_mouse/microbiome/mouse_ileum_deseq.pdf", device = "pdf", dpi = 300, width = 6, height = 4.8)

# feces
j.dat <- sigtabs$jejunum
j.dat$found.in.human <- ifelse(rownames(j.dat) %in% colnames(h.otu.tab), "TRUE", "FALSE")

j.gen.h <- ds.genus.h(j.dat, "Jejunum")
j.gen.h
#ggsave("ms_global_cohort_analysis/figures/gf_mouse/microbiome/mouse_jejunum_deseq.pdf", device = "pdf", dpi = 300, width = 6, height = 4.8)

#### overlap ####
all.overlap <- Reduce(intersect, list(rownames(sigtabs$feces), rownames(sigtabs$ileum), rownames(sigtabs$jejunum))) 
all.overlap # 5 OTUs

# Boxplot differentially abundant OTUs ####
all.dat <- rbind(f.dat, i.dat, j.dat)

all.h <- filter(all.dat, found.in.human == "TRUE")
sig.otus <- c(all.h$otu) 
sig.otus 
unique(sig.otus) # 14 OTUs, 21 unique

taxtab <- data.frame(tax_table(gf.physeq))
taxtab$OTU.num <- rownames(taxtab)
taxtab <- tax_table(as.matrix(taxtab))

gf.physeq <- merge_phyloseq(gf.physeq, taxtab)

# prepare OTU data for plotting
ps.sig <- subset_taxa(gf.physeq, OTU.num %in% sig.otus)# only 9 taxa!!

# get sig CLR data 
otu.sig <- data.frame(t(otu_table(ps.sig)))
sig.taxa <- data.frame(tax_table(ps.sig))

otu.tss <- TSS.transform(as.matrix(otu.sig))
colnames(otu.tss) <- paste(sig.taxa$OTU.num, sig.taxa$Genus)
otu.clr <- clr(otu.tss)

# get metadata
meta <- data.frame(sample_data(gf.physeq))

# set colours, shapes 
cols <- c("SAF" = "#41ab5d", "CAD" = "#ef3b2c")
shapes <- c("SAF" = 22, "CAD" = 24)

# melt otu table for ggplot
df.p <- as.data.frame(otu.clr)
df.p$Sample.ID <- rownames(df.p)

df.m <- melt(df.p, id.vars = "Sample.ID", variable.name = "OTU", value.name = "Abundance")
df.m <- join(df.m, meta, by = "Sample.ID")

# plot higher in saf
df.s <- filter(df.m, OTU %in% c("Otu0010 g__Prevotella", "Otu0185 f__Rikenellaceae_unclassified", "Otu0047 g__Odoribacter", "Otu0033 g__Alistipes"))

ggplot(df.s, aes(x = OTU, y = Abundance, fill = cohort, shape = cohort)) +
  geom_boxplot(outlier.size = NA, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.9, position = position_jitterdodge(dodge.width = 0.6)) +
  facet_grid(~sample.species) +
  theme_bw() +
  scale_fill_manual("Cohort", values = c(cols)) +
  scale_shape_manual("Cohort", values = c(shapes)) +
  xlab("") +
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
        legend.position = "top")

#ggsave("ms_global_cohort_analysis/figures/gf_mouse/microbiome/gf_saf_otus_boxplot.pdf", device = "pdf", dpi = 300, width = 7.8, height = 6)

# plot higher in cad
df.c <- filter(df.m, OTU %in% c("Otu0028 g__Clostridium", "Otu0166 g__Clostridium"))

ggplot(df.c, aes(x = OTU, y = Abundance, fill = cohort, shape = cohort)) +
  geom_boxplot(outlier.size = NA, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.9, position = position_jitterdodge(dodge.width = 0.6)) +
  facet_grid(~sample.species) +
  theme_bw() +
  scale_fill_manual("Cohort", values = c(cols)) +
  scale_shape_manual("Cohort", values = c(shapes)) +
  xlab("") +
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
        legend.position = "top")

#ggsave("ms_global_cohort_analysis/figures/gf_mouse/microbiome/gf_cad_otus_boxplot.pdf", device = "pdf", dpi = 300, width = 6.5, height = 5)
