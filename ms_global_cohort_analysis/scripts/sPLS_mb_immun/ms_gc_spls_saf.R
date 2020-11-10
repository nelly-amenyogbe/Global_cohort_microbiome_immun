####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# sPLS: SAF
####################

# In this script, we perform sPLS for SOUTH AFRICAN children. We then refine the model to only include features that vary across at least 30% of features in the respective data frame, and plot selected examples alongside their univariate correlation strength.

# load packages
library(mixOmics)
library(pheatmap)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(igraph)
library(plotrix)
# helper functions
source("ms_global_cohort_analysis/scripts/functions/function_pairwise_lm.R")
source("ms_global_cohort_analysis/scripts/functions/function_spls_feature_summary.R")
source("ms_global_cohort_analysis/scripts/functions/functions_spls_figures.R")

# load data
spls.dat <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_data.rds")

# colours that will be used for plotting
phy.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_phylum_cols.csv")
stim.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_stim_cols.csv")
fam.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_network_fam_cols.csv")
cohort.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_mb_cols.csv")

# taxonomy table
taxonomy <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_taxonomy_table.csv")
taxonomy$OTU.num <- taxonomy$X
taxonomy$names <- paste(taxonomy$OTU.num, taxonomy$Genus)

# SAF spls ####

X <- spls.dat$SAF$Xotu
Y <- log10(spls.dat$SAF$Ylmx)
n <- dim(X)[1]

# check
length(which(rownames(X) ==  rownames(Y))) # all good

spls.res <- spls(X = X,
                 Y = Y,
                 ncomp = 3,
                 mode = "regression",
                 keepX = c(n, n, n),
                 keepY = c(n, n, n),
                 near.zero.var = TRUE)

# Plot ordinations
plotIndiv(spls.res, comp = c(1,3)) # comp1-comp3 show dispersion

# Correlation circle
cor.circle <- plotVar(spls.res,
                      col = c("blue", "red"),
                      cex = c(2,2),
                      comp = c(1,3),
                      overlap = TRUE,
                      title = "SAF CLR",
                      plot = TRUE)

# make cor circle ####
spls.meta <- get.feature.data(spls.res, taxonomy)

# comp 1-2
cor.comp.12 <- spls.cor.circle(spls.res,
                               x.meta = spls.meta$X.meta,
                               y.meta = spls.meta$Y.meta,
                               comps = c(1,2))

cor.comp.12 + facet_grid(~Data)

# comp 1-3
cor.comp.13 <- spls.cor.circle(spls.res,
                               x.meta = spls.meta$X.meta,
                               y.meta = spls.meta$Y.meta,
                               comps = c(1,3))

cor.comp.13 + facet_grid(~Data)

#### Regression figures ####

# get lm data
lm.dat <- lm.df.dataprep(X, Y)

# run pairwise lm
lm.res <- x.y.lm.scaled(lm.dat)

# get feature summary
sum1 <- spls.feature.summary(spls.res,
                             cor.results = lm.res,
                             comp.no = "comp.1",
                             pv.cutoff = 0.05)

sum2 <- spls.feature.summary(spls.res, 
                             cor.results = lm.res,
                             comp.no = "comp.2", 
                             pv.cutoff = 0.05)

sum3 <- spls.feature.summary(spls.res, 
                             cor.results = lm.res,
                             comp.no = "comp.3", 
                             pv.cutoff = 0.05)

# Make heatmaps of selected features ####
lm.plot <- lm.res
lm.plot$Estimate <- ifelse(lm.res$p.value < 0.05, lm.res$Estimate, 0) # only features with significance over 0.05 are displayed

# comp1
sum1.sX <- sum1$X.sum
sum1.sX <- filter(sum1.sX, percent.x.sig > 30)

sum1.sY <- sum1$Y.sum
sum1.sY <- filter(sum1.sY, percent.y.sig > 30)

hm.est(lm.plot, 
       x.vars = sum1.sX$X.var,
       y.vars = sum1.sY$Y.var,
       title = "S Africa comp1",
       num.display = FALSE,
       spls.meta = spls.meta)

# comp2
sum2.sX <- sum2$X.sum
sum2.sX <- filter(sum2.sX, percent.x.sig > 30)

sum2.sY <- sum2$Y.sum
sum2.sY <- filter(sum2.sY, percent.y.sig > 30)

hm.est(lm.plot, 
       x.vars = sum2.sX$X.var,
       y.vars = sum2.sY$Y.var,
       title = "S Africa comp2",
       num.display = FALSE,
       spls.meta = spls.meta)

# comp3
sum3.sX <- sum3$X.sum
sum3.sX <- filter(sum3.sX, percent.x.sig > 30)

sum3.sY <- sum3$Y.sum
sum3.sY <- filter(sum3.sY, percent.y.sig > 30)

hm.est(lm.plot, 
       x.vars = sum3.sX$X.var,
       y.vars = sum3.sY$Y.var,
       title = "S Africa comp3",
       num.display = FALSE,
       spls.meta = spls.meta)

# hm 1:2
x.feats <- c(sum1.sX$X.var, sum2.sX$X.var, sum3.sX$X.var)
y.feats <- c(sum1.sY$Y.var, sum2.sY$Y.var, sum3.sY$Y.var)

hm.est(lm.plot, 
       x.vars = x.feats,
       y.vars = y.feats,
       title = "S Africa comp1:3",
       num.display = FALSE,
       spls.meta = spls.meta)

# exported as sPLS_mb_immun/saf/ms_gc_spls_hm_saf.png

# save relevant data ####

saf.spls.save <- list("spls.obj" = spls.res, "otu1.select" = sum1.sX$X.var, "otu2.select" = sum2.sX$X.var, "otu3.select" = sum3.sX$X.var, "lmx1.select" = sum1.sY$Y.var, "lmx2.select" = sum2.sY$Y.var, "lmx3.select" = sum3.sY$Y.var)

#saveRDS(saf.spls.save, "ms_global_cohort_analysis/Rdata/R_export/spls_res/saf_spls_res.rds")

# cor circle sig ####

# comp12
cor.circle.sig(spls.res,
               x.features = x.feats,
               y.features = y.feats,
               x.meta = spls.meta$X.meta,
               y.meta = spls.meta$Y.meta,
               comps = c(1,2))

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/saf/saf_cor_circle12_sig.pdf", device = "pdf", dpi = 300, width = 7.2, height = 4.8)

# comp13
cor.circle.sig(spls.res,
               x.features = x.feats,
               y.features = y.feats,
               x.meta = spls.meta$X.meta,
               y.meta = spls.meta$Y.meta,
               comps = c(1,3))

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/saf/saf_cor_circle13_sig.pdf", device = "pdf", dpi = 300, width = 7.2, height = 4.8)

#### Relevance Network ####
# RN Comp1:3 ####
net123 <- spls.rn(lm.res <- lm.plot,
                  x.vars <- x.feats,
                  y.vars <- y.feats,
                  spls.meta <- spls.meta)

# layout
lo <- layout.fruchterman.reingold(net123$net, weights = (1 - 
                                                           abs(E(net123$net)$weight)))

# plot
p <- plot(net123$net, layout = lo, main = "S.Africa Comp1:3")

# plot legend: phyla

net123.fams <- filter(fam.cols, Family %in% net123$col.features)


legend(-4,
       1.6,
       legend = c(as.character(net123.fams$Family)),
       col = c(as.character(net123.fams$col)),
       pch = 19,
       pt.cex = 2,
       title = "Family",
       bty = "n")

# cytokines
net123.stims <- filter(stim.cols, Stim %in% net123$col.features)

legend(-4,
       -0.4,
       legend = c(as.character(net123.stims$Stim)),
       col = c(as.character(net123.stims$stim.cols)),
       pch = 15,
       pt.cex = 2,
       title = "Stimulus",
       bty = "n")

# edges 
color.edge <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))

edge.leg <- seq(from = -1, to = 1, length.out = 7)

edge.legend.cols <- color.edge(length(edge.leg))

color.legend(1.8, 1.4, 1.6, 0.8,
             legend = round(edge.leg, digits = 1),
             rect.col = edge.legend.cols,
             cex = 0.8,
             gradient = "y")

# save network
#write_graph(net123$net, "ms_global_cohort_analysis/figures/sPLS-mb_immun/saf/saf_rn123.graphml", format = "graphml")

# plot selected relationships ####

dat <- lm.dat

dat.p <- filter(dat, X.var %in% c("OTU_1041 O_Clostridiales", "OTU_4077 F_Ruminococcaceae"),
                Y.var %in% c("IL-10_PGN", "IL-12p40_PGN"))

p.saf <- ggplot(dat.p, aes(x = Y.value, y = X.value)) +
  geom_point(size = 3, shape = 18, color = "#41ab5d") +
  geom_smooth(method = "lm", color = "black") +
  facet_grid(~ Y.var + X.var) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  labs(x = "Cytokine concentration", y = "OTU abundance")

p.saf

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/saf/saf_example_plots.pdf", device = "pdf", dpi = 300, width = 11, height = 2.3)

# stat summary
saf.res <- x.y.lm.scaled(dat.p)
saf.res <- saf.res[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(saf.res, "ms_global_cohort_analysis/figures/sPLS-mb_immun/saf/saf_example_plot_annotations.csv")

# END ####