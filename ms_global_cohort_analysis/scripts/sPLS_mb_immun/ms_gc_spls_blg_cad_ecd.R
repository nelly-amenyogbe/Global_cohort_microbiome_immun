####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# sPLS: BLG_CAD_SAF
####################

# In this script, we perform sPLS for BLG, CAD, and ECD children. We then refine the model to only include features that vary across at least 30% of features in the respective data frame, and plot selected examples alongside their univariate correlation strength.

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

# run sPLS ####

# ALL spls ####
aX <- spls.dat$ALL$Xotu
aY <- log10(spls.dat$ALL$Ylmx)

# take out saf
X.rows <- substr(rownames(aX), start = 1, stop = 1)
Y.rows <-  substr(rownames(aY), start = 1, stop = 1)

aXs <- aX[-which(X.rows == "S"),]
aYs <- aY[-which(Y.rows == "S"),]

# check
length(which(rownames(aXs) ==  rownames(aYs))) # all good

all.spls <- spls(X = aXs,
                 Y = aYs,
                 ncomp = 3,
                 mode = "regression",
                 keepX = c(68, 68, 68),
                 keepY = c(67, 67, 67),
                 near.zero.var = TRUE)

# Plot ordinations
plotIndiv(all.spls, comp = c(1,2)) # comp1 - comp3 show dispersion

# Correlation circle
cor.circle <- plotVar(all.spls,
                      col = c("blue", "red"),
                      cex = c(2,2),
                      comp = c(1,2),
                      overlap = TRUE,
                      title = "All CLR",
                      plot = TRUE)

# make cor circle ####
all.spls.meta <- get.feature.data(all.spls, taxonomy)

all.comp.12 <- spls.cor.circle(all.spls,
                               x.meta = all.spls.meta$X.meta,
                               y.meta = all.spls.meta$Y.meta,
                               comps = c(1,2))

all.comp.12 + facet_grid(~Data)

all.comp.13 <- spls.cor.circle(all.spls,
                               x.meta = all.spls.meta$X.meta,
                               y.meta = all.spls.meta$Y.meta,
                               comps = c(1,3))

all.comp.13 + facet_grid(~Data) # nothing promising here

#### Regression figures ####

# get lm data
all.lm.dat <- lm.df.dataprep(aXs, aYs)

# run pairwise lm
all.lm.res <- x.y.lm.scaled(all.lm.dat)

# get feature summary
all.sum1 <- spls.feature.summary(all.spls,
                                 cor.results = all.lm.res,
                                 comp.no = "comp.1",
                                 pv.cutoff = 0.05)

all.sum2 <- spls.feature.summary(all.spls, 
                                 cor.results = all.lm.res,
                                 comp.no = "comp.2", 
                                 pv.cutoff = 0.05)

# Make heatmaps of selected features ####
all.lm.plot <- all.lm.res
all.lm.plot$Estimate <- ifelse(all.lm.res$p.value < 0.05, all.lm.res$Estimate, 0) # only features with significance over 0.05 are displayed

# comp1
all.sum1.sX <- all.sum1$X.sum
all.sum1.sX <- filter(all.sum1.sX, percent.x.sig > 30)

all.sum1.sY <- all.sum1$Y.sum
all.sum1.sY <- filter(all.sum1.sY, percent.y.sig > 30)

hm.est(all.lm.plot, 
       x.vars = all.sum1.sX$X.var,
       y.vars = all.sum1.sY$Y.var,
       title = "BLG CAD ECD",
       num.display = FALSE,
       spls.meta = all.spls.meta)

# saved as "figures/sPLS_mb_immun/blg_cad_ecd/ms_gc_spls_hm_blg_cad_ecd.png"

# comp2
all.sum2.sX <- all.sum2$X.sum
all.sum2.sX <- filter(all.sum2.sX, percent.x.sig > 30)

all.sum2.sY <- all.sum2$Y.sum
all.sum2.sY <- filter(all.sum2.sY, percent.y.sig > 30)

hm.est(all.lm.plot, 
       x.vars = all.sum2.sX$X.var,
       y.vars = all.sum2.sY$Y.var,
       title = "BLG CAD ECD comp2",
       num.display = FALSE,
       spls.meta = all.spls.meta)

# save comp1 result only.

# cor circle sig ####
x.feats <- c(all.sum1.sX$X.var)
y.feats <- c(all.sum1.sY$Y.var)

cor.circle.sig(all.spls,
               x.features = x.feats,
               y.features = y.feats,
               x.meta = all.spls.meta$X.meta,
               y.meta = all.spls.meta$Y.meta,
               comps = c(1,2))

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/blg_cad_ecd/blg_cad_ecd_cor_circle_sig.pdf", device = "pdf", dpi = 300, width = 7.2, height = 4.8)

# save relevant data ###
all.spls.save <- list("spls.obj" = all.spls, "otus.select" = x.feats, "lmx.select" = y.feats, "lm.res" = all.lm.res, "loadings" = all.sum1)

#saveRDS(all.spls.save, "ms_global_cohort_analysis/Rdata/R_export/spls_res/blg_cad_ecd_spls_res.rds")

#### Relevance Network ####
# Comp1
net1 <- spls.rn(lm.res <- all.lm.plot,
                x.vars <- all.sum1.sX$X.var,
                y.vars <- all.sum1.sY$Y.var,
                spls.meta <- all.spls.meta)

# layout
lo <- layout.fruchterman.reingold(net1$net, weights = (1 - 
                                                         abs(E(net1$net)$weight)))


# plot
p <- plot(net1$net, layout = lo, main = "BLG CAD ECD Comp1")

# plot legend: phyla

net1.fams <- filter(fam.cols, Family %in% net1$col.features)

legend(1.2,
       1.8,
       legend = c(as.character(net1.fams$Family)),
       col = c(as.character(net1.fams$col)),
       pch = 19,
       pt.cex = 2,
       title = "Family",
       bty = "n")

# cytokines
net1.stims <- filter(stim.cols, Stim %in% net1$col.features)

legend(1.3,
       -0.6,
       legend = c(as.character(net1.stims$Stim)),
       col = c(as.character(net1.stims$stim.cols)),
       pch = 15,
       pt.cex = 2,
       title = "Stimulus",
       bty = "n")

# edges 
color.edge <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))

edge.leg <- seq(from = -1, to = 1, length.out = 7)

edge.legend.cols <- color.edge(length(edge.leg))

color.legend(-1.8, 1.4, -1.6, 0.8,
             legend = round(edge.leg, digits = 1),
             rect.col = edge.legend.cols,
             cex = 0.8,
             gradient = "y")

# save network
#write_graph(net1$net, "ms_global_cohort_analysis/figures/sPLS-mb_immun/blg_cad_ecd/blg_cad_ecd_rn1.graphml", format = "graphml")

# plot selected relationships ####

# plot select OTUs in all ####
dat <- all.lm.dat
dat$cohort <- substr(dat$sample.id, start = 1, stop = 1)

dat$cohort <- ifelse(dat$cohort == "B", "BLG",
                     ifelse(dat$cohort == "F", "CAD", 
                            ifelse(dat$cohort == "E", "ECD", "uk")))

# filter for relevant features
dat.p <- filter(dat, X.var %in% c("OTU_2215 G_Prevotella", "OTU_1 G_Bacteroides"),
                Y.var %in% c("IP-10_PAM", "IL-10_PAM"))

group.cols <- filter(cohort.cols, Site != "SAF")
g.cols <- as.character(group.cols$Color)
names(g.cols) <- as.character(group.cols$Site)

g.shapes <- group.cols$Shape
names(g.shapes) <- group.cols$Site

p.all <- ggplot(dat.p, aes(x = Y.value, y = X.value)) +
  geom_point(aes(color = cohort, shape = cohort, size = cohort)) +
  geom_smooth(method = "lm", color = "black") +
  facet_grid(~ Y.var + X.var) +
  theme_bw() +
  scale_color_manual("Cohort", values = g.cols) +
  scale_shape_manual("Cohort", values = g.shapes) +
  scale_size_manual("Cohort", values = c(2, 2, 3)) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) +
  labs(x = "Cytokine concentration", y = "OTU abundance")

p.all

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/blg_cad_ecd/blg_cad_ecd_example_plots.pdf", device = "pdf", dpi = 300, width = 11, height = 2.3)

# ECD only
dat.pe <- filter(dat.p, cohort == "ECD")

p.ecd <- ggplot(dat.pe, aes(x = Y.value, y = X.value)) +
  geom_point(aes(color = cohort, shape = cohort, size = cohort)) +
  geom_smooth(method = "lm", color = "black") +
  facet_grid(~Y.var + X.var) +
  theme_bw() +
  scale_color_manual("Cohort", values = "#225ea8") +
  scale_shape_manual("Cohort", values = 18) +
  scale_size_manual("Cohort", values = 3) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) +
  labs(x = "Cytokine concentration", y = "OTU abundance")

p.ecd

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/blg_cad_ecd/all_ecd_example_plots.pdf", device = "pdf", dpi = 300, width = 11, height = 2.3)

# get lm results for plots
all.short.res <- x.y.lm.scaled(dat.p)
all.short.res <- all.short.res[,c("Estimate", "p.value", "X.var", "Y.var")]

#write.csv(all.short.res, "ms_global_cohort_analysis/figures/sPLS-mb_immun/blg_cad_ecd/example_plot_annotations.csv")

ecd.res <- x.y.lm.scaled(dat.pe)
ecd.res <- ecd.res[,c("Estimate", "p.value", "X.var", "Y.var")]
#write.csv(all.short.res, "ms_global_cohort_analysis/figures/sPLS-mb_immun/blg_cad_ecd/ecd_example_plot_annotations.csv")

# END ####