####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# sPLS: hyperheometric test for selected cytokines
####################

# In this script, we perform the hypergeometric test to determine whether any cytokines were enriched for in the sPLS results for each cohort

# load packages
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)

# load data 
all.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/blg_cad_ecd_spls_res.rds")
blg.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/blg_spls_res.rds")
cad.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/cad_spls_res.rds")
ecd.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/ecd_spls_res.rds")
saf.spls <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_res/saf_spls_res.rds")

stim.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_stim_cols.csv")

# run hg test: all ####

# function for hyper test ####

lmx.phyper.test <- function(spls.res, features, cohort, comp){
  
  # list of all input features
  lmx.vars <- data.frame(colnames(spls.res$spls.obj$Y))
  colnames(lmx.vars)[1] <- "lmx.features"
  
  lmx.meta <- lmx.vars %>% separate(lmx.features,
                                    sep = "_",
                                    c("cytokine", "stim"))
  
  # list of selected features
  select <- data.frame(features)
  colnames(select)[1] <- "features"
  
  select.meta <- select %>% separate(features,
                                     sep = "_",
                                     c("cytokine", "stim"))
  
  # set the hypergeometric parameters
  cytokines <- unique(as.character(lmx.meta$cytokine))
  
  all.res <- ldply(as.list(cytokines), function(i){
    
    qt <- length(which(select.meta$cytokine == i)) # success
    mt <- length(which(lmx.meta$cytokine == i)) # number of feats available
    nt <- length(lmx.meta$cytokine) - mt # number of non-stim to choose from
    kt <- length(select.meta$cytokine) # number of features selected
    
    # test for over-representation
    res <- phyper(q = qt, m = mt, n = nt, k = kt, lower.tail = FALSE)
    
    cytokine <- i
    to.return <- c(res, cytokine)
    
  })
  
  colnames(all.res) <- c("p.value", "cytokine")
  all.res$cohort <- cohort
  all.res$component <- comp
  all.res$p.adj <- p.adjust(all.res$p.value, method = "BH")
  all.res
  
}

# do all tests
# all samples ####
all.res <- lmx.phyper.test(all.spls, all.spls$lmx.select, "all", "comp.1")


# BLG ####
# comp 1 and 2
blg1.res <- lmx.phyper.test(blg.spls, blg.spls$lmx1.select, "BLG", "comp.1")

blg2.res <- lmx.phyper.test(blg.spls, blg.spls$lmx2.select, "BLG", "comp.2")

# CAD ####
# comp1 and 2
cad1.res <- lmx.phyper.test(cad.spls, cad.spls$lmx1.select, "CAD", "comp.1")

cad2.res <- lmx.phyper.test(cad.spls, cad.spls$lmx2.select, "CAD", "comp.2")

# ECD ####
# comp 1 and 3
ecd1.res <- lmx.phyper.test(ecd.spls, ecd.spls$lmx1.select, "ECD", "comp.1")

ecd3.res <- lmx.phyper.test(ecd.spls, ecd.spls$lmx3.select, "ECD", "comp.3")

# SAF ####
# comps 1, 2, and 3
saf1.res <- lmx.phyper.test(saf.spls, saf.spls$lmx1.select, "SAF", "comp.1")

saf2.res <- lmx.phyper.test(saf.spls, saf.spls$lmx2.select, "SAF", "comp.2")

saf3.res <- lmx.phyper.test(saf.spls, saf.spls$lmx3.select, "SAF", "comp.3")

# combine and plot #### 
combined.res <- rbind(all.res, blg1.res, blg2.res, cad1.res, cad2.res, ecd1.res, ecd3.res, saf1.res, saf2.res, saf3.res)

combined.res$sig <- ifelse(combined.res$p.adj < 0.05, "sig", "ns")
combined.res$log.p <- log10(as.numeric(combined.res$p.adj))

combined.s <- combined.res[,c("cohort", "cytokine", "sig")]
combined.s <- unique(combined.s)

combined.sig <- filter(combined.s, sig == "sig")

# plot
cols <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666', '#e41a1c','#377eb8','#4daf4a','#984ea3', "#525252")

combined.res <- filter(combined.res, cytokine != "MDC")

ggplot(combined.res, aes(x = cytokine, y = abs(log.p), fill = cytokine)) +
  geom_hline(yintercept = abs(log10(0.05)), color = "#525252", size = 0.8) +
  geom_bar(stat="identity", position = "dodge", color = "black") +
  coord_flip() +
  theme_bw() +
  facet_grid(cohort~component, space = "free") +
  scale_alpha_manual("Significance", values = c(0.3, 1)) +
  scale_fill_manual("Cytokine", values = cols) +
  theme(strip.text = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "white")) +
  labs(x = "Cytokine", y = "Absolute log10(adjusted p-value)")

# summarize sig y/n
res.sum <- combined.res[,c("cytokine", "cohort", "component", "sig")]
res.sum <- unique(res.sum)
res.sum.cast <- dcast(res.sum, cytokine + cohort ~ component, value.var = "sig")
res.sum.cast$sig.yn <- ifelse(res.sum.cast$comp.1 == "sig" |
                                res.sum.cast$comp.2 == "sig" |
                                res.sum.cast$comp.3 == "sig", "sig", "ns")

res.sum.cast$sig.yn <- ifelse(is.na(res.sum.cast$sig.yn), "ns", as.character(res.sum.cast$sig.yn))

res.sum.cast <- res.sum.cast[,-c(3:5)]
res.sum.cast <- unique(res.sum.cast)

res.sum.cast$cohort <- factor(res.sum.cast$cohort, levels = c("SAF", "ECD", "CAD", "BLG", "all"))

# plot
ggplot(filter(res.sum.cast, sig.yn == "sig"), aes(x = cohort, fill = cytokine)) +
  geom_bar(color = "black", stat = "count") +
  coord_flip() +
  theme_classic() +
  scale_fill_manual("Cytokine", values = cols) +
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "", y = "")

#ggsave("ms_global_cohort_analysis/figures/sPLS-mb_immun/gc_cyto_phypher_sig.pdf", device = "pdf", dpi = 300, width = 4.3, height = 2)

# END ####