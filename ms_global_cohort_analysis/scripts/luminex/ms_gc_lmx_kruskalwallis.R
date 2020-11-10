####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Luminex: Kruskal-Wallis test for cytokines that differ across cohort
####################

# In this script, we will perform the kruskal-wallis test for all stimulus-cytokines that passed the fligner-kileen test.  Cytokines that significantly differ across cohort will then be passed to sPLS-DA analysis, to determine the cohort-specific signatures in a multivariate space.

# load packages
library(plyr)
library(dplyr)
library(dunn.test)
library(reshape2)

# load data
lmx <- readRDS("ms_global_cohort_analysis/Rdata/R_export/ms_gc_lmx_filtered_sPLS.rds")

# use lmx data for all subjects
lmx.dat <- lmx$all

# Run Kruskal-Wallis ####
cytokines <- unique(as.character(lmx.dat$cyto.stim))

kw.res <- ldply(as.list(cytokines), function(i){
  
  df <- filter(lmx.dat, cyto.stim == i)
  
  ktest <- kruskal.test(final.concentration ~ Site, data = df)
  ktest.pv <- ktest$p.value
  
  dt <- dunn.test(df$final.concentration, df$Site, kw = TRUE, method = "bh")
  dt.groups <- dt$comparisons
  dt.p <- dt$P.adjusted
  
  res <- data.frame(dt.groups, dt.p)
  res$kw.p <- ktest.pv
  res$cyto.stim <- i
  
  return(res)
  
})

kw.res.cast <- dcast(kw.res, cyto.stim + kw.p ~ dt.groups, value.var = "dt.p") # place pairwise-comparisons in separate columns

# save results ####
#write.csv(kw.res.cast, "ms_global_cohort_analysis/Rdata/R_export/ms_gc_lmx_cohort_kruskalwallis_res.csv")

# END ####