####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Luminex: prepare luminex data for sPLS and sPLS-DA analyses
####################

# In this script, we will prepare the luminex data for use with sPLS and sPLS-DA by selecting only stimulus-cytokine combinations that were significantly produced above baseline.  These will be determine using the Fliger-Kileen test.  Unstimulated values will only be selected if over 70% of values were above the assay threshold.

# load packages
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

# load data
lmx <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_luminex.csv")

# flag anaytes with zeroes
zero <- filter(lmx, final.concentration == 0)
unique(zero$Bead.Name) # IL-23 only

# These will be replaced with 0.01 * the minimum concentration
il23 <- filter(lmx, Bead.Name == "IL-23", final.concentration > 0)
il23.min <- min(il23$final.concentration) # 4
# change zero-values to 0.04
lmx$final.concentration <- ifelse(lmx$final.concentration == 0, 0.04, lmx$final.concentration)

# Prepare data for fligner-kileen test ####
lmx$subject.bead <- paste(lmx$Sample, lmx$Bead.Name)
lmx <- ddply(lmx, .(subject.bead), transform, unstim.value = final.concentration[Stim == "Unstim"])
lmx$cyto.stim <- paste0(lmx$Bead.Name, "_", lmx$Stim)

#### Run Fligner test ####
run.fligner <- function(cohorts){
  x <- filter(lmx, Site %in% cohorts, Stim != "Unstim")
  xm <- melt(x, id.vars = "cyto.stim", measure.vars = c("final.concentration", "unstim.value"))
  cytokines <- unique(as.character(xm$cyto.stim))
  
  df.result <- ldply(as.list(cytokines), function(i){
    
    df.test <- filter(xm, cyto.stim == i)
    test <- fligner.test(df.test$value, df.test$variable)
    
    p.value <- test$p.value
    cyto.stim <- i
    res <- data.frame(cyto.stim, p.value)
    res
  })
  df.result <- df.result[order(df.result$p.value),]
  df.result$p.adj <- 5 * df.result$p.value
  return(df.result)
}

cohorts <- unique(as.character(lmx$Site))

fligner.all <- run.fligner(c(cohorts))

flig.res <- llply(as.list(cohorts), function(i){
  df <- run.fligner(i)
  df
})
names(flig.res) <- cohorts

flig.res2 <- c(flig.res, list(fligner.all))
names(flig.res2) <- c(cohorts, "all")
groups <- c(cohorts, "all")

#### Make list of cytokines to exclude ####

cyto.ex <- llply(flig.res2, function(i){
  
  df.remove <- filter(i, p.value > 0.05)
  cyto.remove <- as.character(df.remove$cyto.stim)
  cyto.remove
  
})

#### Filter baseline data ####
# we will only select unstim variables where over 70% of subjects had a measured response AND the difference between TRUE and FALSE MFI is  > 10

lmx.us <- filter(lmx, Stim == "Unstim")
lmx.us$cutoff <- substr(lmx.us$rawConcentration, start = 1, stop = 1)
lmx.us$cutoff <- ifelse(lmx.us$cutoff == "<" | lmx.us$Bead.Name == "IL-23" & lmx.us$final.concentration == 0.04, "FALSE", "TRUE")

cohort.list <- c(cohorts, list(cohorts))

us.res <- llply(as.list(cohort.list), function(i){
  
  df <- filter(lmx.us, Site %in% c(i))
  cytokines <- unique(as.character(df$cyto.stim))
  
  df.result <- ldply(as.list(cytokines), function(i){
    
    df.test <- filter(df, cyto.stim == i)
    freq.true <- length(which(df.test$cutoff == "TRUE")) / length(df.test$cutoff)
    freq.mfi10 <- length(which(df.test$MFI < 10)) / length(df.test$MFI)
    
    cyto.stim <- i
    res <- data.frame(cyto.stim, freq.true, freq.mfi10)
  })
  
  df.result
})
names(us.res) <- groups

# selet unstim data to remove
us.remove <- llply(us.res, function(i){
  
  res <- i
  res <- filter(res, freq.true <= 0.7 & freq.mfi10 > 0.3)
  rm <- as.character(res$cyto.stim)
  rm
  
})

# Make a list of all data to remove

list.remove <- llply(as.list(1:length(cyto.ex)), function(i){
  
  us.rm <- us.remove[[i]]
  cyto.rm <- cyto.ex[[i]]
  rm <- c(us.rm, cyto.rm)
  rm
  
})
names(list.remove) <- groups

#### Create luminex dataframes with removed cytokines ####
cohorts <- unique(as.character(lmx$Site))
cohort.list <- as.list(cohorts)
cohort.list <- c(cohort.list, list(c(cohorts)))

lmx.trim <- llply(as.list(1:length(cohort.list)), function(i){
  
  cohort <- cohort.list[[i]]
  remove <- list.remove[[i]]
  
  lmx.site <- filter(lmx, Site %in% c(cohort))
  lmx.trim <- lmx.site[-which(lmx.site$cyto.stim %in% c(remove)),]
  lmx.trim
  
})

names(lmx.trim) <- groups

#### Write out filtered luminex data ####
# filtered datasets
#saveRDS(lmx.trim, "ms_global_cohort_analysis/Rdata/R_export/ms_gc_lmx_filtered_sPLS.rds")

# END ####