###############
# Nelly Amenyogbe
# function:  get sPLS feature significance
# This function takes features selected by a mixomics object, and returns a summary data frame for the X and Y features that summarizes the fraction of significant connections in the network.
##############
# this requires the spls.get.features function from spls_plot_functions.R
# cor.results are the result output from x.y.lm or x.y.lm.scaled from function_pairwise_lm.R
# comp.no is written as "comp.1", "comp.2", etc...
#pv.cutoff is the p.value cutoff used to calculate the summary, e.g. 0.05.
# output gives a list for x and y features selected from the specified component, with the raw loadings and absolute value of loadings, in decreasing order of the latter, alongside the percent of p-values below the cutoff for correlations in the y-matrix. 

spls.feature.summary <- function(spls.object, cor.results, comp.no, pv.cutoff){
  
  # get selected features from the sPLS object
  features <- spls.get.features(spls.object, comp.no)
  
  # filter correlation results to only include significant features
  feat.select <- filter(cor.results, X.var %in% features$x.ld$x.feature, Y.var %in% features$y.ld$y.feature)
  
  # calculated proportion of significant values
  feat.select$pv.cut <- pv.cutoff
  
  feat.select <- ddply(feat.select, .(X.var), transform, percent.x.sig = (length(which(p.value < unique(pv.cut))) / length(p.value)) * 100)
  
  feat.select <- ddply(feat.select, .(Y.var), transform, percent.y.sig = (length(which(p.value < unique(pv.cut))) / length(p.value)) * 100)
  
  # get summaries for x and y features
  x.sig <- feat.select[,c("X.var", "percent.x.sig")]
  x.sig <- unique(x.sig)
  
  y.sig <- feat.select[,c("Y.var", "percent.y.sig")]
  y.sig <- unique(y.sig)
  
  # create summary data frames
  # X vars
  x.sum <- features$x.ld
  x.sum$abs.loading <- abs(x.sum$loading)
  colnames(x.sum)[1] <- "X.var"
  x.sum <- join(x.sum, x.sig, by = "X.var")
  x.sum <- x.sum[order(x.sum$abs.loading, decreasing = TRUE),]
  x.sum$n <- c(1:nrow(x.sum))
  
  # Y vars
  y.sum <- features$y.ld
  y.sum$abs.loading <- abs(y.sum$loading)
  colnames(y.sum)[1] <- "Y.var"
  y.sum <- join(y.sum, y.sig, by = "Y.var")
  y.sum <- y.sum[order(y.sum$abs.loading, decreasing = TRUE),]
  y.sum$n <- c(1:nrow(y.sum))
  
  to.return <- list(x.sum, y.sum)
  names(to.return) <- c("X.sum", "Y.sum")
  
  return(to.return)
  
}
