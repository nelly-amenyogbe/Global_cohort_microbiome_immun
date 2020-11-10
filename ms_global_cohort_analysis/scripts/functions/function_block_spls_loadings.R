#########
# Nelly Amenyogbe
# block sPLS: get features
#########

# input loadings data frame extracted from block.spls object, returns data frame with loadings for the component requested ordered by absolute loading

# block.spls.feature.summary
block.spls.get.features <- function(loadings.data, comp){
  
  loadings <- loadings.data
  
  loadings$feature <- rownames(loadings)
  loadings.select <- loadings[,which(colnames(loadings) %in% c("feature", comp))]
  loadings.select$abs.loading <- abs(loadings.select[,1])
  loadings.select <- loadings.select[order(loadings.select$abs.loading, decreasing = TRUE),]
  
  return(loadings.select)

}


