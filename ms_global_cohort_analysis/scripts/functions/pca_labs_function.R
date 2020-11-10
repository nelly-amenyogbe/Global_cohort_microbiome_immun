# Nelly Amenyogbe

# This function returns the percent variance explained in a prcomp PCA object, for the first 3 components.
# Input: output of prcomp()
# Output: list. each item contains a character label for percent variation explainded by each PC.

get.pca.dat <- function(pca.object){
  

  pc.dat <- summary(pca.object)
  pc.dat <- pc.dat$importance
  
  
  pc1.var <- paste("PC1:", 100*pc.dat[2,1], "% expl. var")
  pc2.var <- paste("PC2:", 100*pc.dat[2,2], "% expl. var")
  pc3.var <- paste("PC3:", 100*pc.dat[2,3], "% expl. var")
  
  to.return <- list(pc1.var, pc2.var, pc3.var)
  names(to.return) <- c("pc1", "pc2", "pc3")
  
  return(to.return)
  
}