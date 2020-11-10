#####################
# Nelly Amenyogbe
# sPLS-DA:  useful functions
####################


# mc.get.features ####
# Input either a PLS-DA or sPLS-DA object, and comp as "comp.1", "comp.2" etc...
# Returns a data frame of all features in object, with loadings and absolute value of loadings.  Features are listen in order of highest to lowest absolute loading.

mc.get.features <- function(mixomics.obj, comp){
  
  # get loadings from mixomics object
  loadings <- data.frame(mixomics.obj$loadings$X)
  loadings$feature <- rownames(loadings)
  
  # modify data frame to only include comp of interest
  loadings <- loadings[,c("feature", comp)]
  loadings$raw.loadings <- loadings[,2]
  
  # add absolute loadings and order features by strength
  loadings$abs.loadings <- abs(loadings$raw.loadings)
  loadings <- loadings[order(loadings$abs.loadings, decreasing = TRUE),]
  
  # remove redundant column
  loadings <- loadings[,-2]
  return(loadings)
  
}

# get error data ####
get.error.dat <- function(perf.obj){
  
  # get class error rate
  error <- perf.obj
  er.class <- data.frame(error$error.rate.class$max.dist)
  er.class$group <- rownames(er.class)
  ecm <- melt(er.class, id.vars = "group")
  
  # get overall error rate
  er.all <- data.frame(error$error.rate.all$overall$max.dist)
  er.all$variable <- rownames(er.all)
  er.all$variable <- gsub(" ", ".", er.all$variable)
  er.all$group <- "Overall"
  colnames(er.all)[1] <- "value"
  er.all <- er.all[,c(colnames(ecm))]
  
  # combine class and overall
  er.df <- rbind(ecm, er.all)
  
  # arrage group factor levels to keep overall last
  groups <- unique(as.character(ecm$group))
  groups <- c(groups, "Overall")
  er.df$group <- factor(er.df$group, levels = c(groups))
  
  return(er.df)
}
