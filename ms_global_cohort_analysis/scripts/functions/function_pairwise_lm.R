###################
# Nelly Amenyogbe
# 01 Feb 2018
# Lm:  generic function for multiple pairwise lm tests between two data matrices
##################

# these functions were written to perform a linear regression between fetures in paired data sets, by using the vales as they are in the input data, or scaled.

# This function requires a data frame of all OTUs and cytokines that need to be regressed against each other. 

# lm.df.dataprep ####
# Input: X and Y data matrices as used for mixOmics::spls()
# Output: molten data frame, n "X" features x n "Y" features long

lm.df.dataprep <- function(x.mat, y.mat){
  
  y.mat <- as.data.frame(y.mat) # the response variable
  x.mat <- as.data.frame(x.mat) # the predictor variable
  
  # prepare cytokine df
  y.mat$sample.id <- rownames(y.mat)
  y.melt <- melt(y.mat, id.vars = "sample.id", variable.name = "Y.var", value.name = "Y.value")
  
  # prepare otu df
  x.mat$sample.id <- rownames(x.mat)
  x.melt <- melt(x.mat, id.vars = "sample.id", variable.name = "X.var", value.name = "X.value")
  
  # join
  df <- join(y.melt, x.melt, by = "sample.id")
  df
  
}

# x.y.lm ####
# Lm between each cyto-otu pair
# Input: data frame output from lm.df.dataprep, or a subset thereof
# Output: result table of linear regression of all X features by all Y features. Note that this is equivalent to the pearson correlation coefficient, in terms of estimate and p-value

x.y.lm <- function(df){
  
  y.list <- unique(as.character(df$Y.var))
  x.list <- unique(as.character(df$X.var))
  
  res <- ldply(as.list(y.list), function(i){
    
    df <- filter(df, Y.var == i)
    df.result <- ldply(as.list(as.character(x.list)), function(i){
      
      df <- filter(df, X.var == i)
      fit <- lm(Y.value ~ X.value, data = df)
      x <- summary(fit)
      r2 <- x$r.squared
      y <- as.data.frame(x$coefficients)
      y[,1] <- round(y[,1], digits = 3)
      y[,2] <- round(y[,2], digits = 3)
      y$confint97.5 <- y[,1] + y[,2]
      y$confint2.5 <- y[,1] - y[,2]
      y$R2 <- r2
      y$variable <- rownames(y)
      y$X.var <- i
      colnames(y)[4] <- "p.value"
      y
      
    })
    
    df.result$Y.var <- i
    colnames(df.result)[4] <- "p.value"
    df.result <- filter(df.result, variable != "(Intercept)")
    df.result$p.adj <- p.adjust(df.result$p.value, method = "BH", n = length(x.list))
    df.result
  })
  
}

# x.y.lm.scaled ####
# Scaled Lm between each cyto-otu pair 
# Lm between each cyto-otu pair, computed by scaling the values prior to regression
# Input: data frame output from lm.df.dataprep, or a subset thereof
# Output: result table of linear regression of all X features by all Y features. Note that this is equivalent to the pearson correlation coefficient, in terms of estimate and p-value

x.y.lm.scaled <- function(df){
  
  y.list <- unique(as.character(df$Y.var))
  x.list <- unique(as.character(df$X.var))
  
  res <- ldply(as.list(y.list), function(i){
    
    df <- filter(df, Y.var == i)
    df.result <- ldply(as.list(as.character(x.list)), function(i){
      
      df <- filter(df, X.var == i)
      fit <- lm(scale(Y.value) ~ scale(X.value), data = df)
      x <- summary(fit)
      r2 <- x$r.squared
      y <- as.data.frame(x$coefficients)
      y[,1] <- round(y[,1], digits = 3)
      y[,2] <- round(y[,2], digits = 3)
      y$confint97.5 <- y[,1] + y[,2]
      y$confint2.5 <- y[,1] - y[,2]
      y$R2 <- r2
      y$variable <- rownames(y)
      y$X.var <- i
      colnames(y)[4] <- "p.value"
      y
      
    })
    
    df.result$Y.var <- i
    colnames(df.result)[4] <- "p.value"
    df.result <- filter(df.result, variable != "(Intercept)")
    df.result$p.adj <- p.adjust(df.result$p.value, method = "BH", n = length(x.list))
    df.result
  })
  
}

# spls.get.features ####
# Input: spls.object = output of mixOmics::spls(). comp.no: component of choice (e.g. "comp.1") as a character
# Output: list; data frames for X features and Y features containing loading, absolute loadings, and model component

spls.get.features <- function(spls.object, comp.no){
  
  # get features with loadings
  x.loadings <- data.frame(spls.object$loadings$X)
  x.loadings$x.feature <- rownames(x.loadings)
  
  y.loadings <- data.frame(spls.object$loadings$Y)
  y.loadings$y.feature <- rownames(y.loadings)
  
  x.m <- melt(x.loadings, id.vars = "x.feature")
  colnames(x.m) <- c("x.feature", "comp", "loading")
  
  y.m <- melt(y.loadings, id.vars = "y.feature")
  colnames(y.m) <- c("y.feature", "comp", "loading")
  
  # filter loadings according to preferences
  x.filter <- filter(x.m, comp == comp.no, loading != 0)
  y.filter <- filter(y.m, comp == comp.no, loading != 0)
  
  to.return <- list(x.filter, y.filter)
  names(to.return) <- c("x.ld", "y.ld")
  return(to.return)
  
  
  
}
