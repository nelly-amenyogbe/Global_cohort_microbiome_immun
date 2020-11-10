###############
# Nelly Amenyogbe
# 08-March-2018
# Function:  Run univariate Lms
# Description:  This function will run a regression on every specified variable against a numeric variable that must be called "value"
################

# Input: df: wide data frame (n row = n samples), with numerical value for regression to be saved as "value".
# Output: results table with the estimate, R2 value, and p-value computed from the linear regression

run.lm.univariate <- function(df, vars){

  lm.res <- ldply(as.list(1:length(vars)), function(i){
    
    col.num <- which(colnames(df) == vars[i])
    colnames(df)[col.num] <- "ind.var"
    
    fit <- lm(value ~ ind.var, data = df)
    
    x <- summary(fit)
    r2 <- x$r.squared
    y <- as.data.frame(x$coefficients)
    y[,1] <- round(y[,1], digits = 3)
    y[,2] <- round(y[,2], digits = 3)
    y$confint97.5 <- y[,1] + y[,2]
    y$confint2.5 <- y[,1] - y[,2]
    y$R2 <- r2
    y$variable <- rownames(y)
    y$ind.var <- vars[i]
    colnames(y)[4] <- "p.value"
    y
    
    
  })
  
  return(lm.res)
}
