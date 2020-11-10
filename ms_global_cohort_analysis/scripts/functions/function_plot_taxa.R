####################
# Nelly Amenyogbe
# Functions for plotting abundant species from phyloseq objects
####################


# get_taxbar_data ####
# ps  = phyloseq object with OTUs already transformed to relative abundance (unless you want to aggregate counts)
# num.taxa = the number for top most abundat taxa you want to plot (e.g. 10)
# rank = the taxonomic rank you want to prepare data for (e.g. "Family")

get_taxbar_data <- function(ps, num.taxa, rank){

  p <- plot_bar(ps, fill = rank)
  dat <- p$data
  
  # change the name of the tax.rank of interest to "Taxa
  col.num <- which(colnames(dat) ==  rank)
  colnames(dat)[col.num] <- "Taxa"
  
  # Determine top most abundant taxa in dataset
  tax.ag <- aggregate(Abundance ~ Taxa, data = dat, sum)
  tax.ag <- tax.ag[order(tax.ag$Abundance, decreasing = TRUE),]
  
  # Retain only those that were classified to that level
  tax.ag$max.rank <- substr(tax.ag$Taxa, start = 1, stop = 1)
  
  tax.rank <- substr(rank, start = 1, stop = 1)
  
  tax.ag <- filter(tax.ag, max.rank == tax.rank)
  
  # get top taxa
  top.taxa <- tax.ag$Taxa[1:num.taxa]
  
  # Change the taxon column to contain the top taxa and list the rest as "others"
  
  dat$max.rank <- substr(dat$Taxa, start = 1, stop = 1)
  
  dat$Taxa <- ifelse(dat$Taxa %in% top.taxa, as.character(dat$Taxa),
                     ifelse(dat$max.rank == tax.rank, "Other", "Unknown"))
  
  # Format for plotting
  
  dat$Taxa <- factor(dat$Taxa, levels = c(as.character(top.taxa), "Other", "Unknown"))
  
  return(dat)
}

###################



