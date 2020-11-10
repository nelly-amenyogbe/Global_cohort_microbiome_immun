####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Microbiome: Alpha Diversity
####################

# In this script, we will determine the observed richness and Shannon Diversity of all samples, and apply the Kruskal-Wallis test to determine if either differ across cohorts. Further, we perform linear regression across cohorts to determine if any host factors affect richness and diversity measurements.

# load packages
library(phyloseq)
library(ggplot2)
library(plyr)
library(dplyr)
library(splines)
# helper function
source("ms_global_cohort_analysis/scripts/functions/function_run_lm_univariate.R")

# load data
ps <- readRDS("ms_global_cohort_analysis/Rdata/human_raw_data/gc_physeq.rds")
meta <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_metadata.csv")
bfeed <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_breastfeeding.csv") # breastfeeding duration for Canadian and Ecuadorean children

# add metadata to physeq
ps.meta <- data.frame(sample_data(ps))

meta <- filter(meta, Microbiome == "Y")
rownames(meta) <- meta$Microbiome_ID
colnames(meta)

# create age in months
meta$age <- ifelse(is.na(meta$age.stool.months), meta$age.months, meta$age.stool.months)

# add breastfeeding data
meta <- join(meta, bfeed[,c("SUBJECT_ID", "ever.breastfed", "breastnow", "agesuspendbrfeed.mo")], by = "SUBJECT_ID")

# calculate time since breastfeeding
meta$time.since.bfeed <- ifelse(is.na(meta$breastnow) | meta$breastnow == "N", meta$age - meta$agesuspendbrfeed.mo, 
                                ifelse(meta$breastnow == "Y", 0, NA))

# change negative values to 0
meta$time.since.bfeed <- ifelse(meta$time.since.bfeed < 0, 0, meta$time.since.bfeed)

meta$bfeed.duration <- ifelse(is.na(meta$breastnow) | meta$breastnow == "N", meta$agesuspendbrfeed.mo,
                              ifelse(meta$breastnow == "Y", meta$age, NA))

# select metadata for all cohorts
meta.select <- meta[,c("Microbiome_ID", "ETHNICITY", "COUNTRY", "SEX", "GEST_AGE", "BIRTH_WT", "DELIVERY", "MOM_AGE", "WAZ", "WLZ", "HAZ", "age", "time.since.bfeed", "bfeed.duration")]

# order to be same as metadata
meta.select <- meta.select[match(ps.meta$X.SampleID, meta.select$Microbiome_ID),]
length(which(meta.select$Microbiome_ID == rownames(ps.meta)))
rownames(meta.select) <- meta.select$Microbiome_ID

sampledat <- sample_data(meta.select)
ps <- merge_phyloseq(ps, sampledat)

# get sample depth
ps.meta <- data.frame(sample_data(ps))
depth <- colSums(otu_table(ps))
ps.meta$depth <- depth

# Plot depth ####
gc.cols <- c("BLG" = "#fec44f", "CAD" = "#ef3b2c", "ECD" = "#225ea8", "SAF" = "#41ab5d")
gc.shapes <- c("BLG" = 19, "CAD" = 17, "ECD" = 18, "SAF" = 15)
gc.size <- c(4,4,5,4)

p.depth <- ggplot(ps.meta, aes(x = Cohort, y = depth, color = Cohort)) +
  geom_boxplot(outlier.size = NA) +
  geom_point(position = position_jitter(width = 0.2)) +
  #scale_shape_manual(values = gc.shapes) +
  scale_colour_manual(values = gc.cols) + 
  scale_size_manual(values = gc.size) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(x = "", y = "Number of reads")

p.depth

# rarefy ####
# remove subjet with less than 20000 reads
gc.rar <- prune_samples(sample_sums(ps) > 20000, ps)
gc.rar <- rarefy_even_depth(gc.rar, rngseed = 711)

rar.sub <- colSums(otu_table(gc.rar))
unique(rar.sub) # rarefied to 30,246 reads per sample

# Compute Richness
gc.rich <- plot_richness(gc.rar, measures = "Shannon")
rich.dat <- gc.rich$data

# Plot Richness by cohort ####

p.shannon <- ggplot(rich.dat, aes(x = Cohort, y = value, color = Cohort)) +
  geom_boxplot(outlier.size = NA) +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_colour_manual(values = gc.cols) + 
  scale_size_manual(values = gc.size) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(x = "", y = "Shannon Index")

p.shannon

#ggsave("ms_global_cohort_analysis/figures/microbiome/boxplot_shannon_diversity.pdf", device = "pdf", dpi = 300, width = 2.7, height = 4.0)

# GC: Observed Richness ####

# Compute Richness
gc.rich.obs <- plot_richness(gc.rar, measures = "Observed")
rich.dat.obs <- gc.rich.obs$data

# Plot Richness by cohort ####

p.observed <- ggplot(rich.dat.obs, aes(x = Cohort, y = value, color = Cohort)) +
  geom_boxplot(outlier.size = NA) +
  geom_point(position = position_jitter(width = 0.2)) +
  #scale_shape_manual(values = gc.shapes) +
  scale_colour_manual(values = gc.cols) + 
  scale_size_manual(values = gc.size) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(x = "", y = "Observed Richness")

p.observed

#ggsave("ms_global_cohort_analysis/figures/microbiome/boxplot_observed_richness.pdf", device = "pdf", dpi = 300, width = 2.7, height = 4.0)

# Shannon: LM Statistics ####
all.lm <- lm(value ~ Cohort, data = rich.dat)
summary(all.lm)

# Adj. R2 = 0.01352, p = 0.2352

# Observed richness:  effect of cohort ####
all.lm.obs <- lm(value ~ Cohort, data = rich.dat.obs)
summary(all.lm.obs)

# Adj. R2 = 0.1544, p = 0.0002979, ECD p = 0.0428, NS SAF/ECD

# Shannon Diversity: demographics ####

# Create variables to test for each cohort
cad.vars <- c("SEX", "DELIVERY", "age", "GEST_AGE", "MOM_AGE", "WAZ", "HAZ", "WLZ", "time.since.bfeed", "bfeed.duration")
blg.vars <- c("age", "GEST_AGE", "MOM_AGE", "WAZ", "HAZ", "WLZ")
saf.vars <- c("SEX", "age", "GEST_AGE", "MOM_AGE", "WAZ", "HAZ", "WLZ")
ecd.vars <- c("SEX", "DELIVERY", "age", "GEST_AGE", "MOM_AGE", "WAZ", "HAZ", "WLZ", "time.since.bfeed", "bfeed.duration")

# Shan LM: CAD ####
cad.lm <- run.lm.univariate(filter(rich.dat, Cohort == "CAD"), cad.vars)

# get significant results
cad.sig <- filter(cad.lm, variable != "(Intercept)", p.value < 0.05)
cad.sig

# DELIVERY R2 = 0.1537992; p = 0.02910, MOM_AGE R2 = 0.2575; p = 0.00356


# CAD: Model all significant terms
cad.lm.adj <- lm(value ~ DELIVERY + MOM_AGE, data = filter(rich.dat, Cohort == "CAD"))
summary(cad.lm.adj)

# MomAge p = 0.0141, DeliveryMode p = 0.1192, Model R2 = 0.271

# Conclude that the effect of DeliveryMode is not significant if adjusting for maternal age.

# Shan lm: BLG ####
blg.lm <- run.lm.univariate(filter(rich.dat, Cohort == "BLG"), blg.vars)

# get significant results
blg.sig <- filter(blg.lm, variable != "(Intercept)", p.value < 0.05)
blg.sig # No significant interactions

# Shan lm: SAF ####
saf.lm <- run.lm.univariate(filter(rich.dat, Cohort == "SAF"), saf.vars)

# get significant results
saf.sig <- filter(saf.lm, variable != "(Intercept)", p.value < 0.05) # No significant interactions
saf.sig # no significant interactions

# Shan lm: ECD ####
# Add breastfeeding duration to ECD data ####
ecd.dat <- filter(rich.dat, Cohort == "ECD")

# run LM with first 5 vars
ecd.lm <- run.lm.univariate(ecd.dat, ecd.vars[1:5])

# get significant results
ecd.sig <- filter(ecd.lm, variable != "(Intercept)", p.value < 0.05)
ecd.sig

# MOM_AGE R2 = 0.2219; p = 0.001639, GEST_AGE R2 = 0.1161, p = 0.04852; 

# run LM with first 6:10 vars
ecd.lm2 <- run.lm.univariate(ecd.dat, ecd.vars[6:10])

# get significant results
ecd.sig2 <- filter(ecd.lm2, variable != "(Intercept)", p.value < 0.05)
ecd.sig2

# None

# ECD model all significant terms
ecd.lm.adj <- lm(value ~ MOM_AGE + GEST_AGE, data = filter(rich.dat, Cohort == "ECD"))
summary(ecd.lm.adj)

# MOM_AGE p = 0.00346, GEST_AGE p = 0.08317, Model R2 = 0.2894

# Conclude that the relationship between Maternal Age and richness is independent of Gestational Age, but the effect of GA on diversity is not significant when adjusting for maternal age

# Plot significant findings ####
dm.cols <- c("Caesarean" = "#225ea8", "Vaginal" = "#fc9272")

ce.dat <- filter(rich.dat, Cohort %in% c("CAD", "ECD"))
ce.dat$Cohort <- factor(ce.dat$Cohort, levels = c("ECD", "CAD"))

# boxplot

ggplot(ce.dat, aes(x = Cohort, y = value, color = DELIVERY)) +
  geom_boxplot(outlier.size = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.6)) +
  theme_classic() +
  scale_colour_manual("Delivery Mode", values = dm.cols) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right",
        strip.text.x = element_text(size = 12, face = "bold")) +
  labs(x = "Maternal Age", y = "Shannon Index")

#ggsave("ms_global_cohort_analysis/figures/microbiome/boxplot_momage_diversity.pdf", device = "pdf", dpi = 300, width = 4, height = 4)


# plot linegraph
ggplot(ce.dat, aes(x = MOM_AGE, y = value, color = DELIVERY)) + 
  geom_point() + 
  geom_smooth(color = "black", method = "lm") + 
  facet_grid(~Cohort) +
  theme_classic() +
  scale_colour_manual("Delivery Mode", values = dm.cols) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        strip.text.x = element_text(size = 12, face = "bold")) +
  labs(x = "Maternal Age", y = "Shannon Index")

#ggsave("ms_global_cohort_analysis/figures/microbiome/cad_ecd_momage_diversity.pdf", device = "pdf", dpi = 300, height = 3.3, width = 6.5)

# plot together
ggplot(ce.dat, aes(x = MOM_AGE, y = value, color = Cohort)) + 
  geom_point(size = 2) + 
  geom_smooth(color = "black") + 
  facet_grid(~DELIVERY) +
  theme_classic() +
  scale_colour_manual(values = gc.cols) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right",
        strip.text.x = element_text(size = 12, face = "bold")) +
  labs(x = "Maternal Age", y = "Shannon Index")

# U-shape stats ####
with(ce.dat, scatter.smooth(value ~ MOM_AGE))

# test U-shape
# source:  http://rcompanion.org/rcompanion/e_03.html

# all samples
fit.ma <- lm(value ~ bs(MOM_AGE,
                        knots = 5,
                        degree = 2),
             data = ce.dat)

summary(fit.ma) # Adjusted R2 = 0.09. p = 0.01374

# VD only
vd.dat <- filter(ce.dat, DELIVERY == "Vaginal")
fit.vd <- lm(value ~ bs(MOM_AGE,
                        knots = 5,
                        degree = 2),
             data = vd.dat)

summary(fit.vd) #p=0.0415 on 2 knots.  Adusted R2 = 0.0947, p = 0.03812

# Conclude:  There is a significant quadratic relationship between maternal age and shannon diversity among the Canadian and Ecuadorean infants.  This can be detected among vaginally-delivered infants alone.

# plot all dat
plot(value ~ MOM_AGE,
     data = ce.dat,
     pch = 16,
     cex = 0.8,
     xlab = "Maternal Age (years)",
     ylab = "Shannon Diversity")

i = seq(min(ce.dat$MOM_AGE), max(ce.dat$MOM_AGE), len = 100)
predy = predict(fit.ma, data.frame(MOM_AGE = i))
lines(i, predy,
      lty=1, lwd=2, col="blue")

# export as PDF

# check residuals
hist(residuals(fit.ma),
     col="darkgray") # residuals have a normal distribution.

plot(fitted(fit.ma),
     residuals(fit.ma)) # residuals appear homoskedastic and unbiased.

# Plot ECD GA ####
ggplot(filter(rich.dat, Cohort == "ECD"), aes(x = GEST_AGE, y = value)) + 
  geom_point(size = 2, position = position_jitter(width = 0.1)) + 
  geom_smooth(color = "black", method = "lm") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right",
        strip.text.x = element_text(size = 12, face = "bold")) +
  labs(x = "Gestational age", y = "Shannon Index") +
  ggtitle("Gestational Age and diversity in ECD")

# END ####
