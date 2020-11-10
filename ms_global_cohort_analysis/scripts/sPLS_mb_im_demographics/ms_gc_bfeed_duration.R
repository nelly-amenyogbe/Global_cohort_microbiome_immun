#####################
# Nelly Amenyogbe
# Global Cohort Manuscript Analysis
# Breast feeding duration: associations to select OTUs among Canadian and Ecuadorean children
#####################

# In this script, we visualize the breast feeding duration for Canadian and Ecuadorean children.  We then plot associations identified between time since breast feeding among Canadian children, and investigate whether these trends are evident among Ecuadorean children.


# load packages
library(plyr)
library(dplyr)
library(ggplot2)

# load data
meta.all <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_metadata.csv")
bfeed <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_breastfeeding.csv")
spls.dat <- readRDS("ms_global_cohort_analysis/Rdata/R_export/spls_data.rds")

# colours that will be used for plotting
cohort.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_mb_cols.csv")

# prepare data for plotting ####
# spls data
spls.data <- spls.dat$ALL

# add breastfeeding data to metadata ####
# create age in months
meta.all$age <- ifelse(is.na(meta.all$age.stool.months), meta.all$age.months, meta.all$age.stool.months)

meta.all <- join(meta.all, bfeed[,c("SUBJECT_ID", "ever.breastfed", "breastnow", "agesuspendbrfeed.mo")], by = "SUBJECT_ID")

# calculate time since breastfeeding
meta.all$time.since.bfeed <- ifelse(is.na(meta.all$breastnow) | meta.all$breastnow == "N", meta.all$age - meta.all$agesuspendbrfeed.mo, 
                                    ifelse(meta.all$breastnow == "Y", 0, NA))

# change negative values to 0
meta.all$time.since.bfeed <- ifelse(meta.all$time.since.bfeed < 0, 0, meta.all$time.since.bfeed)

meta.all$bfeed.duration <- ifelse(is.na(meta.all$breastnow) | meta.all$breastnow == "N", meta.all$agesuspendbrfeed.mo,
                                  ifelse(meta.all$breastnow == "Y", meta.all$age, NA))


# prepare metadata ####
spls.subjects <- c(rownames(spls.dat$CAD$Xotu), rownames(spls.dat$ECD$Xotu))

meta <- filter(meta.all, SUBJECT_ID %in% spls.subjects)

rownames(meta) <- meta$SUBJECT_ID
meta <- meta[order(meta$SUBJECT_ID),]

# prepare mb and immun data ####
X <- spls.dat$ALL$Xotu
Xs <- X[which(rownames(X) %in% spls.subjects),]
Xs <- Xs[order(rownames(Xs)),]

# plot Bfeeding ####

# here we plot the distribution of breast feeding duration among CAD and ECD children

bar.dat <- meta[,c("SUBJECT_ID", "Cohort", "time.since.bfeed", "bfeed.duration")]

bar.dat <- bar.dat[order(bar.dat$bfeed.duration, decreasing = TRUE),]

# generate melted data
bar.m <- melt(bar.dat, id.vars = c("SUBJECT_ID", "Cohort"), value.vars = c("time.since.bfeed", "bfeed.duration"))

# set graphing variable names
bar.m$variable <- ifelse(bar.m$variable == "time.since.bfeed", "weaned", ifelse(bar.m$variable == "bfeed.duration", "breastfed", bar.m$variable))

colnames(bar.m)[c(3:4)] <- c("Bfeed.duration", "age.months")

# set factor levels
bar.m$SUBJECT_ID <- factor(bar.m$SUBJECT_ID, levels = c(as.character(bar.dat$SUBJECT_ID)))

bar.m$Bfeed.duration <- factor(bar.m$Bfeed.duration, levels = c("weaned", "breastfed"))

ggplot(bar.m, aes(x = SUBJECT_ID, y = age.months, fill = Bfeed.duration)) + geom_bar(stat = "identity", color = "black") +
  facet_grid(~Cohort, scales = "free", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_fill_manual(values = c("weaned" = "#31a354", "breastfed" = "#ffeda0")) +
  labs(x = "", y = "Age (months)")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad_ecd/bfeed_duration_spls_subjects.pdf", device = "pdf", dpi = 300, width = 6.5, height = 3.2)

# plot b.feed OTUs ####

# prepare colors
gc.cols <- as.character(cohort.cols$Color)
names(gc.cols) <- cohort.cols$Site

# merge OTU data with metadata
otu.dat <- data.frame(Xs)
otu.dat$SUBJECT_ID <- rownames(otu.dat)
otu.dat <- join(otu.dat, meta[,c("SUBJECT_ID", "Cohort", "time.since.bfeed", "bfeed.duration", "breastnow", "time.since.bfeed", "SEX", "DELIVERY")], by = "SUBJECT_ID")

# OTU_48 F_Lachnospiraceae ####

# CAD
fit <- lm(OTU_48.F_Lachnospiraceae ~ time.since.bfeed, data = filter(otu.dat, Cohort == "CAD"))

summary(fit) # p = 0.017, Adj. R2 = 0.249

# ECD
fit <- lm(OTU_48.F_Lachnospiraceae ~ time.since.bfeed, data = filter(otu.dat, Cohort == "ECD"))
summary(fit) # p = 0.3, Adj. R2 = 0.31

ggplot(otu.dat, aes(x = time.since.bfeed, y = OTU_48.F_Lachnospiraceae)) + 
  geom_smooth(method = "lm", color = "black") +
  geom_point(shape = 21, size = 3, position = position_jitter(width = 0.2), aes(fill = Cohort)) +
  scale_fill_manual(values = gc.cols) +
  facet_grid(~Cohort) +
  ggtitle("Lachnospiraceae") +
  theme(legend.position = "bottom")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad_ecd/otu48lachno_bf.pdf", device = "pdf", dpi = 300, width = 5, height = 3.5)

# OTU_30 G_Roseburia

# CAD
fit <- lm(OTU_61.G_Roseburia ~ time.since.bfeed, data = filter(otu.dat, Cohort == "CAD"))
summary(fit) # p = 0.014, Adj. R2 = 0.264

# ECD
fit <- lm(OTU_61.G_Roseburia ~ time.since.bfeed, data = filter(otu.dat, Cohort == "ECD"))
summary(fit) # p = 0.76, Adj. R2 = 0.76

ggplot(otu.dat, aes(x = time.since.bfeed, y = OTU_61.G_Roseburia)) + 
  geom_smooth(method = "lm", color = "black") +
  geom_point(shape = 21, size = 3, position = position_jitter(width = 0.2), aes(fill = Cohort)) +
  scale_fill_manual(values = gc.cols) +
  facet_grid(~Cohort) +
  ggtitle("Roseburia") +
  theme(legend.position = "bottom")

#ggsave("ms_global_cohort_analysis/figures/sPLS_mb_im_demographics/cad_ecd/otu61roseburia_bf.pdf", device = "pdf", dpi = 300, width = 5, height = 3.5)

# END ####