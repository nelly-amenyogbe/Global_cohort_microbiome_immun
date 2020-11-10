####################
# Nelly Amenyogbe
# 21-Sept-2020
# Global Cohort:  Plot time between stool and blood draw
####################

# In this script, we will plot the time difference between blood draw and stool sample collection, in days

# load packages
library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)

# load data
meta <- read.csv("ms_global_cohort_analysis/Rdata/human_raw_data/gc_metadata.csv") 

# Prepare data for plotting

dat <- meta[,c("SUBJECT_ID", "Cohort", "days.between.blood.stool")]
dat <- dat[order(dat$days.between.blood.stool),]
dat$SUBJECT_ID <- factor(dat$SUBJECT_ID, levels = as.character(dat$SUBJECT_ID))

# plot days between blood and stool

breaks <- seq(from = -10, to = 10, by = 2)

ggplot(dat[-which(is.na(dat$days.between.blood.stool)),], aes(x = SUBJECT_ID, y = days.between.blood.stool)) +
  geom_point(size = 2) +
  facet_grid(~Cohort, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 12)) +
  labs(x="Subject", y="Days between collection") +
  scale_y_continuous(limits = c(-10, 10), breaks = breaks) 

#ggsave("ms_global_cohort_analysis/figures/days_between_blood_stool.pdf", device = "pdf", dpi = 300, width = 12, height = 2.5)
