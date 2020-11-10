####################
# Nelly Amenyogbe
# GC Manuscript Analysis
# Germ-Free Mouse: Lactulose-Mannitol ratio
####################

# In this script, we compare the lactulose:mannitol ratio in urine collected from mice gavaged with Canadian or South African stool samples using a t-test, and graph the results

# load packages
library(ggplot2)
library(plyr)
library(dplyr)
library(cowplot)

# load data
dat <- read.csv("ms_global_cohort_analysis/Rdata/mouse_raw_data/gc_gfmouse_lactulose_mannitol.csv")

gc.cols <- read.csv("ms_global_cohort_analysis/Rdata/graph_aesthetics/gc_mb_cols.csv")

# Annotate groups
dat$Animal.ID
# FM IDs are CAD, HEU IDs are SAF
dat$group <- substr(dat$Animal.ID, start=1, stop=1)
dat$Cohort <- ifelse(dat$group == "F", "CAD",
                    ifelse(dat$group == "H", "SAF", dat$group))

# t-test ####
saf <- filter(dat, Cohort == "SAF")
cad <- filter(dat, Cohort == "CAD")

test <- t.test(saf$Lac.Man.ratio, cad$Lac.Man.ratio, alternative = "two.sided")
test$p.value # p = 0.0033

# get mean and standard error
colnames(dat)

# get standard error
st.err <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

# summarize data
dat.sum <- ddply(dat,
                 c("Cohort"),
                 summarise,
                 mean = mean(Lac.Man.ratio),
                 stdev = sd(Lac.Man.ratio),
                 std.err = st.err(Lac.Man.ratio),
                 max = max(Lac.Man.ratio))


# plot
grp.cols <- as.character(gc.cols$Color)
names(grp.cols) <- gc.cols$Site

ggplot(dat.sum, aes(x = Cohort, y = mean, fill = Cohort)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean, ymax = mean + std.err)) +
  scale_fill_manual(values = grp.cols) +
  geom_point(data = dat, shape = 21, size = 4, aes(x = Cohort, y = Lac.Man.ratio), position = position_jitter(width = 0.3)) +
  labs(x = "", y = "Lactulos:Mannitol ratio") +
  theme(legend.position = "none") 

#ggsave("ms_global_cohort_analysis/figures/gf_mouse/microbiome/ms_gc_gfmouse_lmratio.pdf", device = "pdf", dpi = 300, width = 2.8, height = 3.5)

# END ####