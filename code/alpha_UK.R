################################# 
## R script                    ##
## Project: TwinsRA            ##
## Alpha diversity             ##
## Data: Shotgun metagenomics  ##
## Author: JM + KB             ##
## Date: 3/7/21                ##
## Last Updated: 5/9/24        ##
#################################

### Load and save current R script ###
# Load R scripts
# load(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs20/alpha_data.RData")
# load(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs20/alpha_script.RData")

# Save R script
# Do this step prior to closing R
# save.image(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs20/alpha_script.RData")

############################################################################
############################################################################
############################################################################

### Load libraries ###
library(reshape2)
library(phyloseq)
library(vegan)
library(ade4)
library(PMCMRplus)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)

############################################################################
############################################################################
############################################################################

### Statistics functions ###

# All plot statistics: mean, std deviation, median, min value, max value, 10%ile, 25%ile, 75%ile, 90%ile
stats.all = function(x) {
  mean <- mean(x)
  stddev <- sd(x)
  median <- median(x)
  val_min <- min(x)
  val_max <- max(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(mean = mean, sd = stddev, median = median, 
           val_min = val_min, val_max = val_max, 
           per10 = per10, per25 = per25, per75 = per75,  per90 = per90))
}

# Boxplot statistics: median, 25%ile, 75%ile
stats.boxplot <- function(x) {
  m <- median(x)
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  return(c(y = m, ymin = per25, ymax = per75))
}

# Whiskers statistics: median, 10th percentile, 90th percentile
stats.whiskers = function(x) {
  m <- median(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(y = m, ymin = per10, ymax = per90))
}

# Outliers
min.outlier <- function(x) {
  subset(x, quantile(x, prob = c(0.10)) > x)
}

max.outlier <- function(x) {
  subset(x, quantile(x, prob = c(0.90)) < x)
}

############################################################################
############################################################################
############################################################################

# set working dir
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/"

### Alpha Diversity Boxplots ###
df_alpha = read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/metadata_shannon_UK.tsv', 
                      sep = '\t', header = TRUE, row.names = 1, check.names = FALSE,
                      na.strings = "NA")

# convert twin pair to consistent ID
#df_alpha$TwinPair = c('sib_01','sib_01','sib_02','sib_02','sib_03','sib_03',
#                      'sib_04','sib_04','sib_05','sib_05','sib_06','sib_06',
#                      'sib_07','sib_07','sib_08','sib_08')


# create tables for storing wilcoxon and ttest results
stats.table.all <- matrix(data = NA, nrow = 1, ncol = 3)
colnames(stats.table.all) <- c("alpha div", "wilcoxon", "ttest")

# calculate adiv
stats.table.all[1,1] <- colnames(df_alpha)[4]
#stats.table.all[1,2] <- wilcox.test(df_alpha[,4] ~ Diagnosis, data = df_alpha, paired = TRUE)$p.value
stats.table.all[1,2] <- wilcox.test(x=df_alpha[df_alpha$Diagnosis == 'RA',]$shannon_entropy, # convention for direction is x - y
                                    y=df_alpha[df_alpha$Diagnosis == 'Unaffected',]$shannon_entropy,
                                    paired=TRUE)$p.value

#stats.table.all[1,3] <- t.test(df_alpha[,4] ~ Diagnosis, data = df_alpha, paired = TRUE)$p.value
stats.table.all[1,2] <- t.test(x=df_alpha[df_alpha$Diagnosis == 'RA',]$shannon_entropy, # convention for direction is x - y
                                    y=df_alpha[df_alpha$Diagnosis == 'Unaffected',]$shannon_entropy,
                                    paired=TRUE)$p.value

# save
ft.all = paste(dir, "alpha_stats.csv", sep = "")
write.csv(file = ft.all, stats.table.all)


############################################################################
############################################################################
############################################################################

### Analytes plots - all ###

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 24, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 18))+#, color = "black")) +
  theme(axis.title.y = element_text(size = 24, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 18, color = "black")) +
  theme(legend.title = element_text(size = 24, face = "bold", color = "black"))

# choose colors
col1 <- c("#929aab", "#ce2525")
col2 <- c("#f3a333", "#0074e4", "#8f8787")

# choose line types
line1 <- c("solid", "dashed", "dotted")

# variable of interest
a <- 'shannon_entropy'
# create filenames
filename_table = paste(a, "all_table_UK.csv", sep = "_")
filename_box.plot = paste(a, "all_box.plot_UK.pdf", sep = "_")  
filename_line.plot = paste(a, "all_line.plot_UK.pdf", sep = "_")  

# create categories by which to color lineplot: 1) Increase, 2) Decrease, 3) No change
# subset and spread dataset into diagnosis columns
# then calculate delta relative abundance between RA and Unaffected siblings
d.div <- df_alpha %>%
  subset(select = c("TwinPair", "Diagnosis", a)) %>%
  spread(key = "Diagnosis", value = a) %>%
  mutate(diff_sib.pair = (get("RA") - get("Unaffected"))) %>%
  mutate(
    Change.type_sib.pair = case_when(
      sign(diff_sib.pair) == 1 ~ "1_up",
      sign(diff_sib.pair) == -1 ~ "2_down",
      TRUE ~ "3_no.change"
    )
  ) 

# merge with original dataset
d.final <- gather(d.div, "Unaffected", "RA", key = "Diagnosis", value = "shannon_entropy")

# save dataset
ft = paste(dir, filename_table, sep = "")
write.csv(d.final, file = ft)

# rewrite order of factors
d.final$Diagnosis <- factor(d.final$Diagnosis, levels = c("Unaffected", "RA"))

# plot boxplot
p <- ggplot(data = d.final, aes(x = Diagnosis, y = shannon_entropy, fill = Diagnosis)) +
  stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
               color = "black", size = 0.8, width = 0.3) +
  stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
               color = "black", size = 0.5, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.5) +
  scale_x_discrete(labels = c("Unaffected", "RA")) +
  scale_fill_manual(values = col1) +      
  xlab(NULL) +
  ylab("Shannon Entropy") +
  bkg

fpb = paste(dir, filename_box.plot, sep = "")
pdf(file = fpb, height = 4.5, width = 5)
plot(p)
dev.off()

# plot lineplot
p <- ggplot(data = d.final, aes(x = Diagnosis, y = shannon_entropy, fill = Diagnosis)) +
  stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
               color = "black", size = 0.8, width = 0.3) +
  stat_summary(fun.data = stats.boxplot, geom = "crossbar",
               color = "black", size = 0.5, width = 0.5) +
  geom_point(shape = 20, size = 3, color = "black") +
  geom_line(data = subset(d.final), aes(group = TwinPair, color = Change.type_sib.pair, linetype = Change.type_sib.pair), size = 0.5) +
  geom_text_repel(data = subset(d.final, Diagnosis == "Unaffected"), aes(label = TwinPair), 
                  nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
  scale_x_discrete(labels = c("Unaffected", "RA")) +
  scale_fill_manual(values = col1) +
  scale_linetype_manual(values = line1, name = "Entropy Difference", labels = c("Increase", "Decrease", "No change")) +
  scale_color_manual(values = col2, name = "Entropy Difference", labels = c("Increase", "Decrease", "No change")) +      
  xlab(NULL) +
  ylab("Shannon Entropy") +
  bkg

fpl = paste(dir, filename_line.plot, sep = "")
pdf(file = fpl, height = 6, width = 8)
plot(p)
dev.off()
