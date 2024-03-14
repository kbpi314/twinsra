################################# 
## R script                    ##
## Project: RA twins           ##
## Olink                       ##
## Author: KB                  ##
## Date: 3/14/24               ##
#################################

### Load and save current R script ###

# Load R scripts
load(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/RA_twins/16S/code/RA.twins_16S_prep.RData")
load(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/RA_twins/16S/code/RA.twins_16S_metabolon.RData")

# Save R script
# Do this step prior to closing R
save.image(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/RA_twins/16S/code/RA.twins_16S_metabolon.RData")

############################################################################
############################################################################
############################################################################

### Load libraries ###
library(reshape2)
library(phyloseq)
library(vegan)
library(ade4)
library(PMCMR)
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

### olink statistics ###

# directory
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/RA_twins/16S/jobs/13_metabolon/"

# read in data table
d <- read.table(file = paste(dir, "metabolon.tsv", sep = ""),
                header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
                na.strings = "NA")

# create dataframes
d.clean <- as.data.frame(d)
d.all <- as.data.frame(d)

#d <- as.data.frame(d)
#d.all <- d[!d$Sibling_pair == "sib_04",] # remove 506 as no data for 507
#d.clean <- d.all[!d.all$Sibling_pair == "sib_10",] # remove 182, 183

offset = 2
nvar = 10
# create tables for storing wilcoxon and ttest results
stats.table.all <- matrix(data = NA, nrow = nvar, ncol = 3)
colnames(stats.table.all) <- c("cytokine", "wilcoxon", "ttest")

stats.table.clean <- matrix(data = NA, nrow = nvar, ncol = 3)
colnames(stats.table.clean) <- c("cytokine", "wilcoxon", "ttest")

# calculate wilcoxon and ttest between Unaffected/RA siblings for each metabolite
for (i in 1:nvar) {
  stats.table.all[i,1] <- colnames(d.all)[i+offset]
  stats.table.clean[i,1] <- colnames(d.clean)[i+offset]
  
  stats.table.all[i,2] <- wilcox.test(d.all[,i+offset] ~ Diagnosis, data = d.all, paired = TRUE)$p.value
  stats.table.clean[i,2] <- wilcox.test(d.clean[,i+offset] ~ Diagnosis, data = d.clean, paired = TRUE)$p.value
  
  stats.table.all[i,3] <- t.test(d.all[,i+offset] ~ Diagnosis, data = d.all, paired = TRUE)$p.value
  stats.table.clean[i,3] <- t.test(d.clean[,i+offset] ~ Diagnosis, data = d.clean, paired = TRUE)$p.value
}

# save
ft.all = paste(dir, "metabolon_stats_all.csv", sep = "")
write.csv(file = ft.all, stats.table.all)

ft.clean = paste(dir, "metabolon_stats_clean.csv", sep = "")
write.csv(file = ft.clean, stats.table.clean)

############################################################################
############################################################################
############################################################################

### olink plots - all ###

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 8, color = "black")) +
  theme(legend.title = element_text(size = 10, face = "bold", color = "black"))

# directory for storing files
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/RA_twins/16S/jobs/13_metabolon/all/"

# colors
col1 <- c("#929aab", "#ce2525")
col2 <- c("#f3a333", "#0074e4", "#8f8787")

# line types
line1 <- c("solid", "dashed", "dotted")

# for every metabolite in the table, calculate the change in abundance between siblings
for (i in 1:nvar) {
  
  # metabolite of interest
  m <- colnames(d.all)[i+offset]
  
  # create filenames
  filename_table = paste(m, "all_table.csv", sep = "_")
  filename_box.plot = paste(m, "all_box.plot.pdf", sep = "_")  
  filename_line.plot = paste(m, "all_line.plot.pdf", sep = "_")  
  
  # create categories by which to color lineplot: 1) Increase, 2) Decrease, 3) No change
  # subset and spread dataset into diagnosis columns
  # then calculate delta relative abundance between RA and Unaffected siblings
  d.div <- d.all %>%
    subset(select = c("Sibling_pair", "Diagnosis", m)) %>%
    spread(key = "Diagnosis", value = m) %>%
    mutate(Abundance_sib.pair = (get("2_RA") - get("1_Unaffected"))) %>%
    mutate(
      Change.type_sib.pair = case_when(
        sign(Abundance_sib.pair) == 1 ~ "1_up",
        sign(Abundance_sib.pair) == -1 ~ "2_down",
        TRUE ~ "3_no.change"
      )
    ) 
  
  # merge with original dataset
  d.final <- gather(d.div, "1_Unaffected", "2_RA", key = "Diagnosis", value = "Abundance")
  
  # save dataset
  ft = paste(dir, filename_table, sep = "")
  write.csv(d.final, file = ft)
  
  # plot boxplot
  p <- ggplot(data = d.final, aes(x = Diagnosis, y = Abundance, fill = Diagnosis)) +
    stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                 color = "black", size = 0.8, width = 0.3) +
    stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                 color = "black", size = 0.5, width = 0.5) +
    geom_jitter(width = 0.1, size = 1.5) +
    scale_x_discrete(labels = c("Unaffected", "RA")) +
    scale_fill_manual(values = col1) +      
    xlab(NULL) +
    ylab(paste(m, "abundance (ng/ml)")) +
    bkg
  
  fpb = paste(dir, filename_box.plot, sep = "")
  pdf(file = fpb, height = 4.5, width = 5)
  plot(p)
  dev.off()
  
  # plot line plot
  p <- ggplot(data = d.final, aes(x = Diagnosis, y = Abundance, fill = Diagnosis)) +
    stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                 color = "black", size = 0.8, width = 0.3) +
    stat_summary(fun.data = stats.boxplot, geom = "crossbar",
                 color = "black", size = 0.5, width = 0.5) +
    geom_point(shape = 20, size = 3, color = "black") +
    geom_line(data = subset(d.final), aes(group = Sibling_pair, color = Change.type_sib.pair, linetype = Change.type_sib.pair), size = 0.5) +
    geom_text_repel(data = subset(d.final, Diagnosis == "1_Unaffected"), aes(label = Sibling_pair), 
                    nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
    scale_x_discrete(labels = c("Unaffected", "RA")) +
    scale_fill_manual(values = col1) +
    scale_linetype_manual(values = line1, name = "Abundance change", labels = c("Increase", "Decrease", "No change")) +
    scale_color_manual(values = col2, name = "Abundance change", labels = c("Increase", "Decrease", "No change")) +      
    xlab(NULL) +
    ylab(paste(m, "abundance (ng/ml)")) +
    bkg
  
  fpl = paste(dir, filename_line.plot, sep = "")
  pdf(file = fpl, height = 4.5, width = 5)
  plot(p)
  dev.off()
}

############################################################################
############################################################################
############################################################################

### olink plots - clean ###

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 8, color = "black")) +
  theme(legend.title = element_text(size = 10, face = "bold", color = "black"))

# directory for storing files
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/RA_twins/16S/jobs/13_metabolon/clean/"

# colors
col1 <- c("#929aab", "#ce2525")
col2 <- c("#f3a333", "#0074e4", "#8f8787")
col3 <- c("#0074e4")

# line types
line1 <- c("solid", "dashed", "dotted")
line2 <- c("dashed")

# for every metabolite in the table, calculate the change in abundance between siblings
for (i in 1:nvar) {
  
  # metabolite of interest
  m <- colnames(d.clean)[i+offset]
  
  # create filenames
  filename_table = paste(m, "clean_table.csv", sep = "_")
  filename_box.plot = paste(m, "clean_box.plot.pdf", sep = "_")  
  filename_line.plot = paste(m, "clean_line.plot.pdf", sep = "_")  
  
  # create categories by which to color lineplot: 1) Increase, 2) Decrease, 3) No change
  # subset and spread dataset into diagnosis columns
  # then calculate delta relative abundance between RA and Unaffected siblings
  d.div <- d.clean %>%
    subset(select = c("Sibling_pair", "Diagnosis", m)) %>%
    spread(key = "Diagnosis", value = m) %>%
    mutate(Abundance_sib.pair = (get("2_RA") - get("1_Unaffected"))) %>%
    mutate(
      Change.type_sib.pair = case_when(
        sign(Abundance_sib.pair) == 1 ~ "1_up",
        sign(Abundance_sib.pair) == -1 ~ "2_down",
        TRUE ~ "3_no.change"
      )
    ) 
  
  # merge with original dataset
  d.final <- gather(d.div, "1_Unaffected", "2_RA", key = "Diagnosis", value = "Abundance")
  
  # save dataset
  ft = paste(dir, filename_table, sep = "")
  write.csv(d.final, file = ft)
  
  # plot boxplot
  p <- ggplot(data = d.final, aes(x = Diagnosis, y = Abundance, fill = Diagnosis)) +
    stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                 color = "black", size = 0.8, width = 0.3) +
    stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
                 color = "black", size = 0.5, width = 0.5) +
    geom_jitter(width = 0.1, size = 1.5) +
    scale_x_discrete(labels = c("Unaffected", "RA")) +
    scale_fill_manual(values = col1) +      
    xlab(NULL) +
    ylab(paste(m, "abundance (ng/ml)")) +
    bkg
  
  fpb = paste(dir, filename_box.plot, sep = "")
  pdf(file = fpb, height = 4.5, width = 5)
  plot(p)
  dev.off()
  
  # plot line plot
  if(m == "Octanoate") {
    p <- ggplot(data = d.final, aes(x = Diagnosis, y = Abundance, fill = Diagnosis)) +
      stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                   color = "black", size = 0.8, width = 0.3) +
      stat_summary(fun.data = stats.boxplot, geom = "crossbar",
                   color = "black", size = 0.5, width = 0.5) +
      geom_point(shape = 20, size = 3, color = "black") +
      geom_line(data = subset(d.final), aes(group = Sibling_pair, color = Change.type_sib.pair, linetype = Change.type_sib.pair), size = 0.5) +
      geom_text_repel(data = subset(d.final, Diagnosis == "1_Unaffected"), aes(label = Sibling_pair), 
                      nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
      scale_x_discrete(labels = c("Unaffected", "RA")) +
      scale_fill_manual(values = col1) +
      scale_linetype_manual(values = line2, name = "Abundance change", labels = c("Decrease")) +
      scale_color_manual(values = col3, name = "Abundance change", labels = c("Decrease")) +      
      xlab(NULL) +
      ylab(paste(m, "abundance (ng/ml)")) +
      bkg
  } else {
    p <- ggplot(data = d.final, aes(x = Diagnosis, y = Abundance, fill = Diagnosis)) +
      stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
                   color = "black", size = 0.8, width = 0.3) +
      stat_summary(fun.data = stats.boxplot, geom = "crossbar",
                   color = "black", size = 0.5, width = 0.5) +
      geom_point(shape = 20, size = 3, color = "black") +
      geom_line(data = subset(d.final), aes(group = Sibling_pair, color = Change.type_sib.pair, linetype = Change.type_sib.pair), size = 0.5) +
      geom_text_repel(data = subset(d.final, Diagnosis == "1_Unaffected"), aes(label = Sibling_pair), 
                      nudge_x = -0.5, size = 2, color = "black", segment.alpha = 0.3) +
      scale_x_discrete(labels = c("Unaffected", "RA")) +
      scale_fill_manual(values = col1) +
      scale_linetype_manual(values = line1, name = "Abundance change", labels = c("Increase", "Decrease", "No change")) +
      scale_color_manual(values = col2, name = "Abundance change", labels = c("Increase", "Decrease", "No change")) +      
      xlab(NULL) +
      ylab(paste(m, "abundance (ng/ml)")) +
      bkg
  }
    
  fpl = paste(dir, filename_line.plot, sep = "")
  pdf(file = fpl, height = 4.5, width = 5)
  plot(p)
  dev.off()
}
