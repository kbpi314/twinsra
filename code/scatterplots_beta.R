################################# 
## R script                    ##
## Project: TwinsRA            ##
## Beta diversity              ##
## Data: Shotgun metagenomics  ##
## Author: KB                  ##
## Date: 3/7/21                ##
#################################

### Load and save current R script ###

# Load R scripts
#load(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs20/alpha_data.RData")
#load(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs20/alpha_script.RData")

# Save R script
# Do this step prior to closing R
#save.image(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs20/alpha_script.RData")

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

# set working dir
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/"

### Alpha Diversity Boxplots ###
df = read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/df_corr3.tsv',#df_beta_meta.tsv', 
                sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, na.strings = "NA")


# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 8, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 10, face = "bold", color = "black"))

# choose colors
col1 <- c("#929aab", "#ce2525")
col2 <- c("#f3a333", "#0074e4", "#8f8787")

for (i in 1:4) {
  # var of interest
  a <- colnames(df)[i+2]
  
  p <- ggplot(df, aes_string(x='Bray_Curtis', y=paste(a))) + 
    geom_point()+
    geom_smooth(method=lm, color="black")+
    ylab(paste(a)) +
    xlab("Bray-Curtis Dist from UA Twin") + 
    bkg
  
  filename_box.plot = paste(a, "beta_clean_corr.plot.pdf", sep = "_")  
  fpb = paste(dir, filename_box.plot, sep = "")
  
  pdf(file = fpb, height = 3.5, width = 3.5)
  plot(p)
  dev.off()
}

