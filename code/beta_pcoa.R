################################# 
## R script                    ##
## Project: RA twins           ##
## Shotgun                     ##
## Beta diversity pcoa plots   ##
## Author: KB                  ##
## Last Updated: 5/9/24        ##
#################################

### Load and save current R script ###

# Load R scripts
# load(file="/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/code/RA.twins_16S_prep.RData")
# load(file="/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/code/RA.twins_16S_beta.div_pcoa.RData")

# Save R script
# Do this step prior to closing R
# save.image(file="/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/code/RA.twins_16S_beta.div_pcoa.RData")

############################################################################
############################################################################
############################################################################

### Load libraries ###

library(phyloseq)
library(vegan)
library(ade4)
library(PMCMR)
library(PMCMRplus)
library(ggplot2)
library(ggthemes)
library(ggrepel)

############################################################################
############################################################################
############################################################################

### Beta diversity pcoa ###

# background theme
bkg <- theme_bw() +
  theme(axis.text = element_text(size = 24, color = "black")) +
  theme(axis.title = element_text(size = 32, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 18, color = "black"))+ #, face = "bold")) +
  # theme(legend.title = element_blank()) +
  theme(legend.title = element_text(size = 24, face = "bold", color = "black"))
  theme(legend.justification = "right")# +
  #theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5))

# function to specify that axis labels have 2 decimal places
f.dec <- function(x){
  format(round(x, 2), nsmall = 2)
}

# directory for storing files
# dir = "/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/jobs/2_beta.div_pcoa/"
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/"

# list of distance methods
# dists <- c("bray", "unifrac", "wunifrac")
dists <- c('unifrac')

# colors
col1 <- c("#929aab", "#ce2525")

# load data
df <- read.delim(file="/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/bray_curtis_pcoa.tsv",
              row.names=1)

# order factors for legend
df$Diagnosis <- factor(df$Diagnosis, levels=c('Unaffected', 'RA'))

for (j in seq_along(dists)) {
  
  # create filenames
  filename_plot = paste("bdiv", dists[j], "plot.pdf", sep = "_")
  
  # plot beta diversity
  p <- ggplot() + # data=df, aes(x = PC1, y = PC2, fill = Diagnosis)) +
    geom_point(data = df, aes(x = PC1, y = PC2, color = Diagnosis),size=4) +
    scale_color_manual(values = c("Unaffected" = "#929aab", "RA" = "#ce2525")) + #col1, labels = c("Unaffected", "RA")) +
    bkg +
    scale_x_continuous(labels = f.dec) + # 2 decimal places on x-axis
    scale_y_continuous(labels = f.dec)   # 2 decimal places on y-axis
  
  # save plot
  fp = paste(dir, filename_plot, sep = "")
  pdf(file = fp, height = 6, width = 8)
  plot(p)
  dev.off()
}
