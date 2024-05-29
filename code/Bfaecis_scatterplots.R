################################# 
## R script                    ##
## Project: TwinsRA            ##
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

# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 8, color = "black")) +
  theme(axis.text.y = element_text(size = 12, color = "black")) +
  theme(axis.title.x = element_text(size = 16, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(size = 16, color = "black", face = "bold")) +
  # theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.title = element_text(size = 16, color = "black", face = "bold")) +
  theme(legend.text = element_text(size = 12, color = "black"))


dfs = c('bf','bw')
taxa = c('Blautia_faecis','Bilophila_wadsworthia')
vars = list(c('Flt3L','SCF'),c('IL7','TNFB'))

for(i in 1:length(dfs)){
  # grab string
  d = dfs[i]
  df = read.table(paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/df_',d,'.tsv'), 
                  sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, na.strings = "NA")
  
  for (y in vars[[i]]){
    p <- ggplot(df, aes_string(x=taxa[i], y=y,color='Diagnosis')) + 
      geom_point()+scale_color_manual(values=c("#ce2525","#929aab"))+
      geom_smooth(data=subset(df,Diagnosis=='RA'),method=lm, color="#ce2525")+
      geom_smooth(data=subset(df,Diagnosis=='Unaffected'),method=lm, color="#929aab")+
      ylab(y) +
      xlab(gsub("_", " ", taxa[i])) + 
      bkg
    
    filename_box.plot = paste(y,taxa[i],"differential.pdf", sep = "_")  
    fpb = paste(dir, filename_box.plot, sep = "")
    
    pdf(file = fpb, height = 5, width = 5)
    plot(p)
    dev.off()
    
    # UA and RA only
    # for (k in c('UA','RA')){
    #   p <- ggplot(df, aes_string(x=taxa[i], y=y,color='Diagnosis')) + 
    #     geom_point()+scale_color_manual(values=c("#ce2525","#929aab"))+
    #     geom_smooth(data=subset(df,Diagnosis=='RA'),method=lm, color="#ce2525")+
    #     geom_smooth(data=subset(df,Diagnosis=='Unaffected'),method=lm, color="#929aab")+
    #     ylab(y) +
    #     xlab(gsub("_", " ", taxa[i])) + 
    #     bkg
    #   
    #   filename_box.plot = paste(y,taxa[i],"differential.pdf", sep = "_")  
    #   fpb = paste(dir, filename_box.plot, sep = "")
    #   
    #   pdf(file = fpb, height = 5, width = 5)
    #   plot(p)
    #   dev.off()
    #   
    # }
  }
}
