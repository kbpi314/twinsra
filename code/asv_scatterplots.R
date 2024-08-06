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

# c(BF FLT3L, BF SCF) c(BW IL7, BW TNFB)
UA_pvals = list(c('= 0.002  ','= 0.00003'), c('=NA   ','=NA   '))
RA_pvals = list(c('= 0.45  ','= 0.21   '),c('=0.007  ','=0.002  '))
xpos = list(c(0.003,0.0025),c(0.00025, 0.00025))
ypos = list(c(7.5, 8.5),c(-0.375,3.125))
offset = list(c(0.09, 0.05),c(0.095,0.055))

for(i in 1:length(dfs)){
  # grab string
  d = dfs[i]
  df = read.table(paste0('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/df_',d,'.tsv'), 
                  sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, na.strings = "NA")
  
  df$Diagnosis = factor(df$Diagnosis,levels=c('Unaffected','RA'))
  for (j in 1:length(vars[[i]])){
    y = vars[[i]][j]
    p <- ggplot(df, aes_string(x=taxa[i], y=y,color='Diagnosis')) + 
      geom_point()+scale_color_manual(values=c("#929aab","#ce2525"))+
      geom_smooth(data=subset(df,Diagnosis=='RA'),method=lm, color="#ce2525")+
      geom_smooth(data=subset(df,Diagnosis=='Unaffected'),method=lm, color="#929aab")+
      ylab(y) +
      xlab(gsub("_", " ", taxa[i])) + 
      annotate(geom="text",label=sprintf("italic('p')~'%s'",UA_pvals[[i]][j]),parse=TRUE,x=xpos[[i]][j],y=ypos[[i]][j]+offset[[i]][j],size=6,colour = "#929aab",hjust = 0) + 
      annotate(geom="text",label=sprintf("italic('p')~'%s'",RA_pvals[[i]][j]),parse=TRUE,x=xpos[[i]][j],y=ypos[[i]][j]-offset[[i]][j],size=6,colour = "#ce2525",hjust = 0) +
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
