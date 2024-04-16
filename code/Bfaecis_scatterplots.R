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

### B faecis ###
df = read.table('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs02/df_bf.tsv', 
                sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, na.strings = "NA")


# background theme
bkg <- theme_bw() +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
  theme(axis.text.y = element_text(size = 10, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(size = 12, color = "black", face = "bold")) +
  theme(axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 8, color = "black", face = "bold")) +
  theme(legend.title = element_text(size = 10, face = "bold", color = "black"))



#model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE), data=df2)
#coef(model)

#ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
#  geom_point(size=3)+#, color=df$color)+# shape=1 + 
#  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
#  theme(text = element_text(size=16))+theme(legend.position = "none")+
#  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
#  labs(x = X, y = Y) +
#  geom_abline(#data = dfc, 
#    aes(intercept =  10512428.4, slope = 409109.5),
#    #              aes(intercept = `(Intercept)`, slope = dfc),  
#    color = "darkorange2") +
#  geom_abline(#data = dfc, 
#    aes(intercept = 18221970.0, slope = -153429.1 ),
#    #              aes(intercept = `(Intercept)`, slope = dfc),  
#    color = "dodgerblue3")


for (y in c('Flt3L','Hu_IL16','SCF')){
  p <- ggplot(df, aes_string(x='Blautia_faecis', y=y,color='Diagnosis')) + 
    geom_point()+scale_color_manual(values=c("#ce2525","#929aab"))+
    geom_smooth(data=subset(df,Diagnosis=='RA'),method=lm, color="#ce2525")+
    geom_smooth(data=subset(df,Diagnosis=='Unaffected'),method=lm, color="#929aab")+
    ylab(y) +
    xlab("Blautia_faecis") + 
    bkg
  
  filename_box.plot = paste(y,"bfaecis.plot.pdf", sep = "_")  
  fpb = paste(dir, filename_box.plot, sep = "")
  
  pdf(file = fpb, height = 5, width = 5)
  plot(p)
  dev.off()
}

