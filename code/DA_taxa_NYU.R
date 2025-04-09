################################# 
## R script                    ##
## Project: TwinsRA            ##
## LEfSe                       ##
## Data: Shotgun metagenomics  ##
## Author: KB.                 ##
## Last Updated: 5/9/24        ##
#################################

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

bkg <-
  theme(axis.text.x = element_text(size = 18, color = "black")) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
  theme(axis.title.x = element_text(margin = unit(c(0,0,4,0), "mm"))) +
  theme(axis.title.x = element_text(size = 22, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(size = 18, color = "black")) +
  theme(axis.title.y = element_text(size = 24, color = "black")) +
  theme(axis.title.y = element_text(margin = unit(c(0,4,0,0), "mm"))) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 24, color = "black", face = "bold")) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.key.size = unit(0.6, 'cm')) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size=, color= "black", face ="bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(size=24, color= "black"))
####

# RAvH
data <-read.csv("/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/RA_twins/16S/jobs/15_taxa/15_taxa_stats_clean.csv")
plot_data <- subset(data, data$wilcoxon < 0.05)

#plot_data <- subset(data, RawTaxa %in% 'k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfovibrionales.f__Desulfovibrionaceae.g__Bilophila.s__Bilophila_wadsworthia' |
#                          RawTaxa %in% 'k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Lachnospiraceae.g__Blautia.s__Blautia_faecis')
#plot_data <- subset(data, RawTaxa %in% 'k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Lachnospiraceae.g__Blautia.s__Blautia_faecis')

# plot_data$Taxa <- sub(".*_s__", '', plot_data$RawTaxa)

plot_data$logp = -1*log(plot_data$wilcoxon,base=10)

# negative t statistic is enriched in UA
plot_data$logp <- sign(plot_data$ttest.stat) * plot_data$logp

# create diagnosis
plot_data <- plot_data %>%
  mutate(Diagnosis = case_when(
    ttest.stat > 0 ~ "RA",
    ttest.stat < 0 ~ "Unaffected",
    TRUE ~ "zero"
  ))


# set colors and factors
Diagnosis.colors <- c(RA = "#CD3414", Unaffected = "#929aab")
plot_data$Diagnosis <- factor(plot_data$Diagnosis, levels = c('Unaffected', 'RA'))



# RAvH
pdf("/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs06/lefse_pretty_filt_taxa_NYU.pdf", width=18, height=8)

p <- ggbarplot(plot_data, x="feature", y="logp", fill="Diagnosis", width= 1, color = "white", sort.val = "asc", sort.by.Diagnosiss=TRUE) +  
  labs(x = "", y = "-log10(p) * sgn(enrichment)", fill="Diagnosis") + coord_flip() + 
  #scale_fill_manual(name="Legend", values = c("RA", "Unaffected')")) +
  # scale_fill_manual(values=c("#E69F00",'#B3A98C','#605843')) + bkg # flip around as need be
  # scale_fill_manual(values=c("#B3A98C",'#E69F00','#605843')) + bkg # flip around as need be
  scale_fill_manual(values=Diagnosis.colors) + bkg + 
  theme(aspect.ratio = 1/4)
plot(p)
dev.off()
