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
  theme(axis.title.x = element_text(size = 24, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(size = 24, color = "black")) +
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
data <- read.csv("/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/RA_twins/16S/jobs/15_taxa/15_taxa_stats_clean.csv")
plot_data <- subset(data, data$wilcoxon < 0.05)
# plot_data <- subset(data, RawTaxa %in% 'k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Lachnospiraceae.g__Blautia.s__Blautia_faecis')


taxa_strs <- list()
for (raw in plot_data$RawTaxa) {
  split <- as.character(unlist(str_split(raw, "\\.")))
  i <- length(split)
  blanks <- 0
  while (i > 0) {
    if (split[i] == "__") {
      blanks <- blanks + 1
      print(split)
      split <- split[1:i-1]
      print(split)
    }        
    else {
      break
    }
    i <- i - 1
  }
  
  if (length(split) == 1) {
    taxa_str <- split[1]
  }
  else {
    taxa_str <- paste(split[length(split)-1], split[length(split)])
  }
  if (blanks > 0) {
    for (i in 1:blanks) {
      taxa_str <- paste(taxa_str, "__uncl.", sep="")
    }
  }
  taxa_strs <- append(taxa_strs, taxa_str)
}
plot_data$Taxa <- as.character(taxa_strs)

split_string <- function(x) {
  return(strsplit(x, "s__")[[1]][2])
}

plot_data$Taxa <- sapply(plot_data$Taxa,split_string)
plot_data$Taxa <- gsub("_", " ", plot_data$Taxa)

# plot_data$Taxa <- sub(".*_s__", '', plot_data$RawTaxa)
plot_data[plot_data$Diagnosis == "Control",]$LDA <- -1 * plot_data[plot_data$Diagnosis == "Control",]$LDA

# set colors and factors
Diagnosis.colors <- c(RA = "#CD3414", Unaffected = "#929aab")
plot_data$Diagnosis <- factor(plot_data$Diagnosis, levels = c('Unaffected', 'RA'))

# RAvH
pdf("/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs06/lefse_pretty_filt.pdf", width=8, height=4)

p <- ggbarplot(plot_data, x="Taxa", y="LDA", fill="Diagnosis", width= 1, color = "white", sort.val = "asc", sort.by.Diagnosiss=TRUE) +  
  labs(x = "", y = "LDA score", fill="Diagnosis") + coord_flip() + 
  #scale_fill_manual(name="Legend", values = c("RA", "Unaffected')")) +
  # scale_fill_manual(values=c("#E69F00",'#B3A98C','#605843')) + bkg # flip around as need be
  # scale_fill_manual(values=c("#B3A98C",'#E69F00','#605843')) + bkg # flip around as need be
  scale_fill_manual(values=Diagnosis.colors) + bkg + 
  theme(aspect.ratio = 1/4)
plot(p)
dev.off()
