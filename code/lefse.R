library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
#######
bkg <-
  theme(axis.text.x = element_text(size = 24, color = "black")) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
  theme(axis.title.x = element_text(margin = unit(c(0,0,4,0), "mm"))) +
  theme(axis.title.x = element_text(size = 24, color = "black")) +
  theme(axis.text.y = element_text(size = 24, color = "black")) +
  theme(axis.title.y = element_text(size = 24, color = "black")) +
  theme(axis.title.y = element_text(margin = unit(c(0,4,0,0), "mm"))) +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 24)) +
  theme(legend.text = element_text(size = 24)) +
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
data <-read.table("/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs06/lefse_results.res", header = FALSE, sep = "\t")

names(data) <- c("RawTaxa", "X", "Group", "LDA", "pval")

plot_data <- subset(data, !is.na(data$LDA))
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
#plot_data$Taxa <- sub(".*_s__", '', plot_data$RawTaxa)
plot_data[plot_data$Group == "Control",]$LDA <- -1 * plot_data[plot_data$Group == "Control",]$LDA

# set colors and factors
group.colors <- c(RA = "#CD3414", Unaffected = "#1462CD")
plot_data$Group <- factor(plot_data$Group, levels = c('Unaffected', 'RA'))

# RAvH
pdf("/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs06/lefse_pretty.pdf", width=20, height=20)



p <- ggbarplot(plot_data, x="Taxa", y="LDA", fill="Group", width= 1, color = "white", sort.val = "asc", sort.by.groups=TRUE) +  
  labs(x = "", y = "LDA score", fill="Group") + coord_flip() + 
  #scale_fill_manual(name="Legend", values = c("RA", "Unaffected')")) +
  # scale_fill_manual(values=c("#E69F00",'#B3A98C','#605843')) + bkg # flip around as need be
  # scale_fill_manual(values=c("#B3A98C",'#E69F00','#605843')) + bkg # flip around as need be
  scale_fill_manual(values=group.colors) + bkg
plot(p)
dev.off()
