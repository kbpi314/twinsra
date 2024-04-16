################################
## R script                   ##
## Project: ACTIVATE          ##
## FACS data.                 ##
## Author: AB                 ##
## Date started: 09/15/2023   ##
################################

# Project libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
################################################################################################
################################################################################################
### Themes and colors ###
# background theme
bkg <-
  theme(axis.text.x = element_text(size = 9, color = "black")) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.x = element_text(size = 10, color = "black")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.title.y = element_text(size = 10, color = "black")) +
  theme(axis.title.y = element_text(margin = unit(c(0,4,0,0), "mm"))) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size = 8)) +
  theme(legend.key.size = unit(0.4, 'cm')) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size=12, color= "black", face ="bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(size=6, color= "black"))

col1 <- c("#0a4c7d","#18b5c0")
col2 <- c("#e06666","#93c47d")
################################################################################################
################################################################################################

##### Arranging the data and adding necessary columns #####
setwd('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/code/')
# data <- read.csv('All_counts.tsv', header=T, sep='\t')
data <- read.csv('All_counts_KB.txt', header=T, sep='\t')

#calculate statistics
sumdata <- data %>%
  group_by(DeliveryMode, Specimen) %>%
  summarize(n = n(), mean(Total.cells.cm2,na.rm=TRUE)) 
write.csv(sumdata, 'data_statistics.csv')

#Data is in wide format (one column for dead, one column for live). I need to have one column "Viability"
#I will create a column with a unique ID per sample
data$Sample_ID <- paste(data$PID, data$Specimen, sep="_")
#wide to long
data_long <- gather(data, viability, measurement, Live_sorted_pct, Dead_sorted_pct)
#create a new column with Delivery mode and sample type
data_long$Sample_Delivery <- paste(data_long$DeliveryMode, data_long$Specimen, sep="_")
#create a new column with Delivery mode,sample type and viability
data_long$DeliverySampleViability <- paste(data_long$DeliveryMode, data_long$Specimen, data_long$viability, sep="_")

##### Boxplots summarizing bacterial loads #####
#arrange the order of the specimens(value)
data_long$Specimen <- factor(data_long$Specimen,levels=c('Vaginal_swab', 'Gauze', 'Buccal_swab', 'Forehead_swab'))
pdf(file="Bacterial_load.pdf", width=6, height=4)
p <- ggplot(data_long, aes(x=Specimen, y=log(Total.cells.cm2), fill=DeliveryMode)) + 
  geom_boxplot() +  scale_fill_manual(values = col1) + labs(title = '',
                        x = '', y= "log Bacterial load (cells/cm2)") +bkg
plot(p)
dev.off()

#adding lines to connect samples
pdf(file="Bacterial_load_connected.pdf", width=10, height=6)
data_long$PID <- as.factor(data_long$PID)
p1 <- ggplot(data_long, aes(x=Specimen, y=log(Total.cells.cm2), fill=DeliveryMode)) + 
  geom_boxplot() + geom_point() + geom_line(aes(group=PID, colour=PID)) + facet_wrap(~DeliveryMode, scales="free") + scale_fill_manual(values = col1) + labs(title = '',
                                                            x = '', y= "log Bacterial load (cells/cm2)") + bkg
p <- p1 + theme(legend.position = "none")
plot(p)
dev.off()

#compare groups (pairwise.wilcox.test automatically performs a correction for multiple testing, holmes. To change it : p.adjust.method = "bonferroni")
wx_all = pairwise.wilcox.test(data_long$Total.cells.cm2, data_long$Sample_Delivery)
capture.output(wx_all, file = "group_statistics.csv", sep="/", append = TRUE)

#paired analysis to compare changes across families
#First, calculate if the data is normally distributed
res <- aggregate(cbind(P.value=Total.cells.cm2) ~ Specimen + DeliveryMode, data_long, FUN = function(x) shapiro.test(x)$p.value)
#data is not normal, so I cannot run one-way repeated measures ANOVA. I'll use friedman instead
#Compare bacterial load across samples from the same family, in CS and VD separately

data_long$PID <-as.factor(data_long$PID)
data_long$Total.cells.cm2 <-as.numeric(data_long$Total.cells.cm2)

# 'not an unreplicated complete block design' 
# cs.fried <- subset(data_long, DeliveryMode == 'Csection') %>% friedman_test(Total.cells.cm2 ~ Specimen |PID)

# kb_df for vaginal and c section, default behavior is na is removed
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/friedman.test
kb_df = read.csv(file = "Vaginal_wide_df.csv", row.names="PID")
friedman.test(as.matrix(kb_df))

kb_df = read.csv(file = "Csection_wide_df.csv", row.names="PID")
friedman.test(as.matrix(kb_df))

# anova
kb_long = read.csv(file = "Csection_long_df.csv")#, row.names="PID")

res.aov <- anova_test(data=kb_long,
           dv=logLive,
           wid = PID, 
           within = Specimen)

#Friedman doesn't work, because we don't have complete observations (some families have all sample types, while others don't)
#Do a one-way repeated measures ANOVA
data_long$PID <-as.factor(data_long$PID)
res.aov <- anova_test(data = subset(data_long, DeliveryMode == 'Csection'), dv = Total.cells.cm2, wid = PID, within = Specimen)
get_anova_table(res.aov)


##### Boxplots summarizing Viability (%)  #####

pdf(file="Viability2.pdf", width=6, height=4)
p1 <- ggplot(data_long, aes(x=Specimen, y=measurement, fill=viability)) + facet_wrap(~DeliveryMode, scales="free") +
  geom_boxplot() + labs(title = '',
                                                            x = '', y= "% of sorted population") +bkg
p2 <- p1 + scale_fill_manual(name = "Viability", values=col2, labels = c("Dead", "Live"))
plot(p2)
dev.off()

#compare groups
wx_viab = pairwise.wilcox.test(data_long$measurement, data_long$DeliverySampleViability)
capture.output(wx_viab, file = "viability_statistics.csv", sep="/", append = TRUE)

