################################# 
## R script                    ##
## Project: RA twins           ##
## 16S data                    ##
## Prep for analysis           ##
## Author: JM                  ##
## Date: 5/28/19               ##
#################################

### Load and save current R script ###

# Load R script
load(file="/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/code/RA.twins_16S_prep.RData")

# Save R script
# Do this step prior to closing R
save.image(file="/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/code/RA.twins_16S_prep.RData")

############################################################################
############################################################################
############################################################################

### Load libraries ###

library(phyloseq)
library(ggthemes)
library(ggplot2)

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

### Create phyloseq object ###

# Import biom file
b = "/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/outputs/6_biom_for_R/otu_table_RA.twins_16S_all.json"
biom = import_biom(b, taxaPrefix = F)

# Import mapping file
# Must leave A1 cell of mapping file empty for R compatibility
m = "/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/inputs/map/RA.twins_all_R.txt"
map = sample_data(read.table(m, header = TRUE, sep = "\t", row.names = 1))

# Create phyloseq object
ph = phyloseq(otu_table(biom), tax_table(biom), map)

# Provide column names to separate different taxonomic levels
colnames(tax_table(ph)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Load the tree file (use 'unannotated.tree')
t = "/Users/Lyusik/greengenes/gg_13_8_otus/trees/97_otus_unannotated.tree"
tree = import_qiime(treefilename = t) 

# Merge tree with phyloseq object
phy_16S = merge_phyloseq(ph, map, tree)

# Print row and column names
rownames(sample_data(phy_16S))
colnames(sample_data(phy_16S))

# otu_table()   OTU Table:         [ 10393 taxa and 22 samples ]
# sample_data() Sample Data:       [ 22 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 10393 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 10393 tips and 10392 internal nodes ]

## blank had 0 reads

############################################################################
############################################################################
############################################################################

### Phyloseq object manipulation ###

# remove taxa with zero OTUs
phy_16S_z <- subset_taxa(phy_16S, rowSums(otu_table(phy_16S)) > 0)

# otu_table()   OTU Table:         [ 4035 taxa and 22 samples ]
# sample_data() Sample Data:       [ 22 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 4035 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 4035 tips and 4034 internal nodes ]

# transform from counts to relative abundance
rel_abundance = function(x) { x/sum(x) }
phy_16S_zr <- transform_sample_counts(phy_16S_z, rel_abundance)

############################################################################
############################################################################
############################################################################

### Subset phyloseq objects ###

# Controls #
phy_16S_controls <- subset_samples(phy_16S_z, Diagnosis == "Control")
# 2 samples

# Samples - all #
phy_16S_all <- subset_samples(phy_16S_z, Diagnosis != "Control")
# 20 samples

# Samples - clean - no 182, 183 #
phy_16S_clean <- subset_samples(phy_16S_z, Filter == "keep")
# 18 samples

############################################################################
############################################################################
############################################################################

### Rarefy to even depth ###

# Depth 1000
# Sample w/o replacement
phy_16S_even1000 <- rarefy_even_depth(phy_16S_z, sample.size = 1000, rngseed = 711, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
# 2835 OTUs were removed
# 1200 taxa
# 21 samples

# otu_table()   OTU Table:         [ 1200 taxa and 21 samples ]
# sample_data() Sample Data:       [ 21 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 1200 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1200 tips and 1199 internal nodes ]

###

# Depth 1000
# Sample w/o replacement
phy_16S_even10000 <- rarefy_even_depth(phy_16S_z, sample.size = 10000, rngseed = 711, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
# 1286 OTUs were removed
# 2749 taxa
# 20 samples

# otu_table()   OTU Table:         [ 2749 taxa and 20 samples ]
# sample_data() Sample Data:       [ 20 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 2749 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2749 tips and 2748 internal nodes ]

############################################################################
############################################################################
############################################################################

### Subset rarefied objects ###

# Samples - all - depth 1000 #
phy_16S_all_even1000 <- subset_samples(phy_16S_even1000, Diagnosis != "Control")
# 20 samples

# Samples - clean - no 182, 183 - depth 1000 #
phy_16S_clean_even1000 <- subset_samples(phy_16S_even1000, Filter == "keep")
# 18 samples

##########

# Samples - all - depth 10000 #
phy_16S_all_even10000 <- subset_samples(phy_16S_even10000, Diagnosis != "Control")
# 20 samples

# Samples - clean - no 182, 183 - depth 10000 #
phy_16S_clean_even10000 <- subset_samples(phy_16S_even10000, Filter == "keep")
# 18 samples

############################################################################
############################################################################
############################################################################

### Reads ###

## Setup ##

# directory for storing files
dir = "/Users/Lyusik/Dropbox/Julia_MiCRA/2019_Projects/RA_twins/16S/outputs/7_reads/"

# background theme
bkg <- theme_few() +
  theme(axis.text.x = element_text(size = 16, color = "black", face = "bold")) +
  theme(axis.text.y = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 16, color = "black", face = "bold")) +
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"))) +
  theme(legend.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 11, face = "bold"))

####################

## All samples ##

# create dataset for plotting
d <- data.frame(sample_data(phy_16S_all))   

# create filenames
filename_plot = "16S_all_reads.post.otu.picking_plot.pdf" 
filename_plot.stats = "16S_all_reads.post.otu.picking_plot.stats.csv" 
filename_stats = "16S_all_reads.post.otu.picking_stats.txt"

# plot and save
p <- ggplot(data = d, aes(x = Diagnosis, y = Reads_post_otu_picking, fill = Diagnosis)) +
  stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
               color = "black", size = 0.8, width = 0.3) +
  stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
               color = "black", size = 0.5, width = 0.5, fill = "white") +
  geom_jitter(width = 0.1, size = 1.5) +
  scale_x_discrete(labels = c("Unaffected", "RA")) +
  expand_limits(y = 0) +
  xlab(NULL) + ylab("Reads post-OTU picking") +
  guides(fill = FALSE, color = FALSE) + # no legend
  bkg

fp = paste(dir, filename_plot, sep = "")
pdf(file = fp)
plot(p)
dev.off()

# obtain plot statistics and save
s <- aggregate(Reads_post_otu_picking ~ Diagnosis, data = d, stats.all)
s <- t(s)
rownames(s) <- c("treatment", "mean", "sd", "median", "min", "max", 
                 "10%ile", "25%ile", "75%ile", "90%ile")

fps = paste(dir, filename_plot.stats, sep = "")
write.csv(s, file = fps)

# calculate statistics
wilc <- wilcox.test(Reads_post_otu_picking ~ Diagnosis, data = d, paired = TRUE)

# save calculations
fs = paste(dir, filename_stats, sep = "")
cat("Reads:\n\n", file = fs)
cat("Wilcoxon test:\n", file = fs, append = TRUE)
capture.output(wilc, file = fs, append = TRUE)

####################

## Clean samples - no 182, 183 ##

# create dataset for plotting
d <- data.frame(sample_data(phy_16S_clean))   

# create filenames
filename_plot = "16S_clean_reads.post.otu.picking_plot.pdf" 
filename_plot.stats = "16S_clean_reads.post.otu.picking_plot.stats.csv" 
filename_stats = "16S_clean_reads.post.otu.picking_stats.txt"

# plot and save
p <- ggplot(data = d, aes(x = Diagnosis, y = Reads_post_otu_picking, fill = Diagnosis)) +
  stat_summary(fun.data = stats.whiskers, geom = "errorbar", 
               color = "black", size = 0.8, width = 0.3) +
  stat_summary(fun.data = stats.boxplot, geom = "crossbar", 
               color = "black", size = 0.5, width = 0.5, fill = "white") +
  geom_jitter(width = 0.1, size = 1.5) +
  scale_x_discrete(labels = c("Unaffected", "RA")) +
  expand_limits(y = 0) +
  xlab(NULL) + ylab("Reads post-OTU picking") +
  guides(fill = FALSE, color = FALSE) + # no legend
  bkg

fp = paste(dir, filename_plot, sep = "")
pdf(file = fp)
plot(p)
dev.off()

# obtain plot statistics and save
s <- aggregate(Reads_post_otu_picking ~ Diagnosis, data = d, stats.all)
s <- t(s)
rownames(s) <- c("treatment", "mean", "sd", "median", "min", "max", 
                 "10%ile", "25%ile", "75%ile", "90%ile")

fps = paste(dir, filename_plot.stats, sep = "")
write.csv(s, file = fps)

# calculate statistics
wilc <- wilcox.test(Reads_post_otu_picking ~ Diagnosis, data = d, paired = TRUE)

# save calculations
fs = paste(dir, filename_stats, sep = "")
cat("Reads:\n\n", file = fs)
cat("Wilcoxon test:\n", file = fs, append = TRUE)
capture.output(wilc, file = fs, append = TRUE)

