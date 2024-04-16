# DMM Tutorial: https://microbiome.github.io/tutorials/DMM.html

# install.packages('BiocManager')
# library(BiocManager)
# BiocManager::install("microbiome")
# BiocManager::install("DirichletMultinomial")
# install.packages("remotes")
# remotes::install_github("microbiome/microbiome")
# install.packages('Cairo')
# install.packages('svglite')

library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(Cairo)

# Load OTU table
base.path <- "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/"
viz.path <- paste(base.path, 'outputs/jobs05/', sep="")
df <- read.csv(paste(base.path, 'outputs/jobs00/metaphlan_taxa_table_abs.tsv', sep=""), sep='\t')

rownames(df) <- df[,1]
df <- df[-c(1)] 
count_df <- as.matrix(df)

# transpose for metaphlan
count_df <- t(count_df)

# Fit the DMM Model, starts being slow at 7+
fit <- lapply(1:6, dmn, count = count_df, verbose=TRUE)

# Check model fit
lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
aic  <- base::sapply(fit, DirichletMultinomial::AIC) # AIC / BIC / Laplace
bic  <- base::sapply(fit, DirichletMultinomial::BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)
df.fit <- data.frame(cbind(lplc, aic, bic))
write.table(df.fit, "DMM_cluster_results.tsv", sep='\t')

# find replace x = 1:8 with 1:6
options(bitmapType='cairo')
p <- ggplot(data = df.fit) +
    geom_point(aes(x = 1:6, y = lplc, color = 'blue')) + geom_line(aes(x = 1:6, y = lplc, color = 'blue')) +
    geom_point(aes(x = 1:6, y = aic, color = 'green')) + geom_line(aes(x = 1:6, y = aic, color = 'green')) +  
    geom_point(aes(x = 1:6, y = bic, color = 'magenta')) + geom_line(aes(x = 1:6, y = bic, color = 'magenta')) +
    labs(x  ="Maximum Number of Community Types", y = "Quality of fit") +
    scale_color_discrete(name = "Method", labels=c('Laplace', 'AIC', 'BIC')) +
    theme_light() +
    theme(legend.position = c(0.2, 0.75),
        legend.box.background = element_rect(color="grey"),
        legend.text = element_text(size = 10),
        text = element_text(size = 14)
    )
show(p)
ggsave(paste(viz.path, "DMM_bestfit_relabd.svg", sep=""), p, width=5, height=5)

# Pick optimal model
lplc_best <- fit[[which.min(unlist(lplc))]]
aic_best <- fit[[which.min(unlist(aic))]]
bic_best <- fit[[which.min(unlist(bic))]]

# Mixture parameters pi and theta
mixturewt(lplc_best)
mixturewt(aic_best)
mixturewt(bic_best)

# Sample-component assignments
lplc_ass <- apply(mixture(lplc_best), 1, which.max)
aic_ass <- apply(mixture(aic_best), 1, which.max)
bic_ass <- apply(mixture(bic_best), 1, which.max)

# Contribution of each taxonomic group to each component
for (k in seq(ncol(fitted(lplc_best)))) {
  d <- melt(fitted(lplc_best))
  colnames(d) <- c("OTU", "cluster", "value")
  e <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))  
    #filter(abs(value) > value_thresh)  
  
  p <- ggplot(e, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top Drivers: Community Type", k))
  
  ggsave(paste(viz.path, "DMM_relabd_thresh_lplc_grp", k, ".svg", sep=""), p, width=15, height=1+4*dim(e)[1]/24)
}

for (k in seq(ncol(fitted(aic_best)))) {
  d <- melt(fitted(aic_best))
  colnames(d) <- c("OTU", "cluster", "value")
  e <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))  
    #filter(abs(value) > value_thresh)  
  
  p <- ggplot(e, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top Drivers: Community Type", k))
  
  ggsave(paste(viz.path, "DMM_relabd_thresh_aic_grp", k, ".svg", sep=""), p, width=15, height=1+4*dim(e)[1]/24)
}

for (k in seq(ncol(fitted(bic_best)))) {
  d <- melt(fitted(bic_best))
  colnames(d) <- c("OTU", "cluster", "value")
  e <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))  
    #filter(abs(value) > value_thresh)  
  
  p <- ggplot(e, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top Drivers: Community Type", k))
  
  ggsave(paste(viz.path, "DMM_relabd_thresh_bic_grp", k, ".svg", sep=""), p, width=15, height=1+4*dim(e)[1]/24)
}
# Extract which samples belong to each community group (the ass variable)

lplc_ass.df <- data.frame(cbind(names(lplc_ass), matrix(unlist(lplc_ass), nrow=length(lplc_ass), byrow=TRUE)))
colnames(lplc_ass.df) <- c("HostSubjectId", "dmm.grp")
write.table(lplc_ass.df, file = paste(viz.path, 'dmm_lplc_relabd_assign.tsv', sep=""), sep='\t', row.names = FALSE)

aic_ass.df <- data.frame(cbind(names(aic_ass), matrix(unlist(aic_ass), nrow=length(aic_ass), byrow=TRUE)))
colnames(aic_ass.df) <- c("HostSubjectId", "dmm.grp")
write.table(aic_ass.df, file = paste(viz.path, 'dmm_aic_relabd_assign.tsv', sep=""), sep='\t', row.names = FALSE)

bic_ass.df <- data.frame(cbind(names(bic_ass), matrix(unlist(bic_ass), nrow=length(bic_ass), byrow=TRUE)))
colnames(bic_ass.df) <- c("HostSubjectId", "dmm.grp")
write.table(bic_ass.df, file = paste(viz.path, 'dmm_bic_relabd_assign.tsv', sep=""), sep='\t', row.names = FALSE)

