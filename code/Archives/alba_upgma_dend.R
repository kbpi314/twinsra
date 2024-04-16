
# dendrograms
library(dendextend)
library(colorspace)

nclust=2

# PCA and dendrograms
path = '/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/'

#df = read.csv(file.path(path,"inputs/df_asv_norm.tsv"), header=T, sep='\t', row.names=1)
#dd <- dist(df, method = "euclidean")

dd <- read.csv(file.path(path,"outputs/jobs15/distance-matrix.tsv"), header=T,sep='\t', row.names=1)
dd <- as.dist(dd, diag = FALSE, upper = FALSE)

hc <- hclust(dd, method = "average")

dend <- as.dendrogram(hc)
dend <- color_branches(dend, k=nclust) #, groupLabels=iris_species)
labels_colors(dend) <-
  rainbow_hcl(nclust)[sort_levels_values(
    as.numeric(df2[,c('group')])[order.dendrogram(dend)]
  )]

# We shall add the flower type to the labels:
#labels(dend) <- paste(as.character(df[,c('group')])[order.dendrogram(dend)],
#                      "(",labels(dend),")", 
#                      sep = "")

#dend <- hang.dendrogram(dend,hang_height=0.1)
dend <- set(dend, "labels_cex", 0.6)

par(mar = c(3,3,3,7))
pdf(file.path(path,"outputs/jobs15/dend2.pdf"), height=6, width=4)#, units='in')
plot(dend, 
     main = "", 
     horiz =  TRUE,  nodePar = list(cex = .007))
# legend("topleft", legend =  unique(df2[,c('group')]), fill = rainbow_hcl(4))
dev.off()


