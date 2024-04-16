# libraries ---------------------------------------------------------------
library(plyr)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)
library(gridtext)

# Load data ---------------------------------------------------------------
setwd('~/Desktop/clemente_lab/Projects/twinsra/')

dfs = c('olink', 'fa', 'acpa_fecal', 'acpa_plasma', 'rbfa', 'mb', 'plasma', 'quant')
# dfs = c('fa')

for(i in 1:length(dfs)){
  metadf = dfs[i]
  metadata = read.table(file=paste('inputs/df_',metadf,'.tsv',sep=''),sep='\t',header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  
  # Compile Cutie results and calculate corr difference ----------------
  res_dir = '~/Desktop/clemente_lab/Projects/twinsra/outputs/jobs19'
  cutie_res = {}
  for (i in c('UA', 'RA')) {
    fp = sprintf('%s/summary_df_%s.tsv', res_dir, i)
    df = read.table(fp, sep='\t', header=T, 
                    stringsAsFactors=F, check.names=F)
    # filters to keep variables that are in the metadata file; make sure no OTU here
    filt = sapply(df[1:2], function(x) x %in% c(colnames(metadata)))
    
    # keep columns 1, 2 and 4 and keep rows where one of the two entries is in the metadata file list of vars
    df_filt = df[apply(filt, 1, sum)==1, c(1:2, 4)]
    corr_matrix = spread(df_filt, var2, correlations) %>%
      'rownames<-'(.[, 1]) %>% .[, -1]
    cutie_res[[i]] = corr_matrix
  }
  corr_diff = cutie_res[['RA']] - cutie_res[['UA']]
  #colnames(corr_diff)[c(17, 23, 33)] = c('IFN-g', 'IL-1b', 'TNF-a')
  
  # Heatmap -----------------------------------------------------------------
  # Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html
  corr_diff <- corr_diff[,colSums(is.na(corr_diff))<nrow(corr_diff)]
  corr_diff <- corr_diff[rowSums(is.na(corr_diff))<ncol(corr_diff),]
  corr_diff[is.na(corr_diff)] <- 0
  data_matrix = as.matrix(corr_diff) # %>% na.omit
  # is replacing with 0 valid?
  #dend = as.dendrogram(hclust(dist(data_matrix), 
  #                            method='average'))
  
  file_path <- paste('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs19/heatmap_corr_diff_',metadf,'.pdf',sep='')
  pdf(GetoptLong::qq(file_path), width = 8, height=6)  # Image size
  hm <- Heatmap(data_matrix, use_raster = FALSE,
                width = unit(8, "cm"),
                row_dend_width = unit(1.5, "cm"),
                
                clustering_method_rows = "average",  # UPGMA
                clustering_method_columns = "average",
                
                row_title = 'Taxa',
                row_names_gp = gpar(fontsize = 4),
                show_row_names = TRUE,
                
                column_title = 'Quant Vars',
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 4),
                
                
                
                heatmap_legend_param = list(
                  title = "Delta rho (RA - UA)",
                  direction = "horizontal"
                )
  )
  hm = draw(
    hm, background = "transparent",
    heatmap_legend_side="bottom", 
    annotation_legend_side="bottom"
  )
  dev.off()
}

# UA only
for(i in 1:length(dfs)){
  metadf = dfs[i]
  metadata = read.table(file=paste('inputs/df_',metadf,'.tsv',sep=''),sep='\t',header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  
  # Compile Cutie results and calculate corr difference ----------------
  res_dir = '~/Desktop/clemente_lab/Projects/twinsra/outputs/jobs19'
  cutie_res = {}
  for (i in c('UA', 'RA')) {
    fp = sprintf('%s/summary_df_%s.tsv', res_dir, i)
    df = read.table(fp, sep='\t', header=T, 
                    stringsAsFactors=F, check.names=F)
    # filters to keep variables that are in the metadata file; make sure no OTU here
    filt = sapply(df[1:2], function(x) x %in% c(colnames(metadata)))
    
    # keep columns 1, 2 and 4 and keep rows where one of the two entries is in the metadata file list of vars
    df_filt = df[apply(filt, 1, sum)==1, c(1:2, 4)]
    corr_matrix = spread(df_filt, var2, correlations) %>%
      'rownames<-'(.[, 1]) %>% .[, -1]
    cutie_res[[i]] = corr_matrix
  }
  corr_diff = cutie_res[['UA']]
  #colnames(corr_diff)[c(17, 23, 33)] = c('IFN-g', 'IL-1b', 'TNF-a')
  
  # Heatmap -----------------------------------------------------------------
  # Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html
  corr_diff <- corr_diff[,colSums(is.na(corr_diff))<nrow(corr_diff)]
  corr_diff <- corr_diff[rowSums(is.na(corr_diff))<ncol(corr_diff),]
  corr_diff[is.na(corr_diff)] <- 0
  data_matrix = as.matrix(corr_diff) # %>% na.omit
  # is replacing with 0 valid?
  #dend = as.dendrogram(hclust(dist(data_matrix), 
  #                            method='average'))
  
  file_path <- paste('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs19/heatmap_corr_UA_',metadf,'.pdf',sep='')
  pdf(GetoptLong::qq(file_path), width = 6, height=6)  # Image size
  hm <- Heatmap(data_matrix, use_raster = FALSE,
                width = unit(8, "cm"),
                row_dend_width = unit(1.5, "cm"),
                
                clustering_method_rows = "average",  # UPGMA
                clustering_method_columns = "average",
                
                row_title = 'Taxa',
                row_names_gp = gpar(fontsize = 4),
                show_row_names = TRUE,
                
                column_title = 'Quant Vars',
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 4),
                
                
                
                heatmap_legend_param = list(
                  title = "Rho (UA)",
                  direction = "horizontal"
                )
  )
  hm = draw(
    hm, background = "transparent",
    heatmap_legend_side="bottom", 
    annotation_legend_side="bottom"
  )
  dev.off()
}

# RA
for(i in 1:length(dfs)){
  metadf = dfs[i]
  metadata = read.table(file=paste('inputs/df_',metadf,'.tsv',sep=''),sep='\t',header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  
  # Compile Cutie results and calculate corr difference ----------------
  res_dir = '~/Desktop/clemente_lab/Projects/twinsra/outputs/jobs19'
  cutie_res = {}
  for (i in c('UA', 'RA')) {
    fp = sprintf('%s/summary_df_%s.tsv', res_dir, i)
    df = read.table(fp, sep='\t', header=T, 
                    stringsAsFactors=F, check.names=F)
    # filters to keep variables that are in the metadata file; make sure no OTU here
    filt = sapply(df[1:2], function(x) x %in% c(colnames(metadata)))
    
    # keep columns 1, 2 and 4 and keep rows where one of the two entries is in the metadata file list of vars
    df_filt = df[apply(filt, 1, sum)==1, c(1:2, 4)]
    corr_matrix = spread(df_filt, var2, correlations) %>%
      'rownames<-'(.[, 1]) %>% .[, -1]
    cutie_res[[i]] = corr_matrix
  }
  corr_diff = cutie_res[['RA']]
  #colnames(corr_diff)[c(17, 23, 33)] = c('IFN-g', 'IL-1b', 'TNF-a')
  
  # Heatmap -----------------------------------------------------------------
  # Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html
  corr_diff <- corr_diff[,colSums(is.na(corr_diff))<nrow(corr_diff)]
  corr_diff <- corr_diff[rowSums(is.na(corr_diff))<ncol(corr_diff),]
  corr_diff[is.na(corr_diff)] <- 0
  data_matrix = as.matrix(corr_diff) # %>% na.omit
  # is replacing with 0 valid?
  #dend = as.dendrogram(hclust(dist(data_matrix), 
  #                            method='average'))
  
  file_path <- paste('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs19/heatmap_corr_RA_',metadf,'.pdf',sep='')
  pdf(GetoptLong::qq(file_path), width = 6, height=6)  # Image size
  hm <- Heatmap(data_matrix, use_raster = FALSE,
                width = unit(8, "cm"),
                row_dend_width = unit(1.5, "cm"),
                
                clustering_method_rows = "average",  # UPGMA
                clustering_method_columns = "average",
                
                row_title = 'Taxa',
                row_names_gp = gpar(fontsize = 4),
                show_row_names = TRUE,
                
                column_title = 'Quant Vars',
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 4),
                
                
                heatmap_legend_param = list(
                  title = "Rho (RA)",
                  direction = "horizontal"
                )
  )
  hm = draw(
    hm, background = "transparent",
    heatmap_legend_side="bottom", 
    annotation_legend_side="bottom"
  )
  dev.off()
}

