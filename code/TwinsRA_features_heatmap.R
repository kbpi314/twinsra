# load libraries
library(devtools)
library('ComplexHeatmap')
library(grid)
library(caret)
library(circlize)
library(MASS)

# first time install:
# install_github("jokergoo/ComplexHeatmap")

setwd('~/Desktop/clemente_lab/Projects/twinsra/')

dfs = c('olink', 'fa', 'acpa_fecal', 'acpa_plasma', 'mb', 'quant_R', 'asv_filt', 'path_filt')
# dfs = c('olink')
# dfs = c('acpa_fecal')
# df_quant_meta is just a copy of df_quant
for(i in 1:length(dfs)){
  # grab string
  df = dfs[i]
  
  # load df of features w/metadata
  df3 = read.table(file=paste('inputs/df_',df,'_meta.tsv',sep=''),sep='\t',header=T)#,row.names=1)#,header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  
  # df <- read.csv("/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/df_quant.tsv", sep='\t')
  
  # drop sample for FA
  if (df == 'fa'){
    df3 <- df3[-8,]
  }
  
  
  # rename rows as first col
  df2 <- df3[,-1]
  rownames(df2) <- df3[,1]
  
  # extract feature df and convert to numeric matrix (IMPT!) before hclust
  # useRaster must be true for speedy results
  # features only
  df_feat <- df2[4:(length(df2))]
  # metadata only
  df_meta <- df2[1:3]
  # playtesting dataset
  # df_feat <- df2[1:20]
  
  # normalize columns
  df_feat_norm<- scale(df_feat)
  
  # convert to matrix
  mat_feat <- data.matrix(df_feat_norm)

  # grab feature annotations  
  # rename misleading annotations
  df4 = read.table(file=paste('inputs/df_quant_R_labels.tsv',sep=''),sep='\t',header=T)#,row.names=1)#,header=T, row.names=1, stringsAsFactors=F, check.names=F,comment.char='@')
  cols_anno = df4$var_type
  
  # remove nearzerovar
  col_idx <- nearZeroVar(mat_feat)
  print(df)
  print(length(col_idx))
  if (identical(col_idx, integer(0))){
    mat_feat <- mat_feat
  } else {
    mat_feat <- mat_feat[,-c(col_idx)]
    cols_anno <- cols_anno[-c(col_idx)]
  }

  write.matrix(mat_feat,file=paste('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs19/mat_',df,'_meta.tsv',sep=''), sep='\t')
  
  column_ha = columnAnnotation(VarType=cols_anno)# = df_pc$pc1, bar1 = anno_barplot(df_pc$pc2))

  col_fun_SJC = colorRamp2(c(0, 25), c("white", "red"))
  col_fun_TJC = colorRamp2(c(0, 27), c("white", "red"))
  row_ha = rowAnnotation(Diagnosis = df3$Diagnosis,
                         SJC = df_meta$SJC,
                         TJC = df_meta$TJC,   
                         col = list(Diagnosis = c("RA" = "red", "Unaffected" = "grey"),
                                    SJC = col_fun_SJC, #c(1:max(df_meta$SJC)),
                                    TJC = col_fun_TJC)) #c(1:max(df_meta$TJC))))

  file_path <- paste('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs19/heatmap_',df,'_meta.pdf',sep='')
  
  
  pdf(GetoptLong::qq(file_path), width = 8, height=6)  # Image size
  
  if (df == 'fa'){
    hm <- Heatmap(mat_feat, use_raster = FALSE,
            na_col = 'black',
            name = 'Value', 
            right_annotation = row_ha, 
            show_column_names = TRUE, 
            show_row_names = TRUE, # show sample names
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 6),
            row_title = "Samples (n = 15)",
            column_title = "Features",# (|X| = 9717)",#9717 from 9727 after near zero var) 17 for sig, 9436 for df_final, 760 for asv, 8676 for quant
            cluster_rows = TRUE,
            cluster_columns = TRUE)#,
            #clustering_method_columns = 'average',
            #clustering_method_rows = 'average')
  } else if (df == 'rbfa'){
    hm <- Heatmap(mat_feat, use_raster = FALSE,
                  na_col = 'black',
                  name = 'Value', 
                  right_annotation = row_ha, 
                  show_column_names = TRUE, # there are too many features to read
                  show_row_names = TRUE, # show sample names
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 6),
                  row_title = "Samples (n = 16)",
                  column_title = "Features",# (|X| = 9717)",#9717 from 9727 after near zero var) 17 for sig, 9436 for df_final, 760 for asv, 8676 for quant
                  cluster_rows = TRUE,
                  cluster_columns = TRUE)#,
    #clustering_method_columns = 'average',
    #clustering_method_rows = 'average')
    
  } else if (df == 'mb'){
    hm <- Heatmap(mat_feat, use_raster = FALSE,
                  na_col = 'black',
                  name = 'Value', 
                  right_annotation = row_ha, 
                  show_column_names = TRUE, # there are too many features to read
                  show_row_names = TRUE, # show sample names
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 6),
                  row_title = "Samples (n = 16)",
                  column_title = "Features",# (|X| = 9717)",#9717 from 9727 after near zero var) 17 for sig, 9436 for df_final, 760 for asv, 8676 for quant
                  cluster_rows = TRUE,
                  cluster_columns = TRUE)#,
    #clustering_method_columns = 'average',
    #clustering_method_rows = 'average')
    
  } else if (df == 'quant_R'){
    hm <- Heatmap(mat_feat, use_raster = FALSE,
                  na_col = 'black',
                  name = 'Value', 
                  right_annotation = row_ha, 
                  top_annotation = column_ha,
                  show_column_names = FALSE, # there are too many features to read
                  show_row_names = TRUE, # show sample names
                  row_names_gp = gpar(fontsize = 6),
                  row_title = "Samples (n = 16)",
                  column_title = "Features",# (|X| = 9717)",#9717 from 9727 after near zero var) 17 for sig, 9436 for df_final, 760 for asv, 8676 for quant
                  cluster_rows = TRUE,
                  cluster_columns = TRUE)#,
    #clustering_method_columns = 'average',
    #clustering_method_rows = 'average')
    
  } else {
    hm <- Heatmap(mat_feat, use_raster = FALSE,
                  na_col = 'black',
                  name = 'Value', 
                  right_annotation = row_ha, 
                  show_column_names = TRUE, # there are too many features to read
                  show_row_names = TRUE, # show sample names
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 2),
                  row_title = "Samples (n = 16)",
                  column_title = "Features",# (|X| = 9717)",#9717 from 9727 after near zero var) 17 for sig, 9436 for df_final, 760 for asv, 8676 for quant
                  cluster_rows = TRUE,
                  cluster_columns = TRUE)#,
    #clustering_method_columns = 'average',
    #clustering_method_rows = 'average')
    
  }
  
  
  hm = draw(
    hm, background = "transparent"
    #heatmap_legend_side="bottom", 
    #annotation_legend_side="bottom"
  )  
  dev.off()
}
