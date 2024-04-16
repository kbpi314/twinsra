# libraries ---------------------------------------------------------------
library(plyr); library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)
library(gridtext)

# Load data ---------------------------------------------------------------
setwd('~/Clemente/Etc/Alba_heatmap/')
# shotgun = read.table(
#     'data/rel-spp-table-noclades.tsv', sep='\t', header=T, 
#     row.names=1, stringsAsFactors=F, check.names=F
# ) %>% t
metadata = read.table(
    'data/continuous-variables-metadata.txt', sep='\t', header=T, 
    row.names=1, stringsAsFactors=F, check.names=F, 
    comment.char='@'
)

# Create Cutie inputs -----------------------------------------------------
# cutie_data_input = merge(metadata, shotgun, by=0, all=T) %>%
#     .[, -c(2:3, 5:8)]
# colnames(cutie_data_input)[1] = 'SampleID'
# table(cutie_data_input$PsAType)
# 
# # Mixed, Peripheral, Axial
# for (i in unique(cutie_data_input$PsAType)) {
#     print(i)
#     df = cutie_data_input %>% .[.$PsAType==i, ] %>% within(rm(PsAType))
#     fp = sprintf('cutie_data_input_%s.csv', i)
#     write.csv(df, fp, row.names=F)
# }
# # Combined Axial and Mixed
# df = cutie_data_input %>% .[.$PsAType!='Peripheral', ] %>% within(rm(PsAType))
# fp = 'cutie_data_input_Mixed_Axial.csv'
# write.csv(df, fp, row.names=F)

# Compile Cutie results and calculate corr difference ----------------
res_dir = 'cutie/outputs/spearman'
cutie_res = {}
for (i in c('Mixed', 'Peripheral')) {
    fp = sprintf('%s/cutie_summary_%s.txt', res_dir, i)
    df = read.table(fp, sep='\t', header=T, 
                    stringsAsFactors=F, check.names=F)
    filt = sapply(df[1:2], function(x) x %in% colnames(metadata))
        # All taxa are in var1
    df_filt = df[apply(filt, 1, sum)==1, c(1:2, 4)]
    corr_matrix = spread(df_filt, var2, correlations) %>%
        'rownames<-'(.[, 1]) %>% .[, -1]
    cutie_res[[i]] = corr_matrix
}
corr_diff = cutie_res[['Mixed']] - cutie_res[['Peripheral']]
colnames(corr_diff)[c(17, 23, 33)] = c('IFN-g', 'IL-1b', 'TNF-a')

# Heatmap -----------------------------------------------------------------
# Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html
data_matrix = as.matrix(corr_diff) %>% na.omit
dend = as.dendrogram(hclust(dist(data_matrix), 
                            method='average'))

file_path <- 'plots/heatmap_diff_v1.pdf'
pdf(GetoptLong::qq(file_path), width = 6, height=6)  # Image size
hm <- Heatmap(data_matrix, use_raster = FALSE,
    width = unit(8, "cm"),
    row_dend_width = unit(1.5, "cm"),
    
    clustering_method_rows = "average",  # UPGMA
    clustering_method_columns = "average",
    
    row_title = NULL,
    # row_names_gp = gpar(fontsize = 8.5),
    show_row_names = FALSE,
    column_title = NULL,
    column_names_gp = gpar(fontsize = 7),
    
    heatmap_legend_param = list(
        title = "Delta rho (Mixed - Peripheral)",
        direction = "horizontal"
    )
)
hm = draw(
    hm, background = "transparent",
    heatmap_legend_side="bottom", 
    annotation_legend_side="bottom"
)
dev.off()

# Extract dendrogram labels -----------------------------------------------
# Install Shiny
library(devtools)
install_github("jokergoo/InteractiveComplexHeatmap")
library(InteractiveComplexHeatmap)

# Open interactive browser to zoom in on heatmap labels
htShiny(hm)
