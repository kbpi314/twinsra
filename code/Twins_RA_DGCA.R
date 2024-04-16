# load libraries
library(DGCA)
library(ggplot2)

# set path
input_df_A = '/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/df_asv_test.tsv'
input_df_B = '/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/df_quant_test.tsv'
out_dir = '/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs07'
corrType = 'pearson'

df_complete = read.table(input_df, sep='\t', header=TRUE, row.names='SampleID')

# obtain sample names
s_names = rownames(df_complete)
v_names = colnames(df_complete)

# create design matrix
patient = ifelse(df_complete$Diagnosis == 'RA',1,0)
control = ifelse(df_complete$Diagnosis == 'Unaffected',1,0)

df_complete = subset(df_complete, select = -c(Diagnosis))

#convert to df
df_complete = t(df_complete)
df_complete = data.frame(df_complete)
# make numeric if needed
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
# df_complete <- data.frame(lapply(df_complete, as.numeric.factor))

# create design matrix 
desmat = data.frame(patient = as.numeric(patient), control = as.numeric(control))
desmat = data.matrix(desmat)

v_ids = rownames(df_complete)

# create dgca result
# obtain all correlations to a specific bacteria
x = data.frame()
for (v in v_ids){
  dgca = ddcorAll(inputMat = df_complete, design = desmat, corrType = corrType, compare = c('patient', 'control'), 
                  adjust = 'perm', nPerm = 10, splitSet = v)
  x = rbind(x, dgca)
}

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

y = completeFun(x, c('Gene1','Gene2'))
# y <- y[1:1000,]
#y = data.frame(lapply(y, as.numeric))
for (i in 1:nrow(y))
{
  y[i,1:2] = sort(y[i,1:2])
}
y = y[!duplicated(y[,1:2]),]

y$'fdr' = p.adjust(y$pValDiff, method = 'fdr')
write.csv(y,paste(out_dir, 'dgca_results_cohort.csv',sep='/'), row.names = FALSE)

y = read.table(paste(out_dir, 'dgca_results_cohort.csv', sep='/'), sep = ',', header = T)
z = y
z = z[order(z$fdr),]
sig_z = subset(z, fdr <= 0.05)
sig_z2 = subset(z, pValDiff <= 0.05)
dir.create(paste(out_dir, 'plots_nofdr', sep='/')) 
# margin(t = 0, r = 0, b = 0, l = 0, unit = "in")
for (i in 1:nrow(sig_z2)){
  p = plotCors(inputMat = df_complete, design = desmat, compare = c('patient', 'control'),
           geneA = sig_z2[i,]$Gene1,
           geneB = sig_z2[i,]$Gene2,
           xlab = v_names[as.numeric(sig_z2[i,]$Gene1)],
           ylab = v_names[as.numeric(sig_z2[i,]$Gene2)]) + 
          theme(axis.text=element_text(size=10),
           axis.title=element_text(size=16),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16))# theme(plot.margin = unit(c(1,1,1,1), "cm")) # + scale_x_continuous(labels = function(x) format(x, scientific = TRUE))
  ggsave(filename = paste(out_dir, '/plots_nofdr/',sig_z2[i,]$Gene1,'_',sig_z2[i,]$Gene2,'.pdf', sep = ''),
         width=6, height=4, units='in', plot = p)
}

