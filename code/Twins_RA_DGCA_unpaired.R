# load libraries
library(DGCA)
library(ggplot2)

# set path
input_df = '/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/df_final.tsv'
input_df_A = '/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/df_asv_test.tsv'
input_df_B = '/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/df_quant_test.tsv'
out_dir = '/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs07'
corrType = 'pearson'

# process DFs
df_complete_A = read.table(input_df_A, sep='\t', header=TRUE, row.names='SampleID')
df_complete_B = read.table(input_df_B, sep='\t', header=TRUE, row.names='SampleID')
df_complete = read.table(input_df, sep='\t', header=TRUE, row.names='SampleID')

# obtain sample names
s_names_A = rownames(df_complete_A)
s_names_B = rownames(df_complete_B)

# create design matrix
patient = ifelse(df_complete_A$Diagnosis == 'RA',1,0)
control = ifelse(df_complete_A$Diagnosis == 'Unaffected',1,0)
desmat = data.frame(patient = as.numeric(patient), control = as.numeric(control))
desmat = data.matrix(desmat)

# remove diagnosis column
df_complete_A = subset(df_complete_A, select = -c(Diagnosis))
df_complete_B = subset(df_complete_B, select = -c(Diagnosis))
df_complete = subset(df_complete, select = -c(Diagnosis))

# obtain variable names
v_names_A = colnames(df_complete_A)
v_names_B = colnames(df_complete_B)

#convert to df untidy because gene expression world is untidy
df_complete_A = t(df_complete_A)
df_complete_A = data.frame(df_complete_A)
df_complete_B = t(df_complete_B)
df_complete_B = data.frame(df_complete_B)
df_complete = t(df_complete)
df_complete = data.frame(df_complete)


# make numeric if needed
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
# df_complete <- data.frame(lapply(df_complete, as.numeric.factor))

# create DGCA result
x = data.frame()
dgca = ddcorAll(inputMat = df_complete_A, 
                inputMatB = df_complete_B,
                design = desmat, corrType = corrType, compare = c('patient', 'control'), 
                adjust = 'perm', nPerm = 10)
x = rbind(x, dgca)

# function for getting complete cases
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

y = completeFun(x, c('Gene1','Gene2'))

# sort and remove duplicates
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
# if doing fdr
sig_z = subset(z, fdr <= 0.05)
# toggle to ignore fdr correction
sig_z2 = subset(z, pValDiff <= 0.05)
dir.create(paste(out_dir, 'plots_nofdr', sep='/')) 
# margin(t = 0, r = 0, b = 0, l = 0, unit = "in")
for (i in 1:nrow(sig_z2)){
  p = plotCors(inputMat = df_complete, design = desmat, compare = c('patient', 'control'),
           geneA = sig_z2[i,]$Gene1,
           geneB = sig_z2[i,]$Gene2,#+#,
           xlab = sig_z2[i,]$Gene1,
           ylab = sig_z2[i,]$Gene2)+#,
           #xlab = v_names_A[as.numeric(sig_z2[i,]$Gene1)],
           #ylab = v_names_B[as.numeric(sig_z2[i,]$Gene2)]) + 
          theme(axis.text=element_text(size=10),
           axis.title=element_text(size=16),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16))# theme(plot.margin = unit(c(1,1,1,1), "cm")) # + scale_x_continuous(labels = function(x) format(x, scientific = TRUE))
  ggsave(filename = paste(out_dir, '/plots_nofdr/',sig_z2[i,]$Gene1,'_',sig_z2[i,]$Gene2,'.pdf', sep = ''),
         width=12, height=8, units='in', plot = p)
}

