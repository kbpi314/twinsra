# load in alpha div
df_alpha <- read.csv('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs08/alpha_shannon.tsv', 
                     sep='\t')

# load libraries for power analysis
library(pwr)
library(lsr)

# get alpha diffs
alpha_RA = df_alpha[df_alpha$Diagnosis == 'RA',]$shannon_entropy
alpha_UA = df_alpha[df_alpha$Diagnosis == 'Unaffected',]$shannon_entropy
alpha_diffs = as.numeric(alpha_RA) - as.numeric(alpha_UA)

# assume parametric even though is nonparametric setting
cd = cohensD(alpha_diffs)

# above line returns 0.0116015
pwr.t.test(d=cd,sig.level=0.05,power=0.80)

# n=116630 