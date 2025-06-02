# Load required library
library(ggplot2)

# set working dir
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/"

# read in df
df = read.table(paste0(dir, 'inputs/cohens_d_1s_test_df.tsv'), 
                sep = '\t', header = TRUE, row.names = 1, check.names = FALSE,
                na.strings = "NA")

# isolate the subset of data you want to order by
subset_to_order = subset(df, Measure == "Cohens_D_Paired")

# use reorder to reorder the factor
subset_to_order$feature = with(subset_to_order, reorder(feature, Value))

# reorder RA
df$Direction = factor(df$Direction,c('Unaffected','RA'))

# apply the same order to the whole data
df$feature = factor(df$feature, levels = levels(subset_to_order$feature))

col1 <- c("#929aab", "#ce2525")

# flip the bar direction 
df[df$Direction == "Unaffected",]$Value <- -1 * df[df$Direction == "Unaffected",]$Value

p <- ggplot(df, aes(x = Value, y = feature, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cohen's D", y = "Feature") +
  scale_fill_manual(values = col1) +
  
  theme_minimal()

fbp = paste(dir, 'outputs/jobs02/CD.pdf', sep = "")
pdf(file = fbp, height = 5, width = 12)
plot(p)
dev.off()
