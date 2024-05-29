# Load required library
library(ggplot2)

# set working dir
dir = "/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/"

# read in df
df = read.table(paste0(dir, 'inputs/mag_d_1s_test_df.tsv'), 
                sep = '\t', header = TRUE, row.names = 1, check.names = FALSE,
                na.strings = "NA")

# isolate the subset of data you want to order by
subset_to_order = subset(df, Measure == "Mag_D_Paired")

# use reorder to reorder the factor
subset_to_order$feature = with(subset_to_order, reorder(feature, Value))

# apply the same order to the whole data
df$feature = factor(df$feature, levels = levels(subset_to_order$feature))

ggplot(df, aes(x = Value, y = feature, fill = Measure)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Mag's D", y = "Feature") +
  theme_minimal()