#Author: K.T.Shreya Parthasarathi
#Date: 12/06/2024
#Purpose: Boxplots of the ion channels suggestive in prognosis in specific subtype of breast cancer

setwd('path\\to\\transcriptomic\\data_file')
library(ggplot2)
library(tidyr)

# Load the data 
df <- read.csv("ion_channel_expression.txt", sep = '\t', row.names = 1)

# Define sample groups
tumor_groups = c(rep("Normal_like", 30), rep("Adj_normal", 162))

# Transpose the dataframe to get samples as rows
df_t <- as.data.frame(t(df))
df_t$Sample <- rownames(df_t)
df_t$TumorType <- tumor_groups
df_t$TumorType <- factor(df_t$TumorType, levels = c("Normal_like", "Adj_normal"))


# Convert to long format for ggplot
df_long <- df_t %>%
  pivot_longer(cols = -c(Sample, TumorType), names_to = "Gene", values_to = "Expression")

# Create the box plot
ggplot(df_long, aes(x = TumorType, y = Expression, fill = TumorType)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Sample type",
       y = "Gene Expression") +
  theme(legend.position = "none")
