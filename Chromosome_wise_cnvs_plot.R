#Author: K.T.Shreya Parthasarathi
#Date: 12/06/2024
#Purpose: Chromosome wise amplifications and deletions observed in each subtype

# Load necessary libraries
library(ggplot2)
library(svglite)

setwd('path\\to\\transcriptomic\\data_file')

data <- read.csv("cnvs_percentage.txt",sep = '\t')

#List of chromosomes
all_chromosomes <- paste0("chr", c(1:22, "X", "Y"))

#Chromosomes present in the data
present_chromosomes <- intersect(all_chromosomes, unique(data$chrom))

# Factorize 'chrom' with the present levels in chronological order
data$chrom <- factor(data$chrom, levels = present_chromosomes)

# Plot
plot <- ggplot(data, aes(x=chrom, y=del_per, group=genes, fill='Deletion')) + 
  geom_bar(stat="identity", position="dodge") + 
  geom_bar(aes(y=amp_per, fill='Amplification'), stat="identity", position="dodge") +
  scale_fill_manual(values = c('Deletion' = 'blue', 'Amplification' = 'red'))+
  labs(x="Chromosome", y="Percentage of Patients") + 
  theme_bw() + 
  theme(fill='white') +   theme(legend.position = "top",
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank())
ggsave("cnv_plot.svg",plot,width=10,height=6,dpi=300)
