#Author: K.T.Shreya Parthasarathi
#Date: 12/06/2024
#Purpose: Volcano plots for the differentially expressed transcripts

#Load libraries
library("ggplot2")
library(ggrepel)

setwd('path\\to\\transcriptomic\\data_file')

df <- read.table('differential_transcripts.txt', sep = '\t', header = TRUE)

#Set theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
plot = ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  ggtitle('Volcano plot') +
  coord_cartesian(ylim = c(0, 250), xlim = c(-10, 10)) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  scale_color_manual(values = c("#00AFBB", "grey", "#d60808"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))
plot
ggsave("volcano_plot.svg",plot,width=8,height=8,dpi=300)




