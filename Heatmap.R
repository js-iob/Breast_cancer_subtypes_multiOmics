#Author: K.T.Shreya Parthasarathi
#Date: 12/06/2024
#Purpose: Heatmaps for the deregulated transcripts, differentially methylated probes and amplified/deleted genes across breast cancer subtypes

setwd('path\\to\\transcriptomic\\data_file')

#Load libraries
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(ggplot2)
#library(gridSVG)

# Read file
data_matrix <- read.table('deseq2_normalised_transcript_data.csv', sep = ',',row.names=1,header=TRUE)
dim(data_matrix)
dm = t(apply(data_matrix,1,scale))
colnames(dm) = colnames(data_matrix)

# Define the tumor groups
tumor_groups <- c(rep("LumA", 328), rep("LumB", 131), rep("HER2", 40), rep("Basal", 122), rep("Normal-like", 30), rep("Adj_normal",162))

# Annotation data
annotation_data <- data.frame(TumorType = tumor_groups)
rownames(annotation_data) <- colnames(dm)

# Annotation colors
annotation_colors <- list(TumorType = c("LumA" = "blue", "LumB" = "cyan", "HER2" = "yellow", "Basal" = "orange", "Normal-like" = "purple", "Adj_normal" = "pink"))

# Create annotation object
column_annotation <- HeatmapAnnotation(df = annotation_data, col = annotation_colors)

# Set color function
col_fun <- colorRamp2(c(-1, 0, 1), c("green", "white", "red"))

cell_function <- function(j, i, x, y, width, height, fill) {
  grid.rect(x, y, width, height, gp = gpar(col = "grey", fill = NA))
}

# Create heatmaps for each tumor group
tumor1_data <- dm[, tumor_groups == "LumA"]
tumor1_heatmap <- Heatmap(tumor1_data, name = "LumA", top_annotation = column_annotation[1:329], show_column_names = FALSE, cluster_rows = FALSE, col = col_fun, cluster_columns = TRUE,cell_fun = cell_function)

tumor2_data <- dm[, tumor_groups == "LumB"]
tumor2_heatmap <- Heatmap(tumor2_data, name = "LumB", top_annotation = column_annotation[330:460], show_column_names = FALSE, cluster_rows = FALSE, col = col_fun, cluster_columns = TRUE,cell_fun = cell_function)

tumor3_data <- dm[, tumor_groups == "HER2"]
tumor3_heatmap <- Heatmap(tumor3_data, name = "HER2", top_annotation = column_annotation[461:500], show_column_names = FALSE, cluster_rows = FALSE, col = col_fun, cluster_columns = TRUE, cell_fun = cell_function)

tumor4_data <- dm[, tumor_groups == "Basal"]
tumor4_heatmap <- Heatmap(tumor4_data, name = "Basal", top_annotation = column_annotation[501:622], show_column_names = FALSE, cluster_rows = FALSE, col = col_fun, cluster_columns = TRUE, cell_fun = cell_function)

tumor5_data <- dm[, tumor_groups == "Normal-like"]
tumor5_heatmap <- Heatmap(tumor5_data, name = "Normal-like", top_annotation = column_annotation[623:652], show_column_names = FALSE, cluster_rows = FALSE, col = col_fun, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 8), cell_fun = cell_function)

tumor6_data <- dm[, tumor_groups == "Adj_normal"]
tumor6_heatmap <- Heatmap(tumor6_data, name = "Adj_normal", top_annotation = column_annotation[653:814], show_column_names = FALSE, cluster_rows = FALSE, col = col_fun, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 8), cell_fun = cell_function)

# Combine heatmaps
ht_list <- tumor1_heatmap + tumor2_heatmap + tumor3_heatmap + tumor4_heatmap + tumor5_heatmap + tumor6_heatmap

# Draw combined heatmap with annotation
draw(ht_list, annotation_legend_side = "top", heatmap_legend_side = "right")


