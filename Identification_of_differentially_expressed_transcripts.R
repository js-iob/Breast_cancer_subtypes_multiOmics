#Author: K.T.Shreya Parthasarathi
#Date: 12/06/2024
#Purpose: Estimate differentially expressed transcripts in lumA, lumB, HER2, basal and normal-like breast cancer subtypes against adjacent normal samples - (DESeq2)

#Import libraries
library("DESeq2")
library("ggplot2")
setwd('path\\to\\transcriptomic\\data_file')

#Read file
file = read.table("transcriptomic_data_file.txt", sep = '\t', header = TRUE, row.names=1, check.names =TRUE)
dim(file)

#Data transformation to get raw reads required for DESeq2
file = round(2**file)-1
file[is.na(file)] = 0
file = file +1
dim(file)


#Run DESeq2
countdata = as.matrix(file)
condition = factor(c(rep("LumA",378), rep("LumB",131), rep("HER2",40), rep("Basal_like",122), rep("Normal_like",30), rep("Adj_normal",162)))
coldata = data.frame(row.names = colnames(countdata), condition)

ddsFull = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design=~condition)
dds = DESeq(ddsFull)

#for heatmap plotting
#dds_counts = counts(dds, normalized = T)
#write.table(dds_counts, file = 'deseq2_normalized_samples', sep = '\t')

resultsNames(dds)
res_LumA_Adj_normal = results(dds, contrast = c("condition", "LumA","Adj_normal"))
res_LumB_Adj_normal = results(dds, contrast = c("condition", "LumB","Adj_normal"))
res_HER2_Adj_normal = results(dds, contrast = c("condition", "HER2","Adj_normal"))
res_Basal_Adj_normal = results(dds, contrast = c("condition", "Basal_like","Adj_normal"))
res_Normal_Adj_normal = results(dds, contrast = c("condition", "Normal_like","Adj_normal"))

res_LumA_Adj_normal_ordered = res_LumA_Adj_normal[order(res_LumA_Adj_normal$padj),]
res_LumB_Adj_normal_ordered = res_LumB_Adj_normal[order(res_LumB_Adj_normal$padj),]
res_HER2_Adj_normal_ordered = res_HER2_Adj_normal[order(res_HER2_Adj_normal$padj),]
res_Basal_Adj_normal_ordered = res_Basal_Adj_normal[order(res_Basal_Adj_normal$padj),]
res_Normal_Adj_normal_ordered = res_Normal_Adj_normal[order(res_Normal_Adj_normal$padj),]

#Filter and save results
# if log2Foldchange > 0.6 and adjusted pvalue < 0.05, set as "UP" 
# if log2Foldchange < -0.6 and adjusted pvalue < 0.05, set as "DOWN"

#Luminal A
res_LumA_Adj_normal_ordered$diffexpressed <- "NO"
res_LumA_Adj_normal_ordered$diffexpressed[res_LumA_Adj_normal_ordered$log2FoldChange > 0.6 & res_LumA_Adj_normal_ordered$padj < 0.05] <- "UP"
res_LumA_Adj_normal_ordered$diffexpressed[res_LumA_Adj_normal_ordered$log2FoldChange < -0.6 & res_LumA_Adj_normal_ordered$padj < 0.05] <- "DOWN"
res_LumA_Adj_normal_ordered <- cbind(rownames(res_LumA_Adj_normal_ordered), data.frame(res_LumA_Adj_normal_ordered, row.names=NULL))
colnames(res_LumA_Adj_normal_ordered)[1] <- "genes"
write.table(res_LumA_Adj_normal_ordered, file = 'LumA_vs_adj_normal.txt', sep = '\t', row.names=FALSE, col.names=TRUE)

#Luminal B
res_LumB_Adj_normal_ordered$diffexpressed <- "NO"
res_LumB_Adj_normal_ordered$diffexpressed[res_LumB_Adj_normal_ordered$log2FoldChange > 0.6 & res_LumB_Adj_normal_ordered$padj < 0.05] <- "UP"
res_LumB_Adj_normal_ordered$diffexpressed[res_LumB_Adj_normal_ordered$log2FoldChange < -0.6 & res_LumB_Adj_normal_ordered$padj < 0.05] <- "DOWN"
res_LumB_Adj_normal_ordered <- cbind(rownames(res_LumB_Adj_normal_ordered), data.frame(res_LumB_Adj_normal_ordered, row.names=NULL))
colnames(res_LumB_Adj_normal_ordered)[1] <- "genes"
write.table(res_LumB_Adj_normal_ordered, file = 'LumB_vs_adj_normal.txt', sep = '\t', row.names=FALSE, col.names=TRUE)

#HER2
res_HER2_Adj_normal_ordered$diffexpressed <- "NO"
res_HER2_Adj_normal_ordered$diffexpressed[res_HER2_Adj_normal_ordered$log2FoldChange > 0.6 & res_HER2_Adj_normal_ordered$padj < 0.05] <- "UP"
res_HER2_Adj_normal_ordered$diffexpressed[res_HER2_Adj_normal_ordered$log2FoldChange < -0.6 & res_HER2_Adj_normal_ordered$padj < 0.05] <- "DOWN"
res_HER2_Adj_normal_ordered <- cbind(rownames(res_HER2_Adj_normal_ordered), data.frame(res_HER2_Adj_normal_ordered, row.names=NULL))
colnames(res_HER2_Adj_normal_ordered)[1] <- "genes"
write.table(res_HER2_Adj_normal_ordered, file = 'HER2_vs_adj_normal.txt', sep = '\t', row.names=FALSE, col.names=TRUE)

#Basal
res_Basal_Adj_normal_ordered$diffexpressed <- "NO"
res_Basal_Adj_normal_ordered$diffexpressed[res_Basal_Adj_normal_ordered$log2FoldChange > 0.6 & res_Basal_Adj_normal_ordered$padj < 0.05] <- "UP"
res_Basal_Adj_normal_ordered$diffexpressed[res_Basal_Adj_normal_ordered$log2FoldChange < -0.6 & res_Basal_Adj_normal_ordered$padj < 0.05] <- "DOWN"
res_Basal_Adj_normal_ordered <- cbind(rownames(res_Basal_Adj_normal_ordered), data.frame(res_Basal_Adj_normal_ordered, row.names=NULL))
colnames(res_Basal_Adj_normal_ordered)[1] <- "genes"
write.table(res_Basal_Adj_normal_ordered, file = 'Basal_vs_adj_normal.txt', sep = '\t', row.names=FALSE, col.names=TRUE)

#Normal-like
res_Normal_Adj_normal_ordered$diffexpressed <- "NO"
res_Normal_Adj_normal_ordered$diffexpressed[res_Normal_Adj_normal_ordered$log2FoldChange > 0.6 & res_Normal_Adj_normal_ordered$padj < 0.05] <- "UP"
res_Normal_Adj_normal_ordered$diffexpressed[res_Normal_Adj_normal_ordered$log2FoldChange < -0.6 & res_Normal_Adj_normal_ordered$padj < 0.05] <- "DOWN"
res_Normal_Adj_normal_ordered <- cbind(rownames(res_Normal_Adj_normal_ordered), data.frame(res_Normal_Adj_normal_ordered, row.names=NULL))
colnames(res_Normal_Adj_normal_ordered)[1] <- "genes"
write.table(res_Normal_Adj_normal_ordered, file = 'Normal_vs_adj_normal.txt', sep = '\t', row.names=FALSE, col.names=TRUE)






