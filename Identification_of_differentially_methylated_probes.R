#Author: K.T.Shreya Parthasarathi
#Date: 12/06/2024
#Purpose: Estimate differentially methylated probes in lumA, lumB, HER2, basal and normal-like breast cancer subtypes against adjacent normal samples - (ChAMP)

# Load necessary libraries (if not already done)
library(ChAMP)
library(data.table)
library(minfi)
setwd('path\\to\\transcriptomic\\data_file')

#Load beta matrix
myData <- data.frame(fread("methylomic_data_file.csv",sep = '\t', header = TRUE))
row.names(myData) = myData$Composite.Element.REF
myData2 = myData[ , !names(myData) %in% 
                       c("Composite.Element.REF")]
myData_mat = data.matrix(myData2)

#Load phenotype data
phenotype_data <- read.table("sample_ids_phenotype_data.txt",sep = '\t', header = TRUE)
length(phenotype_data)

head(phenotype_data)
colnames(myData_mat)

#Filtering
filtered_myData_mat = champ.filter(
  beta = myData_mat, pd=phenotype_data, autoimpute = TRUE, arraytype = "450K")

#Quality Check
CpG.GUI(CpG=rownames(filtered_myData_mat$beta),arraytype="450K")
qc_results <- champ.QC(beta = filtered_myData_mat$beta, pheno = phenotype_data$Sample_type)
QC.GUI(beta=filtered_myData_mat$beta, pheno = phenotype_data$Sample_type,arraytype="450K")

#Normalization
myNorm_data <- champ.norm(beta = filtered_myData_mat$beta,arraytype="450K")

#Quality check post normalization
QC.GUI(beta=myNorm_data, pheno = phenotype_data$Sample_type,arraytype="450K")
qc_results <- champ.QC(beta = myNorm_data, pheno = phenotype_data$Sample_type)


#Differentially methylated probes
myDMP <- champ.DMP(beta = myNorm_data, pheno = phenotype_data$Sample_type, adjPVal = 0.05, adjust.method = "BH", arraytype = "450K")
DMP = (myDMP[[1]])
write.table(DMP, file = 'DMP_subtype.txt', sep = '\t')

