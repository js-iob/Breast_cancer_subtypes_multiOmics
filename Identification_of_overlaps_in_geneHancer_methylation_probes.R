#Author: K.T.Shreya Parthasarathi
#Date: 12/06/2024
#Purpose: Identify the overlaps between GeneHancer elements and methylation probes

# Load necessary libraries
library(rtracklayer)
library(GenomicRanges)

setwd('path\\to\\transcriptomic\\data_file')

# Load GeneHancer .gff file
genehancer <- read.csv("GeneHancer_file.txt",sep ='\t')

# Load Illumina methylation probes list
probes <- read.csv("illuminaMethyl450_hg38_GDC.txt",sep = '\t')

# Create GRanges objects
genehancer_gr <- GRanges(seqnames = genehancer$chrom,
                         ranges = IRanges(start = genehancer$start, end = genehancer$end),
                         genehancer_id = genehancer$attributes)

probes_gr <- GRanges(seqnames = probes$chrom,
                     ranges = IRanges(start = probes$chromStart, end = probes$chromEnd),
                     probe_id = probes$id)

# Find overlaps
overlaps <- findOverlaps(probes_gr, genehancer_gr)
overlaps

# Annotate probes with GeneHancer information
annotated_probes <- data.frame(probe_id = probes_gr[queryHits(overlaps)]$probe_id,
                               probe_chromStart = probes_gr[queryHits(overlaps)]$probe_ranges,
                               genehancer_id = genehancer_gr[subjectHits(overlaps)]$genehancer_id)
write.table(annotated_probes,'probes_genehancer_overlaps.txt',sep='\t')
