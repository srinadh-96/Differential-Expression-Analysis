install.packages("tidyverse")
install.packages("RColorBrewer")
install.packages("DESeq2")
install.packages("pheatmap")
install.packages("DEGreport")
install.packages("ROCR")


library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ashr)
library(ROCR)


data2 <- read.delim("data/3_750_250.tsv", header=T, row.names=1)
meta2 <- read.delim("meta/3_metadata.tsv", header=T, row.names=1)

all(colnames(data2) %in% rownames(meta2))
all(colnames(data2) == rownames(meta2))


dds2 <- DESeqDataSetFromMatrix(countData=data2, 
                               colData=meta2, 
                               design=~condition)
dds2


dds2 <- estimateSizeFactors(dds2)
sizeFactors(dds2)

normalized_counts <- counts(dds2, normalized=TRUE)


dds2 <- DESeqDataSetFromMatrix(countData = data2, colData = meta2, design = ~ condition)
dds2 <- DESeq(dds2)

contrast_oe2 <- c("condition", "condition1", "condition2")

res2 <- results(dds2, contrast = contrast_oe2, alpha = 0.1, lfcThreshold = 0)
summary(res2)

res2 <- res2[!is.na(res2$padj),]
resSig2 <- res2[ res2$padj < 0.1, ]

up2 <- subset(resSig2, log2FoldChange > 0)
down2 <- subset(resSig2, log2FoldChange < 0)

write.csv(up2, "data/results_up2.csv", row.names=TRUE)

write.csv(down2, "data/results_down2.csv", row.names=TRUE)

write.csv(resSig2, "data/results2.csv", row.names=TRUE)
