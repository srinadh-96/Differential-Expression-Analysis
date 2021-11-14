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


data4 <- read.delim("data/6_500_500.tsv", header=T, row.names=1)
meta4 <- read.delim("meta/6_metadata.tsv", header=T, row.names=1)

all(colnames(data4) %in% rownames(meta4))
all(colnames(data4) == rownames(meta4))


dds4 <- DESeqDataSetFromMatrix(countData=data4, 
                               colData=meta4, 
                               design=~condition)
dds4


dds4 <- estimateSizeFactors(dds4)
sizeFactors(dds4)

normalized_counts <- counts(dds2, normalized=TRUE)


dds4 <- DESeqDataSetFromMatrix(countData = data4, colData = meta4, design = ~ condition)
dds4 <- DESeq(dds4)

contrast_oe2 <- c("condition", "condition1", "condition2")

res4 <- results(dds4, contrast = contrast_oe4, alpha = 0.1, lfcThreshold = 0)
summary(res4)

res4 <- res4[!is.na(res4$padj),]
resSig4 <- res4[ res4$padj < 0.1, ]

up4 <- subset(resSig4, log2FoldChange > 0)
down4 <- subset(resSig4, log2FoldChange < 0)

write.csv(up4, "data/results_up4.csv", row.names=TRUE)

write.csv(down4, "data/results_down4.csv", row.names=TRUE)

write.csv(resSig4, "data/results4.csv", row.names=TRUE)
