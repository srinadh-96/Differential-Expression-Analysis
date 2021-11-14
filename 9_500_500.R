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


data7 <- read.delim("data/9_500_500.tsv", header=T, row.names=1)
meta7 <- read.delim("meta/9_metadata.tsv", header=T, row.names=1)

all(colnames(data7) %in% rownames(meta7))
all(colnames(data7) == rownames(meta7))


dds7 <- DESeqDataSetFromMatrix(countData=data7, 
                               colData=meta7, 
                               design=~condition)
dds7


dds7 <- estimateSizeFactors(dds7)
sizeFactors(dds7)

normalized_counts <- counts(dds7, normalized=TRUE)


dds7 <- DESeqDataSetFromMatrix(countData = data7, colData = meta7, design = ~ condition)
dds7 <- DESeq(dds7)

contrast_oe7 <- c("condition", "condition1", "condition2")

res7 <- results(dds7, contrast = contrast_oe7, alpha = 0.1, lfcThreshold = 0)
summary(res7)

res7 <- res7[!is.na(res7$padj),]
resSig7 <- res7[ res7$padj < 0.1, ]

up7 <- subset(resSig7, log2FoldChange > 0)
down7 <- subset(resSig7, log2FoldChange < 0)

write.csv(up7, "data/results_up7.csv", row.names=TRUE)

write.csv(down7, "data/results_down7.csv", row.names=TRUE)

write.csv(resSig7, "data/results7.csv", row.names=TRUE)
