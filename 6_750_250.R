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


data5 <- read.delim("data/6_750_250.tsv", header=T, row.names=1)
meta5 <- read.delim("meta/6_metadata.tsv", header=T, row.names=1)

all(colnames(data5) %in% rownames(meta5))
all(colnames(data5) == rownames(meta5))


dds5 <- DESeqDataSetFromMatrix(countData=data5, 
                               colData=meta5, 
                               design=~condition)
dds5


dds5 <- estimateSizeFactors(dds5)
sizeFactors(dds5)

normalized_counts <- counts(dds5, normalized=TRUE)


dds5 <- DESeqDataSetFromMatrix(countData = data5, colData = meta5, design = ~ condition)
dds5 <- DESeq(dds5)

contrast_oe5 <- c("condition", "condition1", "condition2")

res5 <- results(dds5, contrast = contrast_oe5, alpha = 0.1, lfcThreshold = 1)
summary(res5)

res5 <- res5[!is.na(res5$padj),]
resSig5 <- res5[ res5$padj < 0.1, ]

up5 <- subset(resSig5, log2FoldChange > 0)
down5 <- subset(resSig5, log2FoldChange < 0)

write.csv(up5, "data/results_up5.csv", row.names=TRUE)

write.csv(down5, "data/results_down5.csv", row.names=TRUE)

write.csv(resSig5, "data/results5.csv", row.names=TRUE)
