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


data3 <- read.delim("data/3_1000_0.tsv", header=T, row.names=1)
meta3 <- read.delim("meta/3_metadata.tsv", header=T, row.names=1)

all(colnames(data3) %in% rownames(meta3))
all(colnames(data3) == rownames(meta3))


dds3 <- DESeqDataSetFromMatrix(countData=data3, 
                               colData=meta3, 
                               design=~condition)
dds3


dds3 <- estimateSizeFactors(dds3)
sizeFactors(dds3)

normalized_counts <- counts(dds3, normalized=TRUE)


dds3 <- DESeqDataSetFromMatrix(countData = data3, colData = meta3, design = ~ condition)
dds3 <- DESeq(dds3)

contrast_oe3 <- c("condition", "condition1", "condition2")

res3 <- results(dds3, contrast = contrast_oe3, alpha = 0.1, lfcThreshold = 0)
summary(res3)

res3 <- res3[!is.na(res3$padj),]
resSig3 <- res3[ res3$padj < 0.1, ]

up3 <- subset(resSig3, log2FoldChange > 0)
down3 <- subset(resSig3, log2FoldChange < 0)

write.csv(up3, "data/results_up3.csv", row.names=TRUE)

write.csv(down3, "data/results_down3.csv", row.names=TRUE)

write.csv(resSig3, "data/results3.csv", row.names=TRUE)

