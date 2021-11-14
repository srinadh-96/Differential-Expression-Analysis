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


data6 <- read.delim("data/6_1000_0.tsv", header=T, row.names=1)
meta6 <- read.delim("meta/6_metadata.tsv", header=T, row.names=1)

all(colnames(data6) %in% rownames(meta6))
all(colnames(data6) == rownames(meta6))


dds6 <- DESeqDataSetFromMatrix(countData=data6, 
                               colData=meta6, 
                               design=~condition)
dds6


dds6 <- estimateSizeFactors(dds6)
sizeFactors(dds6)

normalized_counts <- counts(dds6, normalized=TRUE)


dds6 <- DESeqDataSetFromMatrix(countData = data6, colData = meta6, design = ~ condition)
dds6 <- DESeq(dds6)

contrast_oe6 <- c("condition", "condition1", "condition2")

res6 <- results(dds6, contrast = contrast_oe6, alpha = 0.1, lfcThreshold = 0)
summary(res6)

res6 <- res6[!is.na(res6$padj),]
resSig6 <- res6[ res6$padj < 0.1, ]

up6 <- subset(resSig6, log2FoldChange > 0)
down6 <- subset(resSig6, log2FoldChange < 0)

write.csv(up6, "data/results_up6.csv", row.names=TRUE)

write.csv(down6, "data/results_down6.csv", row.names=TRUE)

write.csv(resSig6, "data/results6.csv", row.names=TRUE)
