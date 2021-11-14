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


data9 <- read.delim("data/9_1000_0.tsv", header=T, row.names=1)
meta9 <- read.delim("meta/9_metadata.tsv", header=T, row.names=1)

all(colnames(data9) %in% rownames(meta9))
all(colnames(data9) == rownames(meta9))


dds9 <- DESeqDataSetFromMatrix(countData=data9, 
                               colData=meta9, 
                               design=~condition)
dds9


dds9 <- estimateSizeFactors(dds9)
sizeFactors(dds9)

normalized_counts <- counts(dds9, normalized=TRUE)


dds9 <- DESeqDataSetFromMatrix(countData = data9, colData = meta9, design = ~ condition)
dds9 <- DESeq(dds9)

contrast_oe9 <- c("condition", "condition1", "condition2")

res9 <- results(dds9, contrast = contrast_oe9, alpha = 0.1, lfcThreshold = 0)
summary(res9)

res9 <- res9[!is.na(res9$padj),]
resSig9 <- res9[ res9$padj < 0.1, ]

up9 <- subset(resSig9, log2FoldChange > 0)
down9 <- subset(resSig9, log2FoldChange < 0)

write.csv(up9, "data/results_up9.csv", row.names=TRUE)

write.csv(down9, "data/results_down9.csv", row.names=TRUE)

write.csv(resSig9, "data/results9.csv", row.names=TRUE)
