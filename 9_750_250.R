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


data8 <- read.delim("data/9_750_250.tsv", header=T, row.names=1)
meta8 <- read.delim("meta/9_metadata.tsv", header=T, row.names=1)

all(colnames(data8) %in% rownames(meta8))
all(colnames(data8) == rownames(meta8))


dds8 <- DESeqDataSetFromMatrix(countData=data8, 
                               colData=meta8, 
                               design=~condition)
dds8


dds8 <- estimateSizeFactors(dds8)
sizeFactors(dds8)

normalized_counts <- counts(dds8, normalized=TRUE)


dds8 <- DESeqDataSetFromMatrix(countData = data8, colData = meta8, design = ~ condition)
dds8 <- DESeq(dds8)

contrast_oe8 <- c("condition", "condition1", "condition2")

res8 <- results(dds8, contrast = contrast_oe8, alpha = 0.1, lfcThreshold = 0)
summary(res8)

res8 <- res8[!is.na(res8$padj),]
resSig8 <- res8[ res8$padj < 0.1, ]

up8 <- subset(resSig8, log2FoldChange > 0)
down8 <- subset(resSig8, log2FoldChange < 0)

write.csv(up8, "data/results_up8.csv", row.names=TRUE)

write.csv(down8, "data/results_down8.csv", row.names=TRUE)

write.csv(resSig8, "data/results8.csv", row.names=TRUE)
