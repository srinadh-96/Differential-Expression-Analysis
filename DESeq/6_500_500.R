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


data <- read.delim("data/6_500_500.tsv", header=T, row.names=1)
meta <- read.delim("meta/6_metadata.tsv", header=T, row.names=1)

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))


dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData=meta, 
                              design= ~ condition)

dds

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)


dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
dds <- DESeq(dds)

contrast_oe <- c("condition", "condition2", "condition1")


# pvalue = 0.05
res <- results(dds, contrast = contrast_oe, alpha = 0.05, lfcThreshold = 0)
summary(res)

res <- res[!is.na(res$padj),]
resSig <- res[ res$padj < 0.05, ]

up <- subset(resSig, log2FoldChange > 0)
down <- subset(resSig, log2FoldChange < 0)

write.csv(down, "data/pvalue5/results_down_6_500_500.csv", row.names=TRUE)

write.csv(up, "data/pvalue5/results_up_6_500_500.csv", row.names=TRUE)

write.csv(resSig, "data/pvalue5/results_6_500_500.csv", row.names=TRUE)


#pvalue = 0.1

res2 <- results(dds, contrast = contrast_oe, alpha = 0.1, lfcThreshold = 0)
summary(res2)

res2 <- res2[!is.na(res2$padj),]
resSig2 <- res2[ res2$padj < 0.1, ]

up2 <- subset(resSig2, log2FoldChange > 0)
down2 <- subset(resSig2, log2FoldChange < 0)

write.csv(down2, "data/pvalue10/results_down_6_500_500.csv", row.names=TRUE)

write.csv(up2, "data/pvalue10/results_up_6_500_500.csv", row.names=TRUE)

write.csv(resSig2, "data/pvalue10/results_6_500_500.csv", row.names=TRUE)


#pvalue = 0.2

res3 <- results(dds, contrast = contrast_oe, alpha = 0.2, lfcThreshold = 0)
summary(res3)

res3 <- res3[!is.na(res3$padj),]
resSig3 <- res3[ res3$padj < 0.2, ]

up3 <- subset(resSig3, log2FoldChange > 0)
down3 <- subset(resSig3, log2FoldChange < 0)

write.csv(down3, "data/pvalue20/results_down_6_500_500.csv", row.names=TRUE)

write.csv(up3, "data/pvalue20/results_up_6_500_500.csv", row.names=TRUE)

write.csv(resSig3, "data/pvalue20/results_6_500_500.csv", row.names=TRUE)

