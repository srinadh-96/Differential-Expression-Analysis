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

data <- read.delim("data/3_500_500.tsv", header=T, row.names=1)
meta <-  read.delim("meta/Book1.txt", header=T, row.names=1)

data
View(meta)
class(data)
class(meta)



all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData=meta, 
                              design=~condition)
dds

# Generate the normalised counts
dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)



### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
### Plot PCA
plotPCA(rld, intgroup="sampletype")


# hierarchical clustering

### Extract the rlog matrix from the object

rld_mat <- assay(rld)  ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

pheatmap(rld_cor)

# Dispersion estimates
dds <- DESeq(dds)
plotDispEsts(dds)


dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData=meta, 
                              design=condition)

dds <- DESeq(dds)
plotDispEsts(dds)


## DE analysis

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~condition)
# DE analysis itself
dds <- DESeq(dds)


# Wald test

BiocManager::install("apeglm")

contrast = c("condition", "condition1", "condition2")

res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.3)

res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken,
                         type="ashr")

install.packages('ashr')

## MA plot

plotMA(res_tableOE_unshrunken, ylim=c(-2,2))

plotMA(res_tableOE, ylim=c(-2,2))


class(res_tableOE)

mcols(res_tableOE, use.names=T)

res_tableOE %>% data.frame() %>% View()

summary(res_tableOE)


##set thresholds
padj.cutoff <- 0.3
lfc.cutoff <- 0.58

res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

sigOE <- res_tableOE_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sigOE  

summary(sigOE)

res <- results(dds, contrast = contrast_oe, alpha = 0.1, lfcThreshold = 0.2)
summary(res)
View(res)

res <- results(dds, contrast = contrast, alpha = 0.05, lfcThreshold = 0)
summary(res)

res_tableOE_tb <- res_tableOE_tb %>%
  mutate(threshold_OE = padj < 0.1 & abs(log2FoldChange) >= 0.2)

ggplot(res_tableOE_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("control") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))


res_tableOE_tb <- res_tableOE_tb %>% arrange(padj) %>% mutate(genelabels = "")

res_tableOE_tb$genelabels[1:10] <- res_tableOE_tb$gene[1:10]


View(res_tableOE_tb)

summary(res_tableOE_tb)

ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(colour = threshold_OE)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("control") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))


res

resSig <- res[ res$padj < 0.1, ]

summary(res)
res1 = res[ which(res$log2FoldChange > 0 & res$padj < 0.1),]
res <- res[!is.na(res$padj),]
up <- subset(resSig, log2FoldChange > 0)
down <- subset(resSig, log2FoldChange < 0)
resSig
write.csv(resSig, "data/results.csv", row.names=TRUE)
up
write.csv(up, "data/results_up.csv", row.names=TRUE)
down
write.csv(down, "data/results_down.csv", row.names=TRUE)
resSig %>% data.frame() %>% View()


write.csv(resSig, "data/results.csv", row.names=TRUE)

View(dds)

res_data <- as.data.frame(resSig)
res_data

res_data %>% mutate(diff_expr =
                     case_when(log2FoldChange < 0 ~ 0, 
                               log2FoldChange > 0 ~ 1)
)

res_data<-transform(res_data,upregulation=ifelse(log2FoldChange > 0,1 ,0))

res_data<-transform(res_data,downregulation=ifelse(log2FoldChange <= 0,1 ,0))
res_data<-transform(res_data,differential.expression=ifelse(upregulation | downregulation == 1 ,1 ,0))



df = subset(res_data, select = c(upregulation, downregulation, differential.expression))
write.csv(df, "data/results.csv", row.names=TRUE)
GT <- read.delim("meta/3_500_500_meta.tsv", header=T, row.names=1) 
GT1 = subset(GT, select = -c(upregulation, downregulation, differential.expression) )
GT1
write.csv(GT1, "data/3_500_500.csv", row.names=TRUE)
