library(edgeR)

# Load data 
data <- read.delim("data/9_1000_0.tsv", header=T, row.names=1)
head(data)

# Create DGE list object 
group <- factor(c('condition1','condition1', 'condition1', 'condition1','condition1', 'condition1', 'condition1','condition1', 'condition1', 'condition2','condition2','condition2', 'condition2','condition2','condition2', 'condition2','condition2','condition2'))
DGE <- DGEList(counts=data,group=group)

# Total counts per sample 
apply(DGE$counts, 2, sum)

# Normalise
DGE <- calcNormFactors(DGE, method="TMM")

# Plot samples
plotMDS(DGE)

# Create model 
design <- model.matrix(~0+group, data=DGE$samples)
colnames(design) <- levels(DGE$samples$group)
design

# Dispersion 
DGE <- estimateGLMCommonDisp(DGE, design=design)
DGE <- estimateGLMTrendedDisp(DGE, design=design)
DGE <- estimateGLMTagwiseDisp(DGE, design=design)
plotBCV(DGE)


# De Analysis 
et12 <- exactTest(DGE, pair=c(1,2)) # compare 1 and 2
topTags <- topTags(et12, n=500) 
topTags

de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.1)
summary(de1)

# Plot DE genes 
de1tags12 <- rownames(DGE)[as.logical(de1)]
plotSmear(et12, de.tags=de1tags12)
abline(h=c(-2,2), col = "blue")


write.csv(de1,"/Users/kariprivat/Documents/results_9_1000_0.csv", row.names = TRUE)
write.csv(et12,"/Users/kariprivat/Documents/exact_results_9_1000_0.csv", row.names = TRUE)

