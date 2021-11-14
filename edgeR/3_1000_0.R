library(edgeR)

# Load data 
data <- read.delim("data/3_1000_0.tsv", header=T, row.names=1)
head(data)

# Create DGE list object 
group <- factor(c('condition1','condition1', 'condition1', 'condition2','condition2','condition2'))
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

fit <- glmQLFit(DGE, design)
res <- glmQLFTest(fit)
pVals <- res$table[,4]
names(pVals) <- rownames(res$table)

pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)

# Dispersion 
DGE <- estimateGLMCommonDisp(DGE, design=design)
DGE <- estimateGLMTrendedDisp(DGE, design=design)
DGE <- estimateGLMTagwiseDisp(DGE, design=design)
plotBCV(DGE)


# De Analysis 
et12 <- exactTest(DGE, pair=c(1,2)) # compare 1 and 4
topTags <- topTags(et12, n=500) 
topTags

de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.1)
summary(de1)

# Print DE genes, upregulated and downregulated
write.csv(DE.up,"/Users/kariprivat/Documents/results_3_1000_0_up.csv", row.names = TRUE)
write.csv(DE.dowm,"/Users/kariprivat/Documents/results_3_1000_0_down.csv", row.names = TRUE)
write.csv(DE,"/Users/kariprivat/Documents/results_3_1000_.csv", row.names = TRUE)

