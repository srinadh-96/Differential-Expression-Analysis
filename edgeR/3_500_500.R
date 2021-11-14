library(edgeR)
install.packages("writexl")
library("writexl")

# Load data 
data <- read.delim("data/3_500_500.tsv", header=T, row.names=1)
head(data)
groundTruth <- read.table("ground_truth/3_500_500_meta.tsv")


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

# Dispersion 
DGE <- estimateGLMCommonDisp(DGE, design=design)
DGE <- estimateGLMTrendedDisp(DGE, design=design)
DGE <- estimateGLMTagwiseDisp(DGE, design=design)
plotBCV(DGE)


# De Analysis 
et12 <- exactTest(DGE, pair=c(1,2)) # compare group 1 and 2
topTags <- topTags(et12, n=1000) 
topTags

de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.1)
summary(de1)

# Plot DE genes 
de1tags12 <- rownames(DGE)[as.logical(de1)]
plotSmear(et12, de.tags=de1tags12)
abline(h=c(-2,2), col = "blue")

out <- topTags(et12, n = Inf)$table
DE.up <- out[out$logFC>0.0, ]
DE.down <- out[out$logFC<0.0, ]
DE <- out[]

write.csv(DE.up,"/Users/kariprivat/Documents/results_3_500_500_up.csv", row.names = TRUE)
write.csv(DE.dowm,"/Users/kariprivat/Documents/results_3_500_500_down.csv", row.names = TRUE)
write.csv(DE,"/Users/kariprivat/Documents/results_3_500_500.csv", row.names = TRUE)

sigDE <- names(de1)[de1 < 0.05]
length(sigDE) 

