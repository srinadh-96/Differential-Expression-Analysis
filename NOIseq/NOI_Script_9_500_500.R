library(NOISeq)

data_1 <- read.delim("data/9_500_500.tsv", header = T,row.names = 1)
head(data_1)
data_2<- read.delim("data/Book3.tsv", header = T,row.names = 1)
head(data_2)
mycounts <- data_1
myfactors <- data_2
mydata <-  readData(data = mycounts, factors = myfactors)


str(mydata)
head(assayData(mydata)$exprs)
head(pData(mydata))
head(featureData(mydata)@data)

set.seed(123)
mycounts2 = mycounts
mycounts2[, 1:4] = mycounts2[, 1:4] + runif(nrow(mycounts2) * 4, 3, 5)
myfactors = data.frame(myfactors, batch = c(rep(1, 3), rep(2, 3)))
mydata2 = readData(mycounts2, factors = myfactors)


myPCA = dat(mydata2, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA, factor = "sampletype")
explo.plot(myPCA, factor = "condition")

QCreport(mydata, samples = NULL, factor = "condition", norm = false)

myTMM = tmm(assayData(mydata)$exprs, long = 1000, lc = 0)
head(myTMM[, 1:4])


myfilt = filtered.data(mycounts, factor = myfactors$sampletype, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")

mydata2corr1 = ARSyNseq(mydata2, factor = "condition", batch = TRUE, norm = "rpkm",logtransf = FALSE)
myPCA = dat(mydata2corr1, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA, factor = "sampletype")
explo.plot(myPCA, factor = "condition")


mydata2corr2 = ARSyNseq(mydata2, factor = "sampletype", batch = FALSE, norm = "rpkm",logtransf = FALSE)
myPCA = dat(mydata2corr2, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA, factor = "sampletype")
explo.plot(myPCA, factor = "condition")



mynoiseq.tmm = noiseq(mydata, k = 0, norm = "tmm", factor = "sampletype",conditions = c("control", "treat"), lc = 0, replicates = "technical")
head(mynoiseq.tmm@results[[1]])

mynoiseq.deg = degenes(mynoiseq.tmm, q = 0.7, M = NULL)
mynoiseq.deg
write.csv(mynoiseq.deg,"Data/r_results/results_9_500_0.7.csv", row.names = TRUE)


mynoiseq.deg1 = degenes(mynoiseq.tmm, q = 0.7, M = "up")
write.csv(mynoiseq.deg1,"Data/r_results/results_9_500_0.7_up.csv", row.names = TRUE)


mynoiseq.deg2 = degenes(mynoiseq.tmm, q = 0.7, M = "down")
write.csv(mynoiseq.deg2,"Data/r_results/results_9_500_0.7_down.csv", row.names = TRUE)


