```R
## Load the necessary packages
library(DESeq2)
library(org.Mm.eg.db)

## Process RNA-seq data
# Import RNA-seq data
RNA <- read.delim("Data/RNAseq/RNA.txt", header=T)

# Perform differential expression analysis using DESeq2
MetaData <- data.frame(Condition = c("D0","H4","D1","D2","D7","D0","H4","D1","D2","D7"), Replicate = c(rep("a",5), rep("b",5)))
rownames(MetaData) <- colnames(RNA)[c(6:15)]
rownames(RNA) <- RNA$RefSeq
DGE <- DESeqDataSetFromMatrix(as.matrix(RNA[,c(6:15)]), colData = MetaData, design = ~ Condition)
DGE <- estimateSizeFactors(DGE)
DGE <- estimateDispersions(DGE)
DGE <- nbinomWaldTest(DGE, betaPrior=F)

# Extract p-values and do correction for multiple testing
Pvalues <- as.matrix(mcols(DGE)[,grep("WaldPvalue_Condition", colnames(mcols(DGE)))])
Padj <- as.data.frame(apply(as.matrix(Pvalues),2,function(X) {p.adjust(X, method='BH')}))
colnames(Padj) <- c("Padj_One","Padj_Two","Padj_Seven","Padj_Four")
Padj$RefSeq <- names(DGE)

# Get variance stabilized data 
VSD <- as.data.frame(assay(varianceStabilizingTransformation(DGE, blind=F)))
VSD$RefSeq <- rownames(VSD)

# Merge all of the data and remove NAs (All zero counts)
RNA <- merge(RNA[,c(1:5)], VSD, by="RefSeq", sort=F)
RNA <- merge(RNA, Padj, by="RefSeq", sort=F)
RNA <- na.omit(RNA)

# Calculate mean expression across replicates
RNA$D0 <- rowMeans(2^RNA[,c("D0_E1","D0_E2")])
RNA$H4 <- rowMeans(2^RNA[,c("H4_E1","H4_E2")])
RNA$D1 <- rowMeans(2^RNA[,c("D1_E1","D1_E2")])
RNA$D2 <- rowMeans(2^RNA[,c("D2_E1","D2_E2")])
RNA$D7 <- rowMeans(2^RNA[,c("D7_E1","D7_E2")])

# Do clustering only for genes significant between D0 and either H4 or D2
rownames(RNA) <- RNA$RefSeq
ClusterData <- as.matrix(RNA[ RNA$Padj_Two < 0.01 | RNA$Padj_Four < 0.01,c(6,7,9,11,12,14)])
set.seed(13)
Clusters <- kmeans(t(scale(t(ClusterData), center=T, scale=F)), 7, iter.max=500)
Clusters <- data.frame(Cluster = Clusters$cluster)
Clusters$RefSeq <- rownames(Clusters)
RNA <- merge(RNA, Clusters, by="RefSeq", all.x=T)
RNA[ is.na(RNA$Cluster), "Cluster"] <- 0

# Import transcript length information and merge it
Lengths <- read.delim("Data/RNAseq/Annotation/Gene.lengths", header=T)
RNA <- merge(RNA, Lengths, by="RefSeq")

# Normalize RNA-seq counts to txLength
RNA$D0 <- RNA$D0 / (RNA$txLength/1000)
RNA$H4 <- RNA$H4 / (RNA$txLength/1000)
RNA$D2 <- RNA$D2 / (RNA$txLength/1000)

# Extract mapping to official symbols
Convert <- suppressMessages(select(org.Mm.eg.db, as.character(RNA$RefSeq), columns = c("SYMBOL","REFSEQ"), "REFSEQ"))
colnames(Convert) <- c("RefSeq","Symbol")
RNA <- merge(RNA, Convert, by="RefSeq")

## Plot the expression of the examples
par(mfcol=c(1,3))
plot(seq(1,3,by=1),as.numeric(RNA[ RNA$Symbol %in% "Pparg",c("D0","H4","D2")]), ylim=c(0,4000), las=1, type="o", col="green3", lwd=2, ylab="Normalized gene expression", xlab="", xaxt="n", pch=18, cex=2, main="Pparg")
axis(1, at=c(1,2,3), labels=c("D0","4h","D2"))
plot(seq(1,3,by=1),as.numeric(RNA[ RNA$Symbol %in% "Hmga2",c("D0","H4","D2")]), ylim=c(0,2000), las=1, type="o", col="blue3", lwd=2, ylab="Normalized gene expression", xlab="", xaxt="n", pch=18, cex=2, main="Hmga2")
axis(1, at=c(1,2,3), labels=c("D0","4h","D2"))
plot(seq(1,3,by=1),as.numeric(RNA[ RNA$Symbol %in% "Maged2",c("D0","H4","D2")]), ylim=c(0,7000), las=1, type="o", col="red3", lwd=2, ylab="Normalized gene expression", xlab="", xaxt="n", pch=18, cex=2, main="Maged2")
axis(1, at=c(1,2,3), labels=c("D0","4h","D2"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 4](../Links/Figure4.md)