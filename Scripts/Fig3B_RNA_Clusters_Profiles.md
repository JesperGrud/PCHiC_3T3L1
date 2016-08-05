```R
## Load libraries
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

# Rename clusters for intuitive grouping during plots
RNA[ RNA$Cluster == 2, "Cluster"] <- 8
RNA[ RNA$Cluster == 4, "Cluster"] <- 9
RNA[ RNA$Cluster == 6, "Cluster"] <- 10
RNA[ RNA$Cluster == 7, "Cluster"] <- 11
RNA[ RNA$Cluster == 8, "Cluster"] <- 4
RNA[ RNA$Cluster == 9, "Cluster"] <- 2
RNA[ RNA$Cluster == 10, "Cluster"] <- 7
RNA[ RNA$Cluster == 11, "Cluster"] <- 6

# Plot the line plots for all clusters
par(mfcol=c(3,3))
color_vec=c("#90909030","#0000CD30", "#0000CD30","#CD000030","#00640030","#CD000030","#00640030","#CD000030")
for (i in 0:7) {
tmp <- RNA[ RNA$Cluster == i, c(6,7,9,11,12,14)]
tmp <- (tmp[,1:(ncol(tmp)/2)]+tmp[,((ncol(tmp)/2)+1):ncol(tmp)])/2
tmp <- as.matrix(t(scale(t(tmp),center=T, scale=F)))
plot(0,0,pch=' ', xlim=c(1,3), ylim=c(min(tmp),max(tmp)), xlab="", ylab="", yaxt="n", xaxt="n", main=paste("Cluster",i,sep=" "))
mtext("Relative expression (log2)", side=2, line=2.5, cex=1)
mtext(c("Timepoint during differentiation"), side=1, line=2, cex=1)
for(j in 1:dim(tmp)[1]){ lines(1:3, tmp[j,], col=color_vec[(i+1)]) }
lines(1:3, apply(tmp,2, median), lwd=2)
axis(2, at=c(seq(-10,10,by=1)), lab=c(seq(-10,10,by=1)), cex.axis=1.3, las=2)
axis(1, at=c(1:3), lab=c("D0", "4h", "D2"))
}
```

[Back to start](../README.md)<br>
[Back to overview of Figure 3](../Links/Figure3.md)