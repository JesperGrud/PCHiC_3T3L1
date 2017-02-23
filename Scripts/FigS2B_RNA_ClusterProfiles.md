```R
# Load libraries
library(DESeq2)
library(org.Mm.eg.db)
library(e1071)

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

# Import transcript length information and merge it
Lengths <- read.delim("Data/RNAseq/Annotation/Gene.lengths", header=T)
RNA <- merge(RNA, Lengths, by="RefSeq")

# Normalize RNA-seq counts to txLength
RNA$D0 <- RNA$D0 / (RNA$txLength/1000)
RNA$H4 <- RNA$H4 / (RNA$txLength/1000)
RNA$D1 <- RNA$D1 / (RNA$txLength/1000)
RNA$D2 <- RNA$D2 / (RNA$txLength/1000)
RNA$D7 <- RNA$D7 / (RNA$txLength/1000)

# Calculate maximal expression
RNA$Max <- apply(RNA[,c("D0","H4","D1","D2","D7")],1,FUN="max")

# Calculate log2 fold changes
RNA$logFC_H4 <- RNA$H4 / RNA$D0
RNA$logFC_H4 <- log2(RNA$logFC_H4)
RNA$logFC_D1 <- RNA$D1 / RNA$D0
RNA$logFC_D1 <- log2(RNA$logFC_D1)
RNA$logFC_D2 <- RNA$D2 / RNA$D0
RNA$logFC_D2 <- log2(RNA$logFC_D2)
RNA$logFC_D7 <- RNA$D7 / RNA$D0
RNA$logFC_D7 <- log2(RNA$logFC_D7)

# Define genes to cluster
RNA$Include <- 0
RNA[ (abs(RNA$logFC_H4) > 0 & RNA$Padj_Four <= 0.01) | (abs(RNA$logFC_D2) > 0 & RNA$Padj_Two <= 0.01), "Include"] <- 1

# Do clustering only for genes significant between D0 and either H4 or D2
rownames(RNA) <- RNA$RefSeq
ClusterData <- as.data.frame(RNA[ ,c("D0","H4","D2")])
ClusterData <- log2(ClusterData)
ClusterData <- as.data.frame(t(scale(t(ClusterData))))
ClusterData <- ClusterData[ rownames(ClusterData) %in% RNA[ RNA$Include == 1,"RefSeq"],]

# Set fuzzy-like parameter
Mest <- 2

# Do the clustering
NoClusters <- 8
set.seed(15)
Clusters <- cmeans(ClusterData, NoClusters, m = Mest, method = "cmeans")
Centers <- Clusters$centers

# Calculate membership for each replicate and the mean within each cluster
Rep1 <- as.matrix(RNA[ RNA$Include == 1 ,c("D0_E1","H4_E1","D2_E1")])
Rep1 <- t(scale(t(Rep1)))
dm <- sapply(seq_len(nrow(Rep1)),function(i) apply(Centers, 1, function(v) sqrt(sum((Rep1[i, ]-v)^2))))
Rep1Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Rep1Membership <- as.data.frame(Rep1Membership)
rownames(Rep1Membership) <- RNA[ RNA$Include == 1,c("RefSeq")]

Rep2 <- as.matrix(RNA[ RNA$Include == 1,c("D0_E2","H4_E2","D2_E2")])
Rep2 <- t(scale(t(Rep2)))
dm <- sapply(seq_len(nrow(Rep2)),function(i) apply(Centers, 1, function(v) sqrt(sum((Rep2[i, ]-v)^2))))
Rep2Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Rep2Membership <- as.data.frame(Rep2Membership)
rownames(Rep2Membership) <- RNA[ RNA$Include == 1,c("RefSeq")]

All <- as.matrix(RNA[ RNA$Include == 1,c("D0","H4","D2")])
All <- log2(All)
All <- t(scale(t(All), center = T, scale = F))
dm <- sapply(seq_len(nrow(All)),function(i) apply(Centers, 1, function(v) sqrt(sum((All[i, ]-v)^2))))
Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Membership <- as.data.frame(Membership)
rownames(Membership) <- RNA[ RNA$Include == 1,c("RefSeq")]

# Define clusters based on membership using replicates
Rep1Membership$Max <- apply(Rep1Membership,1,FUN="max")
Rep1Membership$Include <- 0
Rep1Membership[ Rep1Membership$Max >= 0.2,"Include"] <- 1
Rep1Membership$Which <- 0
for (i in 1:nrow(Rep1Membership)) { Rep1Membership[i,"Which"] <- which.max(Rep1Membership[i,c(1:(NoClusters))]) }
Rep1Membership$Which <- Rep1Membership$Which * Rep1Membership$Include
Rep1Membership$RefSeq <- rownames(Rep1Membership)

Rep2Membership$Max <- apply(Rep2Membership,1,FUN="max")
Rep2Membership$Include <- 0
Rep2Membership[ Rep2Membership$Max >= 0.2,"Include"] <- 1
Rep2Membership$Which <- 0
for (i in 1:nrow(Rep2Membership)) { Rep2Membership[i,"Which"] <- which.max(Rep2Membership[i,c(1:(NoClusters))]) }
Rep2Membership$Which <- Rep2Membership$Which * Rep2Membership$Include
Rep2Membership$RefSeq <- rownames(Rep2Membership)

Membership$Max <- apply(Membership,1,FUN="max")
Membership$Include <- 0
Membership[ Membership$Max >= 0.3,"Include"] <- 1
Membership$Which <- 0
for (i in 1:nrow(Membership)) { Membership[i,"Which"] <- which.max(Membership[i,c(1:(NoClusters))]) }
Membership$Which <- Membership$Which * Membership$Include
Membership$RefSeq <- rownames(Membership)

Clusters <- merge(Membership[,c("Max","RefSeq","Which")], Rep1Membership[,c("RefSeq","Which")], by="RefSeq")
Clusters <- merge(Clusters, Rep2Membership[,c("RefSeq","Which")], by="RefSeq")
Clusters <- Clusters[ Clusters$Which.x != 0,]
Clusters <- Clusters[ Clusters$Which.x == Clusters$Which.y,]
Clusters <- Clusters[ Clusters$Which.x == Clusters$Which,]
colnames(Clusters)[3] <- "Cluster"
colnames(Clusters)[2] <- "Membership"
RNA <- merge(RNA[,c(1:29)], Clusters[,c(1,2,3)], by="RefSeq", all.x=T)
RNA[ is.na(RNA$Cluster), "Cluster"] <- 0
RNA[ is.na(RNA$Membership), "Membership"] <- 0

# Filter the results by expression
RNA[ RNA$Include == 0 & RNA$Cluster < (NoClusters), "Cluster"] <- 0

# Remove too small clusters
RNA[ RNA$Cluster == 5, "Cluster"] <- 0

# Rename clusters for intuitive grouping during plots
RNA[ RNA$Cluster == 4, "Cluster"] <- 10
RNA[ RNA$Cluster == 8, "Cluster"] <- 11
RNA[ RNA$Cluster == 3, "Cluster"] <- 12
RNA[ RNA$Cluster == 7, "Cluster"] <- 13
RNA[ RNA$Cluster == 2, "Cluster"] <- 14
RNA[ RNA$Cluster == 1, "Cluster"] <- 15
RNA[ RNA$Cluster == 6, "Cluster"] <- 16
RNA[ RNA$Cluster == 10, "Cluster"] <- 1
RNA[ RNA$Cluster == 11, "Cluster"] <- 2
RNA[ RNA$Cluster == 12, "Cluster"] <- 3
RNA[ RNA$Cluster == 13, "Cluster"] <- 4
RNA[ RNA$Cluster == 14, "Cluster"] <- 5
RNA[ RNA$Cluster == 15, "Cluster"] <- 6
RNA[ RNA$Cluster == 16, "Cluster"] <- 7

# Plot the line plots for all clusters
par(mfcol=c(1,4))
color_vec=c("#0000CD30", "#00640030","#CD000030")
for (i in c(4:7)) {
tmp <- RNA[ RNA$Cluster == i, c("D0","H4","D2")]
tmp <- log2(tmp)
tmp <- t(scale(t(tmp)))
plot(0,0,pch=' ', xlim=c(1,3), ylim=c(-2.5,2.5), xlab="", ylab="", yaxt="n", xaxt="n", main=paste("Cluster",i,sep=" "))
mtext("Scaled expression", side=2, line=2.5, cex=1)
mtext(c("Timepoint during differentiation"), side=1, line=2, cex=1)
for(j in 1:dim(tmp)[1]){ lines(1:3, tmp[j,], col=color_vec[(i)]) }
lines(1:3, apply(tmp,2, mean), lwd=2)
axis(2, at=c(seq(-10,10,by=1)), lab=c(seq(-10,10,by=1)), cex.axis=1.3, las=2)
axis(1, at=c(1:3), lab=c("D0", "4h", "D2"))
}
```

[Back to start](../README.md)<br>
[Back to overview of Figure S2](../Links/FigureS2.md)
