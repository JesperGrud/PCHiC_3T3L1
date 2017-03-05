```R
# Load the necessary packages
library(DESeq2)
library(e1071)

## Process super-enhancer data
# Import the data
Counts <- read.delim("Data/ChIPseq/Counts/Superenhancers_MED1.count")

# Setup to run DESeq2
rownames(Counts) <- Counts$PeakID

# Setup to run DESeq2
rownames(Counts) <- Counts$PeakID
Design <- data.frame(Condition = c("D0","h4","D2","D4","D7","D0","h4","D2","D4","D7"), Replicate = c(rep("a",5),rep("b",5)))
rownames(Design) <- colnames(Counts)[c(5:14)]

# Initialize DESeq2
DDS <- DESeqDataSetFromMatrix(countData = as.matrix(Counts[,c(5:14)]), colData = Design, design = ~Condition)

# Inject custom size factor (HOMER-style normalization: Total tags / 10 million)
sizeFactors(DDS)=c(1.2620968, 1.3888037, 1.3180341, 1.2340109, 1.2289277, 1.2801847, 1.4696396, 1.4138558, 1.6553577, 1.1341176)

# Estimate dispersions and extract normalized tag counts
DDS=estimateDispersions(DDS)
Normalized <- varianceStabilizingTransformation(DDS, blind=F)

# Perform WaldTests (all-vs-all) and FDR-adjust the result
colData(DDS)$Condition<-relevel(colData(DDS)$Condition, 'D0')
D0 <- nbinomWaldTest(DDS, betaPrior=F)
D0 <- apply(as.matrix(mcols(D0)[,grep("WaldPvalue_Condition", colnames(mcols(D0)))]),2,function(X) {p.adjust(X, method='BH')})

colData(DDS)$Condition<-relevel(colData(DDS)$Condition, 'h4')
H4 <- nbinomWaldTest(DDS, betaPrior=F)
H4 <- apply(as.matrix(mcols(H4)[,grep("WaldPvalue_Condition", colnames(mcols(H4)))]),2,function(X) {p.adjust(X, method='BH')})

colData(DDS)$Condition<-relevel(colData(DDS)$Condition, 'D2')
D2 <- nbinomWaldTest(DDS, betaPrior=F)
D2 <- apply(as.matrix(mcols(D2)[,grep("WaldPvalue_Condition", colnames(mcols(D2)))]),2,function(X) {p.adjust(X, method='BH')})

colData(DDS)$Condition<-relevel(colData(DDS)$Condition, 'D4')
D4 <- nbinomWaldTest(DDS, betaPrior=F)
D4 <- apply(as.matrix(mcols(D4)[,grep("WaldPvalue_Condition", colnames(mcols(D4)))]),2,function(X) {p.adjust(X, method='BH')})

colData(DDS)$Condition<-relevel(colData(DDS)$Condition, 'D7')
D7 <- nbinomWaldTest(DDS, betaPrior=F)
D7 <- apply(as.matrix(mcols(D7)[,grep("WaldPvalue_Condition", colnames(mcols(D7)))]),2,function(X) {p.adjust(X, method='BH')})

# Combine the FDR-adjusted p-values and normalized counts
Result=cbind(Counts[,c(1:4,15:19)],assay(Normalized),D0,H4,D2,D4,D7)

# Identify differentially occupied super-enhancers (FDR-adjusted p-values <= 0.05)
Result$Differential <- 0
Result[ apply(Result[,c(20:39)],1,FUN="min") <= 0.05,"Differential"] <- 1

# Calculate average occupancy in replicates
Result$Count_D0 <- rowMeans(2^Result[,c("MED1_D0_exp1","MED1_D0_exp2")])
Result$Count_H4 <- rowMeans(2^Result[,c("MED1_4h_exp1","MED1_4h_exp2")])
Result$Count_D2 <- rowMeans(2^Result[,c("MED1_D2_exp1","MED1_D2_exp2")])
Result$Count_D4 <- rowMeans(2^Result[,c("MED1_D4_exp1","MED1_D4_exp2")])
Result$Count_D7 <- rowMeans(2^Result[,c("MED1_D7_exp1","MED1_D7_exp2")])

# Partition data into 5 cluster using fuzzy clustering for the differential super-enhancers (on both replicates)
rownames(Result) <- Result$PeakID
ClusterData <- Result[ Result$Differential == 1 & (Result$D0 == 1 | Result$D2 == 1 | Result$H4 == 1 | Result$D4 == 1 | Result$D7 == 1), c("Count_D0","Count_H4","Count_D2","Count_D4","Count_D7")]
ClusterData <- log2(ClusterData)
ClusterData <- t(scale(t(ClusterData)))

# Set fuzzy-like parameter
Mest <- 2

# Perform the clustering
NoClusters <- 5
set.seed(30)
Clusters <- cmeans(ClusterData, NoClusters)
Centers <- Clusters$centers

## Filter based on the full dataset
# Calculate membership for each replicate and the mean within each cluster
Rep1 <- as.matrix(Result[ Result$Differential == 1 ,c("MED1_D0_exp1","MED1_4h_exp1","MED1_D2_exp1","MED1_D4_exp1","MED1_D7_exp1")])
Rep1 <- t(scale(t(Rep1)))
dm <- sapply(seq_len(nrow(Rep1)),function(i) apply(Centers, 1, function(v) sqrt(sum((Rep1[i, ]-v)^2))))
Rep1Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Rep1Membership <- as.data.frame(Rep1Membership)
rownames(Rep1Membership) <- Result[ Result$Differential == 1,c("PeakID")]

Rep2 <- as.matrix(Result[ Result$Differential == 1 ,c("MED1_D0_exp2","MED1_4h_exp2","MED1_D2_exp2","MED1_D4_exp2","MED1_D7_exp2")])
Rep2 <- t(scale(t(Rep2)))
dm <- sapply(seq_len(nrow(Rep2)),function(i) apply(Centers, 1, function(v) sqrt(sum((Rep2[i, ]-v)^2))))
Rep2Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Rep2Membership <- as.data.frame(Rep2Membership)
rownames(Rep2Membership) <- Result[ Result$Differential == 1,c("PeakID")]

All <- as.matrix(Result[ Result$Differential == 1 ,c("Count_D0","Count_H4","Count_D2","Count_D4","Count_D7")])
All <- log2(All)
All <- t(scale(t(All)))
dm <- sapply(seq_len(nrow(All)),function(i) apply(Centers, 1, function(v) sqrt(sum((All[i, ]-v)^2))))
Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Membership <- as.data.frame(Membership)
rownames(Membership) <- Result[ Result$Differential == 1,c("PeakID")]

# Define clusters based on membership using replicates
Rep1Membership$Max <- apply(Rep1Membership,1,FUN="max")
Rep1Membership$Include <- 0
Rep1Membership[ Rep1Membership$Max >= 0.3,"Include"] <- 1
Rep1Membership$Which <- 0
for (i in 1:nrow(Rep1Membership)) { Rep1Membership[i,"Which"] <- which.max(Rep1Membership[i,c(1:(NoClusters))]) }
Rep1Membership$Which <- Rep1Membership$Which * Rep1Membership$Include
Rep1Membership$PeakID <- rownames(Rep1Membership)

Rep2Membership$Max <- apply(Rep2Membership,1,FUN="max")
Rep2Membership$Include <- 0
Rep2Membership[ Rep2Membership$Max >= 0.3,"Include"] <- 1
Rep2Membership$Which <- 0
for (i in 1:nrow(Rep2Membership)) { Rep2Membership[i,"Which"] <- which.max(Rep2Membership[i,c(1:(NoClusters))]) }
Rep2Membership$Which <- Rep2Membership$Which * Rep2Membership$Include
Rep2Membership$PeakID <- rownames(Rep2Membership)

Membership$Max <- apply(Membership,1,FUN="max")
Membership$Include <- 0
Membership[ Membership$Max >= 0.6,"Include"] <- 1
Membership$Which <- 0
for (i in 1:nrow(Membership)) { Membership[i,"Which"] <- which.max(Membership[i,c(1:(NoClusters))]) }
Membership$Which <- Membership$Which * Membership$Include
Membership$PeakID <- rownames(Membership)

Clusters <- merge(Membership[,c("PeakID","Which","Max")],Rep1Membership[,c("PeakID","Which")], by="PeakID")
Clusters <- merge(Clusters, Rep2Membership[,c("PeakID","Which")], by="PeakID")
Clusters <- Clusters[ Clusters$Which.x != 0,]
Clusters <- Clusters[ Clusters$Which.x == Clusters$Which.y,]
Clusters <- Clusters[ Clusters$Which.x == Clusters$Which,]
colnames(Clusters)[2] <- "Cluster"
colnames(Clusters)[3] <- "Membership"
Result <- merge(Result[,c(1:45)], Clusters[,c(1,2,3)], by="PeakID", all.x=T)
Result[ is.na(Result$Cluster), "Cluster"] <- 0
Result[ is.na(Result$Membership), "Membership"] <- 0
Result[ Result$Differential == 0, "Cluster"] <- 0

# Rename clusters for intuitive cluster profiles
Result[ Result$Cluster == 4, "Cluster"] <- 7
Result[ Result$Cluster == 1, "Cluster"] <- 8
Result[ Result$Cluster == 3, "Cluster"] <- 9
Result[ Result$Cluster == 5, "Cluster"] <- 10
Result[ Result$Cluster == 2, "Cluster"] <- 11
Result[ Result$Cluster == 7, "Cluster"] <- 1
Result[ Result$Cluster == 8, "Cluster"] <- 2
Result[ Result$Cluster == 9, "Cluster"] <- 3
Result[ Result$Cluster == 10, "Cluster"] <- 4
Result[ Result$Cluster == 11, "Cluster"] <- 5

# Plot the cluster profiles
par(mfrow=c(2,3),pty='s')
color_vec=c("#0000CD30", "#0000CD30","#00640030","#00640030","#CD000030")
rownames(Result) <- Result$PeakID
ClusterData <- Result[ Result$Differential == 1, c("Count_D0","Count_H4","Count_D2","Count_D4","Count_D7")]
ClusterData <- log2(ClusterData)
ClusterData <- t(scale(t(ClusterData)))

for (i in 1:NoClusters) {
	tmp <- ClusterData[ rownames(ClusterData) %in% Result[ Result$Cluster == i,"PeakID"],]
	plot(0,0,pch=' ', xlim=c(1,5), ylim=c(-2.6,2.6), xlab="", ylab="", yaxt="n", xaxt="n", main=paste("Cluster", i,sep = " "))
	mtext("Scaled binding", side=2, line=2.5, cex=1)
	mtext(c("Timepoint during differentiation"), side=1, line=2, cex=1)
	for(j in 1:dim(tmp)[1]){ lines(1:5, tmp[j,c(1:5)], col=color_vec[i]) }
	lines(1:5, apply(tmp[,c(1:5)],2, mean), lwd=2)
	axis(2, at=c(-2,-1,0,1,2), lab=c(-2,-1,0,1,2), cex.axis=1.3, las=2)
	axis(1, at=c(1:5), lab=c("D0", "4h", "D2","D4","D7"))
	}
```

[Back to start](../README.md)<br>
[Back to overview of Figure S3](../Links/FigureS3.md)
