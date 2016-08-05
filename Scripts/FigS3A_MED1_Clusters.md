```R
## Load the necessary packages
library(DESeq2)

## Process enhancer data
# Import the data
Counts <- read.delim("Data/ChIPseq/Counts/MED1.count")

# Keep only distal sites (>= 2kb away from closest TSS)
Counts <- Counts[ abs(Counts$DistanceToNearestGene) >= 2000,]

# Setup to run DESeq2
rownames(Counts) <- Counts$PeakID
Design <- data.frame(Condition = c("D0","h4","D2","D4","D7","D0","h4","D2","D4","D7"), Replicate = c(rep("a",5),rep("b",5)))
rownames(Design) <- colnames(Counts)[c(6:15)]

# Initialize DESeq2
DDS <- DESeqDataSetFromMatrix(countData = as.matrix(Counts[,c(6:15)]), colData = Design, design = ~Condition)

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
Result=cbind(Counts[,c(1:4)],assay(Normalized),D0,H4,D2,D4,D7)

# Identify differentially occupied enhancers (FDR-adjusted p-values <= 0.05)
Result$Differential <- 0
Result[ apply(Result[,c(15:34)],1,FUN="min") <= 0.05,"Differential"] <- 1

# Calculate average occupancy in replicates
Result$Count_D0 <- rowMeans(Result[,c("MED1_D0_exp1","MED1_D0_exp2")])
Result$Count_H4 <- rowMeans(Result[,c("MED1_4h_exp1","MED1_4h_exp2")])
Result$Count_D2 <- rowMeans(Result[,c("MED1_D2_exp1","MED1_D2_exp2")])
Result$Count_D4 <- rowMeans(Result[,c("MED1_D4_exp1","MED1_D4_exp2")])
Result$Count_D7 <- rowMeans(Result[,c("MED1_D7_exp1","MED1_D7_exp2")])

# Center the averaged counts
Result[,c(36:40)] <- t(scale(t(Result[,c(36:40)]), center=T, scale=F))

# Partition data into 5 cluster using k-means for the differential super-enhancers (on both replicates)
rownames(Result) <- Result$PeakID
set.seed(13)
Clusters <- kmeans(t(scale(t(as.matrix(Result[ Result$Differential == 1,c(5:14)])), center=T, scale=F)), 7, iter.max=500)
Clustering <- data.frame(Cluster = Clusters$cluster)
Clustering$PeakID <- rownames(Clustering)

# Rename clusters for intuitive cluster profiles
Clustering[ Clustering$Cluster == 1, "Cluster"] <- 8
Clustering[ Clustering$Cluster == 2, "Cluster"] <- 9
Clustering[ Clustering$Cluster == 3, "Cluster"] <- 10
Clustering[ Clustering$Cluster == 4, "Cluster"] <- 11
Clustering[ Clustering$Cluster == 6, "Cluster"] <- 12
Clustering[ Clustering$Cluster == 7, "Cluster"] <- 13
Clustering[ Clustering$Cluster == 8, "Cluster"] <- 6
Clustering[ Clustering$Cluster == 9, "Cluster"] <- 3
Clustering[ Clustering$Cluster == 10, "Cluster"] <- 7
Clustering[ Clustering$Cluster == 11, "Cluster"] <- 1
Clustering[ Clustering$Cluster == 12, "Cluster"] <- 4
Clustering[ Clustering$Cluster == 13, "Cluster"] <- 2

# Plot the cluster profiles
par(mfcol=c(2,4))
color_vec=c("#0000CD30","#0000CD30","#0000CD30", "#00640030","#00640030","#00640030","#CD000030")
for (i in 1:7) {
	tmp <- Result[ Result$Differential == 1 ,c(1,36:38)]
	tmp <- tmp[ tmp$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"], ]
	plot(0,0,pch=' ', xlim=c(1,3), ylim=c(-2.6,2.6), xlab="", ylab="", yaxt="n", xaxt="n", main=paste("Cluster", i,sep = " "))
	mtext("Relative binding (log2)", side=2, line=2.5, cex=1)
	mtext(c("Timepoint during differentiation"), side=1, line=2, cex=1)
	for(j in 1:dim(tmp)[1]){ lines(1:3, tmp[j,c(2:4)], col=color_vec[i]) }
	lines(1:3, apply(tmp[,c(2:4)],2, median), lwd=2)
	axis(2, at=c(-2,-1,0,1,2), lab=c(-2,-1,0,1,2), cex.axis=1.3, las=2)
	axis(1, at=c(1:3), lab=c("D0", "4h", "D2"))
	}
```

[Back to start](../README.md)<br>
[Back to overview of Figure S3](../Links/FigureS3.md)
