```R
## Load the necessary packages
library(GenomicRanges)
library(DESeq2)

## Process super-enhancer data
# Import the data
Counts <- read.delim("Data/ChIPseq/Counts/Superenhancers_MED1.count")

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
Result$Count_D0 <- rowMeans(Result[,c("MED1_D0_exp1","MED1_D0_exp2")])
Result$Count_H4 <- rowMeans(Result[,c("MED1_4h_exp1","MED1_4h_exp2")])
Result$Count_D2 <- rowMeans(Result[,c("MED1_D2_exp1","MED1_D2_exp2")])
Result$Count_D4 <- rowMeans(Result[,c("MED1_D4_exp1","MED1_D4_exp2")])
Result$Count_D7 <- rowMeans(Result[,c("MED1_D7_exp1","MED1_D7_exp2")])

# Center the averaged counts
Result[,c(41:45)] <- t(scale(t(Result[,c(41:45)]), center=T, scale=F))

# Partition data into 5 cluster using k-means for the differential super-enhancers (on both replicates)
rownames(Result) <- Result$PeakID
set.seed(13)
Clusters <- kmeans(t(scale(t(as.matrix(Result[ Result$Differential == 1,c(10:19)])), center=T, scale=F)), 5, iter.max=500)
Clustering <- data.frame(Cluster = Clusters$cluster)
Clustering$PeakID <- rownames(Clustering)

# Rename clusters for intuitive cluster profiles
Clustering[ Clustering$Cluster == 1, "Cluster"] <- 6
Clustering[ Clustering$Cluster == 3, "Cluster"] <- 7
Clustering[ Clustering$Cluster == 7, "Cluster"] <- 1
Clustering[ Clustering$Cluster == 6, "Cluster"] <- 3

## Process SMC1 binding data
# Setup to run DESeq2
rownames(Design) <- colnames(Counts)[c(20:29)]

# Initialize DESeq2
DDS <- DESeqDataSetFromMatrix(countData = as.matrix(Counts[,c(20:29)]), colData = Design, design = ~Condition)

# Inject custom size factor (HOMER-style normalization: Total tags / 10 million)
sizeFactors(DDS)=c(1.3281913, 1.3934452, 1.1364366, 1.2030321, 1.2917667, 1.3695474, 2.4996968, 1.2890919, 1.3021922, 1.1995943)

# Estimate dispersions and extract normalized tag counts
DDS=estimateDispersions(DDS)
Normalized <- varianceStabilizingTransformation(DDS, blind=F)
Normalized <- as.data.frame(assay(Normalized))
Normalized$PeakID <- rownames(Normalized)

# Average replicates
Normalized$SMC1_D0 <- rowMeans(2^Normalized[,c("SMC1_D0_exp1","SMC1_D0_exp2")])
Normalized$SMC1_H4 <- rowMeans(2^Normalized[,c("SMC1_4h_exp1","SMC1_4h_exp2")])
Normalized$SMC1_D2 <- rowMeans(2^Normalized[,c("SMC1_D2_exp1","SMC1_D2_exp2")])
Normalized$SMC1_D4 <- rowMeans(2^Normalized[,c("SMC1_D4_exp1","SMC1_D4_exp2")])
Normalized$SMC1_D7 <- rowMeans(2^Normalized[,c("SMC1_D7_exp1","SMC1_D7_exp2")])

# Plot the result
par(mfcol=c(1,5))
color_vec=c("#0000CD", "#0000CD","#006400","#006400","#CD0000")
for (clust in 1:5) {
	boxplot(Normalized[ Normalized$PeakID %in% Clustering[ Clustering$Cluster == clust, "PeakID"],c(12:16)], outline=F, names=c("D0","4h","D2","D4","D7"), las=1, ylab="SMC1 binding", xlab="Timepoint", col=color_vec[clust])
	}
```

[Back to start](../README.md)<br>
[Back to overview of Figure S3](../Links/FigureS3.md)