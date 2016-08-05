```R
## Load the necessary packages
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
Result[ apply(Result[,c(15:34)],1,FUN="min") <= 0.05,"Differential"] <- 1

# Calculate average occupancy in replicates
Result$Count_D0 <- rowMeans(Result[,c(5,10)])
Result$Count_H4 <- rowMeans(Result[,c(6,11)])
Result$Count_D2 <- rowMeans(Result[,c(7,12)])
Result$Count_D4 <- rowMeans(Result[,c(8,13)])
Result$Count_D7 <- rowMeans(Result[,c(9,14)])

# Center the averaged counts (only D0, H4 and D2)
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

## Process super-enhancer constitutents
Peaks <- read.delim("Data/ChIPseq/Peak_lists/MED1.bed", header=TRUE)

# Overlap MED1 peaks and super-enhancers
PeakRanges <- GRanges(seqnames = Peaks$Chr, ranges = IRanges(start = Peaks$Start, end = Peaks$End), strand = rep("+",nrow(Peaks)), mcols = Peaks$PeakID)
SERanges <- GRanges(seqnames = Result$Chr, ranges = IRanges(start = Result$Start, end = Result$End), strand = rep("+",nrow(Result)), mcols = Result$PeakID)
Overlap <- suppressWarnings(findOverlaps(PeakRanges, SERanges))
Overlap <- as.data.frame(cbind(as.character(as.data.frame(PeakRanges[as.matrix(Overlap)[,1]])[,6]),as.character(as.data.frame(SERanges[as.matrix(Overlap)[,2]])[,6])))

# Define peaks based on tag counts and fold over input
Peaks$MED_D0 <- 0
Peaks$MED_H4 <- 0
Peaks$MED_D2 <- 0
Peaks$MED_D4 <- 0
Peaks$MED_D7 <- 0
for (i in 6:10) { Peaks[ (Peaks[,i] / Peaks[,5] >= 4 & Peaks[,i] >= 24) & (Peaks[,(i+5)] / Peaks[,5] >= 4 & Peaks[,(i+5)] >= 24),((i-5)+15)] <- 1 }

# Make a list to capture the counts
ResultList <- list()

# For each cluster, count the number of super-enhancer regions and 'regular regions' at each timepoint
for (i in 1:5) {
	tmp <- data.frame(matrix(ncol=5, nrow=2))
	rownames(tmp) <- c("SE","RR")
	colnames(tmp) <- c("D0","H4","D2","D4","D7")
	tmp[1,1] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$D0 == 1,])
	tmp[1,2] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$H4 == 1,])
	tmp[1,3] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$D2 == 1,])
	tmp[1,4] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$D4 == 1,])
	tmp[1,5] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$D7 == 1,])
	tmp[2,1] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$D0 == 0 & Result$PeakID %in% Overlap[ Overlap[,1] %in% Peaks[ Peaks$MED_D0 == 1, "PeakID"],2],])
	tmp[2,2] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$H4 == 0 & Result$PeakID %in% Overlap[ Overlap[,1] %in% Peaks[ Peaks$MED_H4 == 1, "PeakID"],2],])
	tmp[2,3] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$D2 == 0 & Result$PeakID %in% Overlap[ Overlap[,1] %in% Peaks[ Peaks$MED_D2 == 1, "PeakID"],2],])
	tmp[2,4] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$D4 == 0 & Result$PeakID %in% Overlap[ Overlap[,1] %in% Peaks[ Peaks$MED_D4 == 1, "PeakID"],2],])
	tmp[2,5] <- nrow(Result[ Result$Differential == 1 & Result$PeakID %in% Clustering[ Clustering$Cluster == i, "PeakID"] & Result$D7 == 0 & Result$PeakID %in% Overlap[ Overlap[,1] %in% Peaks[ Peaks$MED_D7 == 1, "PeakID"],2],])
	ResultList[[i]] <- tmp
	}

# Plot the results
Limits <- c(350,200,300,200,100)
par(mfcol=c(1,5))
for (i in 1:5) {
	tmp <- ResultList[[i]]
	barplot(as.matrix(tmp), col=c("blue4","grey"), las=1, ylim=c(0,Limits[i]), ylab="# regions")
	}
```

[Back to start](../README.md)<br>
[Back to overview of Figure S3](../Links/FigureS3.md)