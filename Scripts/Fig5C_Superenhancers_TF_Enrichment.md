```R
## Load the necessary packages
library(GenomicRanges)
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
Result[ Result[,20] <= 0.05 | Result[,23] <= 0.05 | Result[,25] <= 0.05,"Differential"] <- 1

# Calculate average occupancy in replicates
Result$Count_D0 <- rowMeans(2^Result[,c("MED1_D0_exp1","MED1_D0_exp2")])
Result$Count_H4 <- rowMeans(2^Result[,c("MED1_4h_exp1","MED1_4h_exp2")])
Result$Count_D2 <- rowMeans(2^Result[,c("MED1_D2_exp1","MED1_D2_exp2")])
Result$Count_D4 <- rowMeans(2^Result[,c("MED1_D4_exp1","MED1_D4_exp2")])
Result$Count_D7 <- rowMeans(2^Result[,c("MED1_D7_exp1","MED1_D7_exp2")])

# Partition data into 5 cluster using fuzzy clustering for the differential super-enhancers (on both replicates)
rownames(Result) <- Result$PeakID
ClusterData <- Result[ Result$Differential == 1 & (Result$D0 == 1 | Result$D2 == 1 | Result$H4 == 1), c("Count_D0","Count_H4","Count_D2")]
ClusterData <- ClusterData
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
Rep1 <- as.matrix(Result[ Result$Differential == 1 ,c("MED1_D0_exp1","MED1_4h_exp1","MED1_D2_exp1")])
Rep1 <- t(scale(t(Rep1)))
dm <- sapply(seq_len(nrow(Rep1)),function(i) apply(Centers, 1, function(v) sqrt(sum((Rep1[i, ]-v)^2))))
Rep1Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Rep1Membership <- as.data.frame(Rep1Membership)
rownames(Rep1Membership) <- Result[ Result$Differential == 1,c("PeakID")]

Rep2 <- as.matrix(Result[ Result$Differential == 1 ,c("MED1_D0_exp2","MED1_4h_exp2","MED1_D2_exp2")])
Rep2 <- t(scale(t(Rep2)))
dm <- sapply(seq_len(nrow(Rep2)),function(i) apply(Centers, 1, function(v) sqrt(sum((Rep2[i, ]-v)^2))))
Rep2Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Rep2Membership <- as.data.frame(Rep2Membership)
rownames(Rep2Membership) <- Result[ Result$Differential == 1,c("PeakID")]

All <- as.matrix(Result[ Result$Differential == 1 ,c("Count_D0","Count_H4","Count_D2")])
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
Rep1Membership[ Rep1Membership$Max >= 0.2,"Include"] <- 1
Rep1Membership$Which <- 0
for (i in 1:nrow(Rep1Membership)) { Rep1Membership[i,"Which"] <- which.max(Rep1Membership[i,c(1:(NoClusters))]) }
Rep1Membership$Which <- Rep1Membership$Which * Rep1Membership$Include
Rep1Membership$PeakID <- rownames(Rep1Membership)

Rep2Membership$Max <- apply(Rep2Membership,1,FUN="max")
Rep2Membership$Include <- 0
Rep2Membership[ Rep2Membership$Max >= 0.2,"Include"] <- 1
Rep2Membership$Which <- 0
for (i in 1:nrow(Rep2Membership)) { Rep2Membership[i,"Which"] <- which.max(Rep2Membership[i,c(1:(NoClusters))]) }
Rep2Membership$Which <- Rep2Membership$Which * Rep2Membership$Include
Rep2Membership$PeakID <- rownames(Rep2Membership)

Membership$Max <- apply(Membership,1,FUN="max")
Membership$Include <- 0
Membership[ Membership$Max >= 0.3,"Include"] <- 1
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

## Filter-based on the hue of cluster and the three time points
# Calculate the hue of the individual super-enhancers
Result$V <- apply(Result[,c("Count_D0","Count_H4","Count_D2")],1,FUN="max")

# Calculate the S values
Result$Minimum <- apply(Result[,c("Count_D0","Count_H4","Count_D2")],1,FUN="min")
Result$S <- 1 - (Result$Minimum / Result$V)

# Calculate the H values 
Result$H <- Result$Count_D0 + Result$Count_H4 - Result$Count_D2 - Result$V
Result$H <- Result$H / (Result$V * Result$S)
Result$H <- Result$H + 2
Result$H <- Result$H * sign(Result$Count_H4 - Result$Count_D0)
Result$H <- Result$H * 60

# Calculate the hue of the centers
Centers <- as.data.frame(Centers)
Centers$V <- apply(Centers[,c("Count_D0","Count_H4","Count_D2")],1,FUN="max")

# Calculate the S values
Centers$Minimum <- apply(Centers[,c("Count_D0","Count_H4","Count_D2")],1,FUN="min")
Centers$S <- 1 - (Centers$Minimum / Centers$V)

# Calculate the H values 
Centers$H <- Centers$Count_D0 + Centers$Count_H4 - Centers$Count_D2 - Centers$V
Centers$H <- Centers$H / (Centers$V * Centers$S)
Centers$H <- Centers$H + 2
Centers$H <- Centers$H * sign(Centers$Count_H4 - Centers$Count_D0)
Centers$H <- Centers$H * 60

# Filter based on the hue
for (i in 1:NoClusters) {
  Hue <- Centers[i,"H"]  
  Tmp <- Result[ Result$Cluster == i,]
  Tmp$Distance <- 0
  for (q in 1:nrow(Tmp)) { Tmp[q,"Distance"] <- min(c(abs(Hue - Tmp[q,"H"]),(180-abs(Tmp[q,"H"]))+(180-abs(Hue))))  }
  Result[ Result$PeakID %in% Tmp[ Tmp$Distance > 20,"PeakID"],"Cluster"] <- 0
}

# Rename clusters for intuitive cluster profiles
Result[ Result$Cluster == 1, "Cluster"] <- 7
Result[ Result$Cluster == 3, "Cluster"] <- 8
Result[ Result$Cluster == 2, "Cluster"] <- 9
Result[ Result$Cluster == 5, "Cluster"] <- 10
Result[ Result$Cluster == 4, "Cluster"] <- 11
Result[ Result$Cluster == 7, "Cluster"] <- 1
Result[ Result$Cluster == 8, "Cluster"] <- 2
Result[ Result$Cluster == 9, "Cluster"] <- 3
Result[ Result$Cluster == 10, "Cluster"] <- 4
Result[ Result$Cluster == 11, "Cluster"] <- 5

## Process super-enhancer constitutents
Peaks <- read.delim("Data/ChIPseq/Peak_lists/MED1.bed", header=TRUE)

# Overlap MED1 peaks and super-enhancers
PeakRanges <- GRanges(seqnames = Peaks$Chr, ranges = IRanges(start = Peaks$Start, end = Peaks$End), strand = rep("+",nrow(Peaks)), mcols = Peaks$PeakID)
SERanges <- GRanges(seqnames = Result$Chr, ranges = IRanges(start = Result$Start, end = Result$End), strand = rep("+",nrow(Result)), mcols = Result$PeakID)
Overlap <- suppressWarnings(findOverlaps(PeakRanges, SERanges))
Overlap <- as.data.frame(cbind(as.character(as.data.frame(PeakRanges[as.matrix(Overlap)[,1]])[,6]),as.character(as.data.frame(SERanges[as.matrix(Overlap)[,2]])[,6])))

# Read the data
C1 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Result[ Result$Cluster == 1, "PeakID"],1],]
C2 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Result[ Result$Cluster == 2, "PeakID"],1],]
C3 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Result[ Result$Cluster == 3, "PeakID"],1],]
C4 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Result[ Result$Cluster == 4, "PeakID"],1],]
C5 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Result[ Result$Cluster == 5, "PeakID"],1],]

# Define unique PeakID
C1$PeakID <- paste("C1_",seq(1,nrow(C1),by=1),sep="")
C2$PeakID <- paste("C2_",seq(1,nrow(C2),by=1),sep="")
C3$PeakID <- paste("C3_",seq(1,nrow(C3),by=1),sep="")
C4$PeakID <- paste("C4_",seq(1,nrow(C4),by=1),sep="")
C5$PeakID <- paste("C5_",seq(1,nrow(C5),by=1),sep="")

# Place into a list
Constituents <- list(C1, C2, C3, C4, C5)

## Process ChIP-seq data
# Import all peaks
CEBPB_Early <- read.table("Data/ChIPseq/Peak_lists/CEBPb_4h.bed", quote="\"")
VDR_Early <- read.table("Data/ChIPseq/Peak_lists/VDR_4h.bed", quote="\"")
STAT5A_Early <- read.table("Data/ChIPseq/Peak_lists/Stat5a_4h.bed", quote="\"")
JUNB_Early <- read.table("Data/ChIPseq/Peak_lists/JunB_4h.bed", quote="\"")
KLF5_Early <- read.table("Data/ChIPseq/Peak_lists/KLF5_4h.bed", quote="\"")
KLF4_Early <- read.table("Data/ChIPseq/Peak_lists/KLF4_4h.bed", quote="\"")
GR_Early <- read.table("Data/ChIPseq/Peak_lists/GR_4h.bed", quote="\"")
KLF4_Preadip <- read.table("Data/ChIPseq/Peak_lists/KLF4_D0.bed", quote="\"")
CEBPB_Preadip <- read.table("Data/ChIPseq/Peak_lists/CEBPb_D0.bed", quote="\"")
PPARG_Late <- read.table("Data/ChIPseq/Peak_lists/PPARg_D6.bed", quote="\"")
CEBPA_Late <- read.table("Data/ChIPseq/Peak_lists/CEBPa_Late_SRR1233891.bed", quote="\"")

# Place all peaks into a list
Factors <- list(CEBPB_Preadip,KLF4_Preadip,CEBPB_Early,GR_Early,JUNB_Early,KLF4_Early,KLF5_Early,STAT5A_Early,VDR_Early,PPARG_Late,CEBPA_Late)
Names <- c("CEBPb","KLF4","CEBPb","GR","JunB","KLF4","KLF5","STAT5A","VDR","PPARg","CEBPa")

# Make a list to save the results in
Enrichments <- list()

# Make a GRange object for all super-enhancer constituents
AllConstituents <- Peaks[ Peaks$PeakID %in% Overlap[,1],]
AllConstituentRanges <- GRanges(seqnames = AllConstituents[,1], IRanges(start = AllConstituents[,2], end = AllConstituents[,3]), strand = rep("+",nrow(AllConstituents)), mcols = AllConstituents[,4])

# Overlap all constituents as well as each cluster of constituents with each factor
for (consti in 1:5) {
	# Make a data.frame to capture the overlap
	TmpResult <- data.frame(matrix(ncol=4,nrow=11))
	# Grab constituent sites
	ConstituentSites <- Constituents[[consti]]
	# Make a GRanges object
	ConstituentRanges <- GRanges(seqnames = ConstituentSites[,1], IRanges(start = ConstituentSites[,2], end = ConstituentSites[,3]), strand = rep("+",nrow(ConstituentSites)), mcols = ConstituentSites[,4])
	for (fact in 1:11) {
		# Grab factor sites
		FactorSites <- Factors[[fact]]
		# Make a GRange objects
		FactorRanges <- GRanges(seqnames = FactorSites[,1], IRanges(start = FactorSites[,2], end = FactorSites[,3]), strand = rep("+",nrow(FactorSites)), mcols = FactorSites[,4])
		# Get the TF sites that overlap with a constituent in the constituent cluster
		Overlap <- suppressWarnings(findOverlaps(ConstituentRanges, FactorRanges))
		Overlap <- as.data.frame(Overlap)
		Overlap <- as.data.frame(FactorRanges[Overlap[,2]])
		Overlap <- Overlap[ duplicated(Overlap$mcols)==F,]
		TmpResult[fact,1] <- nrow(Overlap)
		# Get the TF sites that overlap with a constituent in all the super-enhancer constituents
		Overlap <- suppressWarnings(findOverlaps(AllConstituentRanges, FactorRanges))
		Overlap <- as.data.frame(Overlap)
		Overlap <- as.data.frame(FactorRanges[Overlap[,2]])
		Overlap <- Overlap[ duplicated(Overlap$mcols)==F,]
		TmpResult[fact,2] <- nrow(Overlap)
		# Record the number of constituents in each group
		TmpResult[fact,3] <- nrow(as.data.frame(ConstituentRanges))
		TmpResult[fact,4] <- nrow(as.data.frame(AllConstituentRanges))
		}
	# Calculte enrichment
	TmpResult[,5] <- log2((TmpResult[,1]/TmpResult[,3])/(TmpResult[,2]/TmpResult[,4]))
	# Place into list
	Enrichments[[consti]] <- TmpResult
	}
	
# Plot the results
par(mfcol=c(1,5))
for (i in 1:5) {
	TmpResult <- Enrichments[[i]]
	Min <- round(min(TmpResult[,5]),0)-0.5
	Max <- round(max(TmpResult[,5]),0)+0.5
	barplot(TmpResult[,5], horiz = T, xlim=c(Min,Max), col=c(rep("red3",2),rep("blue3",7), rep("green3",2)), names=Names, las=2, main=paste("Cluster",i, sep=" "))
}
```

[Back to start](../README.md)<br>
[Back to overview of Figure 5](../Links/Figure5.md)