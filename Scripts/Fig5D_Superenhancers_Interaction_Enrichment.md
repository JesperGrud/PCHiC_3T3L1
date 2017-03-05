```R
# Load the necessary packages
library(edgeR)
library(GenomicRanges)
library(DESeq2)
library(circlize)
library(e1071)
library(fdrtool)

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

## Process the interaction data
# Import the interactions
Interactions <- read.delim("Data/Interactions/PCHiC/Interactions.txt", header = T)

# Remove interactions from or to blacklisted fragments
Blacklist <- read.delim("Data/Interactions/PCHiC/Annotation/Blacklisted_Fragments.txt", header=F)
FilteredInteractions <- Interactions[ !(Interactions$oeID %in% Blacklist[,1]),]
FilteredInteractions <- FilteredInteractions[ !(FilteredInteractions$baitID %in% Blacklist[,1]),]

# Filter away interactions in trans
FilteredInteractions <- FilteredInteractions[ FilteredInteractions$oeChr == FilteredInteractions$baitChr,]

# Filter away interactions spanning more than 1 megabase
FilteredInteractions <- FilteredInteractions[ abs(FilteredInteractions$dist) <= 1000000,]

# Remove baited fragments that are not genes, but ultra-conserved elements
UCE <- read.delim("Data/Interactions/PCHiC/Annotation/UCE_Fragments.txt", header=F)
FilteredInteractions <- FilteredInteractions[ !(FilteredInteractions$baitID %in% UCE[,1]),]

# Remove non-baited fragments that contain any type of gene
Genes <- read.delim("Data/Interactions/PCHiC/Annotation/Genes_Fragments.txt", header=T)
FilteredInteractions <- FilteredInteractions[ !(FilteredInteractions$oeID %in% Genes$Fragment),] 

# Remove baited fragments that does not contain a protein-coding gene
Genes <- Genes[ Genes$Type == "protein_coding",]
FilteredInteractions <- FilteredInteractions[ FilteredInteractions$baitID %in% Genes$Fragment,] 

# Normalize the read counts using edgeR
rownames(FilteredInteractions) <- FilteredInteractions$ID
D <- DGEList(as.matrix(FilteredInteractions[,c("Zero_R","Four_R","Two_R")]))
D <- calcNormFactors(D)
Normalized_Counts <- cpm(D, normalized.lib.sizes = T, log = T)
Normalized_Counts <- 2^Normalized_Counts
Normalized_Counts <- as.data.frame(Normalized_Counts)
colnames(Normalized_Counts) <- c("Zero_Norm","Four_Norm","Two_Norm")
Normalized_Counts$ID <- rownames(Normalized_Counts)
FilteredInteractions <- merge(FilteredInteractions, Normalized_Counts, by="ID")

# Calculate log2 fold changes between all samples using an additional pseudo-count of 0.25
FilteredInteractions$logFC_Four_Zero <- log2((FilteredInteractions$Four_Norm + 0.25) / (FilteredInteractions$Zero_Norm + 0.25))
FilteredInteractions$logFC_Two_Zero <- log2((FilteredInteractions$Two_Norm + 0.25) / (FilteredInteractions$Zero_Norm + 0.25))
FilteredInteractions$logFC_Two_Four <- log2((FilteredInteractions$Two_Norm + 0.25) / (FilteredInteractions$Four_Norm + 0.25))

# Split the interactions into three groups depending on the maximal absolute log2 fold change
Group1 <- FilteredInteractions[ abs(FilteredInteractions$logFC_Four_Zero) >= abs(FilteredInteractions$logFC_Two_Zero) & abs(FilteredInteractions$logFC_Four_Zero) >= abs(FilteredInteractions$logFC_Two_Four), ]
Group2 <- FilteredInteractions[ abs(FilteredInteractions$logFC_Two_Zero) >= abs(FilteredInteractions$logFC_Four_Zero) & abs(FilteredInteractions$logFC_Two_Zero) >= abs(FilteredInteractions$logFC_Two_Four), ]
Group3 <- FilteredInteractions[ abs(FilteredInteractions$logFC_Two_Four) >= abs(FilteredInteractions$logFC_Two_Zero) & abs(FilteredInteractions$logFC_Two_Four) >= abs(FilteredInteractions$logFC_Four_Zero), ]

# Set the true log2 fold change into new column for each group with different maxima
Group1$logFC_max <- Group1$logFC_Four_Zero
Group2$logFC_max <- Group2$logFC_Two_Zero
Group3$logFC_max <- Group3$logFC_Two_Four

# Combine the 3 groups into a new data.frame
Combined <- rbind(Group1, Group2, Group3)

# Select one maximal fold change for entries with duplicated maximal absolute log2 fold changes
set.seed(9999)
Combined$Random <- sample(c(1:nrow(Combined)), nrow(Combined), replace = T)
Combined <- Combined[ order(Combined$ID, Combined$Random),]
Combined <- Combined[ duplicated(Combined$ID)==F,]

# Merge maximal log2 fold change and the remaining data
FilteredInteractions <- merge(FilteredInteractions, Combined[,c("ID","logFC_max")], by="ID", all.x=T)

# Calculate the V values
FilteredInteractions$V <- apply(FilteredInteractions[,c("Zero_Norm","Four_Norm","Two_Norm")],1,FUN="max")

# Calculate the S values
FilteredInteractions$Minimum <- apply(FilteredInteractions[,c("Zero_Norm","Four_Norm","Two_Norm")],1,FUN="min")
FilteredInteractions$S <- 1 - (FilteredInteractions$Minimum / FilteredInteractions$V)

# Calculate the H values 
FilteredInteractions$H <- FilteredInteractions$Zero_Norm + FilteredInteractions$Four_Norm - FilteredInteractions$Two_Norm - FilteredInteractions$V
FilteredInteractions$H <- FilteredInteractions$H / (FilteredInteractions$V * FilteredInteractions$S)
FilteredInteractions$H <- FilteredInteractions$H + 2
FilteredInteractions$H <- FilteredInteractions$H * sign(FilteredInteractions$Four_Norm - FilteredInteractions$Zero_Norm)
FilteredInteractions$H <- FilteredInteractions$H * 60

# Devide V into equally sized groups between 0 and 1
Steps <- data.frame(Value = seq(1,0, by=-0.01), Step = 1 )
for (i in 1:101) { Steps[i,2] <- ceiling(nrow(FilteredInteractions)/101) * i }
FilteredInteractions <- FilteredInteractions[ order(-FilteredInteractions$V),]
First <- 1
for (i in 1:101) {
Last <- Steps[i,2]
FilteredInteractions[ c(First:Last),"V"] <- Steps[i,1]
First <- Last + 1
}

# Devide S into equally sized groups between 0 and 1
FilteredInteractions <- FilteredInteractions[ order(-FilteredInteractions$S),]
First <- 1
for (i in 1:101) {
Last <- Steps[i,2]
FilteredInteractions[ c(First:Last),"S"] <- Steps[i,1]
First <- Last + 1
}

# Calculate the color for each interaction
FilteredInteractions$color <- hsv(h = (FilteredInteractions$H + 180)/360, s = FilteredInteractions$S, v= FilteredInteractions$V, alpha = FilteredInteractions$V)

# Import all MED1 peaks
Peaks <- read.delim("Data/ChIPseq/Peak_Lists/MED1.bed", header=TRUE)

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
C1$PeakID2 <- paste("C1_",seq(1,nrow(C1),by=1),sep="")
C2$PeakID2 <- paste("C2_",seq(1,nrow(C2),by=1),sep="")
C3$PeakID2 <- paste("C3_",seq(1,nrow(C3),by=1),sep="")
C4$PeakID2 <- paste("C4_",seq(1,nrow(C4),by=1),sep="")
C5$PeakID2 <- paste("C5_",seq(1,nrow(C5),by=1),sep="")

## Overlap SE constituents and otherEnds
# Make a GRanges object from otherEnds
otherEnds <- FilteredInteractions[ duplicated(FilteredInteractions$oeID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$oeChr, IRanges(start = otherEnds$oeStart, end = otherEnds$oeEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$oeID)
	
# Make a list with all constituents
Constituents <- list(C1, C2, C3, C4, C5)

# Define bin size for calculating enrichments
Binsize <- 40
WindowSize <- 1
Length <- length(seq(-180,180-WindowSize,by=WindowSize))

# Make a list for capturing the results
ResultList <- list()

# Loop through all sets and plot
for (consti in 1:5) {
	# Grab the data from the list
	tmp <- Constituents[[consti]]
	# Create GRanges object for enhancers and otherEnds
	ConstRanges <- GRanges(seqnames = tmp$Chr, IRanges(start = tmp$Start, end = tmp$End), strand = rep("+",nrow(tmp)), mcols = tmp$PeakID2)
	# Get overlap and width of overlap, merge into a single data.frame and collaps
	ConstOverlap <- suppressWarnings(findOverlaps(otherEndRanges, ConstRanges))
	ConstWidth <- as.data.frame(ranges(ConstOverlap, ranges(otherEndRanges), ranges(ConstRanges)))
	ConstOverlap <- as.data.frame(ConstOverlap)
	# Convert overlap to a single data frame
	ConstOverlap <- cbind(as.data.frame(otherEndRanges[ConstOverlap[,1]]),as.data.frame(ConstRanges[ConstOverlap[,2]]),ConstWidth)
	# Collaps overlap to only the largest overlap
	ConstOverlap <- ConstOverlap[ order(ConstOverlap[,12], -ConstOverlap[,15]),]
	ConstOverlap <- ConstOverlap[ duplicated(ConstOverlap[,12])==F,]
	# Set column names
	colnames(ConstOverlap)[6] <- "oeID"
	# Make a data frame to capture the enrichment results
	Enrichment <- data.frame(matrix(ncol=2, nrow= Length))
	# Set column names and define upper and lower bounds
	colnames(Enrichment) <- c("Enrichment","Pvalue")
	HueLower <- seq(-180,180-Binsize,by=WindowSize)
  HueUpper <- seq(-180+Binsize,180,by=WindowSize)
	# Calculate enrichment and P-values in a loop
	for (i in 1:length(HueLower)) {
	  # Define the hue range
	  Lower <- HueLower[i]
		Upper <- HueUpper[i]
	  # Devide interaction into a group within the hue range and a group outside
	  Within <- FilteredInteractions[ FilteredInteractions$H > Lower & FilteredInteractions$H < Upper,]
		Outside <- FilteredInteractions[ FilteredInteractions$H < Lower | FilteredInteractions$H > Upper,]
    # Get the interactions within the hue range overlapping and not overlapping constituents
		WithinOverlap <- Within[ Within$oeID %in% ConstOverlap$oeID,]
		WithinNoOverlap <- Within[ !(Within$oeID %in% ConstOverlap$oeID),]
		# Get the interactions outside the hue range overlapping and not overlapping constituents
		OutsideOverlap <- Outside[ Outside$oeID %in% ConstOverlap$oeID,]
		OutsideNoOverlap <- Outside[ !(Outside$oeID %in% ConstOverlap$oeID),]
		# Calculate enrichment and p-values
		Enrichment[i,"Enrichment"] <- log2(((nrow(WithinOverlap)+1)/(nrow(WithinOverlap)+nrow(WithinNoOverlap)+1))/((nrow(OutsideOverlap)+1)/(nrow(OutsideOverlap)+nrow(OutsideNoOverlap)+1)))
		Enrichment[i,"Pvalue"] <- binom.test(x = nrow(WithinOverlap), n = nrow(WithinOverlap)+nrow(WithinNoOverlap), p = nrow(OutsideOverlap)/(nrow(OutsideOverlap)+nrow(OutsideNoOverlap)))$p.value
	}
  for (q in 1:(Binsize-WindowSize)/WindowSize) {
    i <- length(HueLower) + q
    Lower <- max(HueLower) + (q * WindowSize)
	  Upper <- -180 + (q * WindowSize)
	  Within <- FilteredInteractions[ FilteredInteractions$H > Lower | FilteredInteractions$H < Upper,]
	  Outside <- FilteredInteractions[ FilteredInteractions$H < Lower & FilteredInteractions$H > Upper,]
	  # Get interactions inside and outside gene clusters
	  WithinOverlap <- Within[ Within$oeID %in% ConstOverlap$oeID,]
	  WithinNoOverlap <- Within[ !(Within$oeID %in% ConstOverlap$oeID),]
		# Get the interactions outside the hue range overlapping and not overlapping constituents
		OutsideOverlap <- Outside[ Outside$oeID %in% ConstOverlap$oeID,]
		OutsideNoOverlap <- Outside[ !(Outside$oeID %in% ConstOverlap$oeID),]
		# Calculate enrichment and p-values
		Enrichment[i,"Enrichment"] <- log2(((nrow(WithinOverlap)+1)/(nrow(WithinOverlap)+nrow(WithinNoOverlap)+1))/((nrow(OutsideOverlap)+1)/(nrow(OutsideOverlap)+nrow(OutsideNoOverlap)+1)))
		Enrichment[i,"Pvalue"] <- binom.test(x = nrow(WithinOverlap), n = nrow(WithinOverlap)+nrow(WithinNoOverlap), p = nrow(OutsideOverlap)/(nrow(OutsideOverlap)+nrow(OutsideNoOverlap)))$p.value
	 }
	# Place the result into a new list
	ResultList[[consti]] <- Enrichment
	print(consti)
	}

# Perform FDR correction
for (i in 1:5) {
PBinom[i,] <- fdrtool(as.vector(as.matrix(PBinom[i,])), statistic = "pvalue", plot=F, verbose=F)$qval
}

## Plot the enrichments in a loop
par(mfcol=c(2,3))
Limits <- c(1.75,2,1.5,1.5,3,1.5)
for (q in 1:5) {
	# Grab the data
	tmp <- ResultList[[q]]
	# Set min and max
	Max <- Limits[q]
	Min <- -1 * Max
	# Initialize the plot
	circos.clear()
	circos.par(gap.degree=0, start.degree=-90)
	circos.initialize("A", xlim=c(-179,179))
  # Plot the raw enrichments that dont cross
	circos.trackPlotRegion("A", ylim=c(Min,Max), track.height=0.6, bg.border=NA)	
	for (i in 1:length(HueLower)) {circos.rect(HueLower[i]+19.5, 0, HueUpper[i]-19.5, tmp[i,1], col = ifelse(tmp[i,2] <= 0.1, ifelse(tmp[i,1] < 0, "darkred","darkgreen"),"grey"), border = NA) }
	for (m in 1:(Binsize-WindowSize)/WindowSize) {
  Lower <- max(HueLower) + (m * WindowSize)+19.5
	Upper <- max(HueLower) + (m * WindowSize)+20.5
  circos.rect(Lower, 0, Upper, tmp[(length(HueLower)+m),1], col = ifelse(tmp[(length(HueLower)+m),2] <= 0.1,ifelse(tmp[(length(HueLower)+m),1] < 0, "darkred","darkgreen"),"grey"), border = NA)
	}
	circos.lines(c(-180,180),c(0,0))
	circos.lines(c(-180,180),c(Max,Max),lty=2)
	circos.lines(c(-180,180),c(Min,Min),lty=2)
	circos.lines(c(max(Result[ Result$Cluster == q, "H"])+10,max(Result[ Result$Cluster == q, "H"])+10),c(Min,Max), col="blue", lwd=2)
	 circos.lines(c(min(Result[ Result$Cluster == q, "H"])-10,min(Result[ Result$Cluster == q, "H"])-10),c(Min,Max), col="blue", lwd=2)
	circos.trackPlotRegion("A", ylim=c(Min,Max), track.height=0.2, bg.border=NA)
	circos.text(0,-3,paste("Cluster",q,sep=" "))
	}
```

[Back to start](../README.md)<br>
[Back to overview of Figure 5](../Links/Figure5.md)
