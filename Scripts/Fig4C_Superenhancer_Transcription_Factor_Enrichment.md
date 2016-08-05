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

## Process super-enhancer constitutents
Peaks <- read.delim("Data/ChIPseq/Peak_lists/MED1.bed", header=TRUE)

# Overlap MED1 peaks and super-enhancers
PeakRanges <- GRanges(seqnames = Peaks$Chr, ranges = IRanges(start = Peaks$Start, end = Peaks$End), strand = rep("+",nrow(Peaks)), mcols = Peaks$PeakID)
SERanges <- GRanges(seqnames = Result$Chr, ranges = IRanges(start = Result$Start, end = Result$End), strand = rep("+",nrow(Result)), mcols = Result$PeakID)
Overlap <- suppressWarnings(findOverlaps(PeakRanges, SERanges))
Overlap <- as.data.frame(cbind(as.character(as.data.frame(PeakRanges[as.matrix(Overlap)[,1]])[,6]),as.character(as.data.frame(SERanges[as.matrix(Overlap)[,2]])[,6])))

# Read the data
C1 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Clustering[ Clustering$Cluster == 1, "PeakID"],1],]
C2 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Clustering[ Clustering$Cluster == 2, "PeakID"],1],]
C3 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Clustering[ Clustering$Cluster == 3, "PeakID"],1],]
C4 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Clustering[ Clustering$Cluster == 4, "PeakID"],1],]
C5 <- Peaks[ Peaks$PeakID %in% Overlap[ Overlap[,2] %in% Clustering[ Clustering$Cluster == 5, "PeakID"],1],]

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
CEBPB_Preadip <- read.table("data/CEBPb_D0.Siers.TD.findPeaks.bed", quote="\"")
PPARG_Late <- read.table("Data/ChIPseq/Peak_lists/PPARg_D6.bed", quote="\"")
CEBPA_Late <- read.table("Data/ChIPseq/Peak_lists/CEBPa_Late_SRR1233891.bed", quote="\"")

# Place all peaks into a list
Factors <- list(CEBPB_Preadip,KLF4_Preadip,CEBPB_Early,GR_Early,JUNB_Early,KLF4_Early,KLF5_Early,STAT5A_Early,VDR_Early,PPARG_Late,CEBPA_Late)
Names <- c("CEBPb","KLF4","CEBPb","GR","JunB","KLF4","KLF5","STAT5A","VDR","PPARg","CEBPa")

# Make a list to save the results in
Enrichments <- list()

# Make a GRange object for all super-enhancer constituents
AllConstituents <- rbind(C1, C2, C3, C4, C5)
AllConstituentRanges <- GRanges(seqnames = AllConstituents[,1], IRanges(start = AllConstituents[,2], end = AllConstituents[,3]), strand = rep("+",nrow(AllConstituents)), mcols = AllConstituents[,4])

# Overlap all constituents as well as each cluster of constituents with each factor
for (consti in 1:5) {
	# Make a data.frame to capture the overlap
	Result <- data.frame(matrix(ncol=4,nrow=11))
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
		Result[fact,1] <- nrow(Overlap)
		# Get the TF sites that overlap with a constituent in all the super-enhancer constituents
		Overlap <- suppressWarnings(findOverlaps(AllConstituentRanges, FactorRanges))
		Overlap <- as.data.frame(Overlap)
		Overlap <- as.data.frame(FactorRanges[Overlap[,2]])
		Overlap <- Overlap[ duplicated(Overlap$mcols)==F,]
		Result[fact,2] <- nrow(Overlap)
		# Record the number of constituents in each group
		Result[fact,3] <- nrow(as.data.frame(ConstituentRanges))
		Result[fact,4] <- nrow(as.data.frame(AllConstituentRanges))
		}
	# Calculte enrichment
	Result[,5] <- log2((Result[,1]/Result[,3])/(Result[,2]/Result[,4]))
	# Place into list
	Enrichments[[consti]] <- Result
	}
	
# Plot the results
par(mfcol=c(1,5))
for (i in 1:5) {
	Result <- Enrichments[[i]]
	Min <- round(min(Result[,5]),0)-0.5
	Max <- round(max(Result[,5]),0)+0.5
	barplot(Result[,5], horiz = T, xlim=c(Min,Max), col=c(rep("red3",2),rep("blue3",7), rep("green3",2)), names=Names, las=2)
}
```

[Back to start](../README.md)<br>
[Back to overview of Figure 4](../Links/Figure4.md)