```R
## Load the necessary libraries
library(GenomicRanges)
library(pROC)
library(edgeR)

## Process ChIP-seq data
# Read count data
Broad <- read.delim("Data/ChIPseq/Counts/Prediction_Histone_Marks.count")
Point <- read.delim("Data/ChIPseq/Counts/Prediction_Cofactors.count")

# Extract peak location and the names of peak-detected factors at each site
Peaks <- Point[,c(1:4,7)]
colnames(Peaks) <- c("PeakID","Chr","Start","End","Factors")

# Keep only PeakID and counts for both broad and point marks
Broad <- Broad[,c(1,8:19)]
colnames(Broad) <- c("PeakID","H3K27ac_4h_E1","H3K27ac_4h_E2","H3K27ac_D0_E1","H3K27ac_D0_E2","H3K4me1_4h_E1","H3K4me1_4h_E2","H3K4me1_D0_E1","H3K4me1_D0_E2","H3K4me2_4h_E1","H3K4me2_4h_E2","H3K4me2_D0_E1","H3K4me2_D0_E2")
Point <- Point[,c(1,8:31)]
colnames(Point) <- c("PeakID","MED1_D0_E1","MED1_D0_E2","MED1_4h_E1","MED1_4h_E2","SMC1_D0_E1","SMC1_D0_E2","SMC1_4h_E1","SMC1_4h_E2","P300_4h_E1","P300_4h_E2","P300_D0_E1","P300_D0_E2","NCoR_4h_E1","NCoR_4h_E2","NCoR_D0_E1","NCoR_D0_E2","HDAC2_4h_E1","HDAC2_4h_E2","HDAC2_D0_E1","HDAC2_D0_E2","HDAC3_4h_E1","HDAC3_4h_E2","HDAC3_D0_E1","HDAC3_D0_E2")

# Merge count data and clean up
Counts <- merge(Point, Broad, by="PeakID")
rm(list=c("Broad","Point"))

# Take average for replicates
Counts$HDAC3_4h <- rowMeans(Counts[,c("HDAC3_4h_E1","HDAC3_4h_E2")])
Counts$HDAC3_D0 <- rowMeans(Counts[,c("HDAC3_D0_E1","HDAC3_D0_E2")])
Counts$HDAC2_4h <- rowMeans(Counts[,c("HDAC2_4h_E1","HDAC2_4h_E2")])
Counts$HDAC2_D0 <- rowMeans(Counts[,c("HDAC2_D0_E1","HDAC2_D0_E2")])
Counts$NCoR_4h <- rowMeans(Counts[,c("NCoR_4h_E1","NCoR_4h_E2")])
Counts$NCoR_D0 <- rowMeans(Counts[,c("NCoR_D0_E1","NCoR_D0_E2")])
Counts$SMC1_4h <- rowMeans(Counts[,c("SMC1_4h_E1","SMC1_4h_E2")])
Counts$SMC1_D0 <- rowMeans(Counts[,c("SMC1_D0_E1","SMC1_D0_E2")])
Counts$MED1_4h <- rowMeans(Counts[,c("MED1_4h_E1","MED1_4h_E2")])
Counts$MED1_D0 <- rowMeans(Counts[,c("MED1_D0_E1","MED1_D0_E2")])
Counts$P300_4h <- rowMeans(Counts[,c("P300_4h_E1","P300_4h_E2")])
Counts$P300_D0 <- rowMeans(Counts[,c("P300_D0_E1","P300_D0_E2")])
Counts$H3K27ac_4h <- rowMeans(Counts[,c("H3K27ac_4h_E1","H3K27ac_4h_E2")])
Counts$H3K27ac_D0 <- rowMeans(Counts[,c("H3K27ac_D0_E1","H3K27ac_D0_E2")])
Counts$H3K4me1_4h <- rowMeans(Counts[,c("H3K4me1_4h_E1","H3K4me1_4h_E2")])
Counts$H3K4me1_D0 <- rowMeans(Counts[,c("H3K4me1_D0_E1","H3K4me1_D0_E2")])
Counts$H3K4me2_4h <- rowMeans(Counts[,c("H3K4me2_4h_E1","H3K4me2_4h_E2")])
Counts$H3K4me2_D0 <- rowMeans(Counts[,c("H3K4me2_D0_E1","H3K4me2_D0_E2")])

# Keep only columns with replicate averages
Counts <- Counts[,c(1,38:55)]

# Calculate log2 fold change for each mark using 0.25 as a pseudocount
Counts$P300 <- log2((Counts$P300_4h+0.25)/(Counts$P300_D0+0.25))
Counts$HDAC2 <- log2((Counts$HDAC2_4h+0.25)/(Counts$HDAC2_D0+0.25))
Counts$HDAC3 <- log2((Counts$HDAC3_4h+0.25)/(Counts$HDAC3_D0+0.25))
Counts$NCoR <- log2((Counts$NCoR_4h+0.25)/(Counts$NCoR_D0+0.25))
Counts$SMC1 <- log2((Counts$SMC1_4h+0.25)/(Counts$SMC1_D0+0.25))
Counts$MED1 <- log2((Counts$MED1_4h+0.25)/(Counts$MED1_D0+0.25))
Counts$H3K27ac <- log2((Counts$H3K27ac_4h+0.25)/(Counts$H3K27ac_D0+0.25))
Counts$H3K4me1 <- log2((Counts$H3K4me1_4h+0.25)/(Counts$H3K4me1_D0+0.25))
Counts$H3K4me2 <- log2((Counts$H3K4me2_4h+0.25)/(Counts$H3K4me2_D0+0.25))

# Keep only columns with log2 fold changes
Counts <- Counts[,c(1,20:28)]

## Process PCHiC data
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

## Overlap interaction and ChIP-seq data
# Create GRanges object for enhancers and otherEnds
EnhancerRanges <- GRanges(seqnames = Peaks$Chr, IRanges(start = Peaks$Start, end = Peaks$End), strand = rep("+",nrow(Peaks)), mcols = Peaks$PeakID)
otherEnds <- FilteredInteractions[ duplicated(FilteredInteractions$oeID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$oeChr, IRanges(start = otherEnds$oeStart, end = otherEnds$oeEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$oeID)

# Get overlap and width of overlap, merge into a single data.frame and collaps
Overlap <- suppressWarnings(findOverlaps(otherEndRanges, EnhancerRanges))
Width <- as.data.frame(ranges(Overlap, ranges(otherEndRanges), ranges(EnhancerRanges)))
Overlap <- as.data.frame(Overlap)

# Convert to data frames
otherEnds <- as.data.frame(otherEndRanges[Overlap[,1]])
Enhancers <- as.data.frame(EnhancerRanges[Overlap[,2]])

# Set column names
colnames(otherEnds)[6] <- "oeID"
colnames(Enhancers)[6] <- "PeakID"

# Make an ID for merging
otherEnds$Merge <- seq(1, nrow(otherEnds), by=1)
Enhancers$Merge <- seq(1, nrow(Enhancers), by=1)
Width$Merge <- seq(1, nrow(Width), by=1)

# Merge the overlap results
Overlap <- merge(otherEnds[,c(6,7)], Enhancers[,c(6,7)])
Overlap <- merge(Overlap, Width[,c(3,4)], by="Merge")

# Collaps overlap to only the largest overlap
Overlap <- Overlap[ order(Overlap$PeakID, -Overlap$width),]
Overlap <- Overlap[ duplicated(Overlap$PeakID)==F,]
Overlap <- Overlap[,c(2,3)]

# Merge interaction information onto overlap
Overlap <- merge(Overlap, FilteredInteractions[,c("oeID","ID","logFC_Four_Zero")], by="oeID")

## Filter interactions by distance
FilteredInteractions <- FilteredInteractions[ abs(FilteredInteractions$dist) <= 250000,]

## Define a response variable (|FC| >= 1.5: -1 = Repressed interactions, 1 = Induced interactions)
## Only include otherEnds where all interactions are dynamic in the same direction or unchanged
# Define count variable on all interactions, summarize and set column names
FilteredInteractions$Count <- 1
Table1 <- aggregate(FilteredInteractions$Count, by=list(FilteredInteractions$oeID), FUN="sum")
colnames(Table1) <- c("oeID","Count")

# Define a variable for induced interactions (log2(FC) >= log2(1.5)), summarize and set column names
FilteredInteractions$Response <- 0
FilteredInteractions[ FilteredInteractions$logFC_Four_Zero >= log2(1.5),"Response"] <- 1
Table2 <- aggregate(FilteredInteractions$Response, by=list(FilteredInteractions$oeID), FUN="sum")
colnames(Table2) <- c("oeID","Induced")

# Define a variable for repressed interactions (log2(FC) <= -log2(1.5)), summarize and set column names
FilteredInteractions$Response <- 0
FilteredInteractions[ FilteredInteractions$logFC_Four_Zero <= -log2(1.5),"Response"] <- 1
Table3 <- aggregate(FilteredInteractions$Response, by=list(FilteredInteractions$oeID), FUN="sum")
colnames(Table3) <- c("oeID","Repressed")

# Merge the tables
Table <- merge(Table1, Table2, by ="oeID")
Table <- merge(Table, Table3, by ="oeID")

# Setup to define response on overlapped interactions
Overlap$Response <- 0

# Set all of the interactions from otherEnds where the total number of interactions is equal to the number of dynamic interactions to 1
Overlap[ Overlap$oeID %in% Table[ Table$Count > 0 & Table$Induced > 0 & Table$Repressed == 0,"oeID"],"Response"] <- 1

# Set all of the interactions from otherEnds where the total number of interactions is equal to -1 multiplied by the number of dynamic interactions to -1
Overlap[ Overlap$oeID %in% Table[ Table$Count > 0 & Table$Induced == 0 & Table$Repressed > 0,"oeID"],"Response"] <- -1

# Set the response of all interactions with a log2 fold change less than 1.5 to zero
Overlap[ abs(Overlap$logFC_Four_Zero) < log2(1.5),"Response"] <- 0

# Get only peaks bound by any factor in two replicates
TestPeaks <- Peaks
TestPeaks2 <- data.frame()
for (i in c(1:9)) {
	Factor <- colnames(Counts)[(i+1)]
	TestPeaks$Zero_1 <- regexpr(paste(as.character(Factor),"_D0_E1",sep=""), TestPeaks$Factors)
	TestPeaks$Zero_2 <- regexpr(paste(as.character(Factor),"_D0_E2",sep=""), TestPeaks$Factors)
	TestPeaks$Four_1 <- regexpr(paste(as.character(Factor),"_4h_E1",sep=""), TestPeaks$Factors)
	TestPeaks$Four_2 <- regexpr(paste(as.character(Factor),"_4h_E2",sep=""), TestPeaks$Factors)
	Tmp <- TestPeaks[ (TestPeaks$Zero_1 > -1 & TestPeaks$Zero_2 > -1) | (TestPeaks$Four_1 > -1 & TestPeaks$Four_2 > -1),]
	TestPeaks2 <- rbind(TestPeaks2, Tmp)
	}
TestPeaks2 <- TestPeaks2[ duplicated(TestPeaks2$PeakID) == F,]

# Get interactions with peaks
Overlap2 <- Overlap[ Overlap$PeakID %in% TestPeaks2$PeakID,]

# For each factor calculate the fractions
Direction <- data.frame(matrix(ncol=4,nrow=9))
for (i in c(1:9)) {
	Factor <- colnames(Counts)[(i+1)]
	TestPeaks$Zero_1 <- regexpr(paste(as.character(Factor),"_D0_E1",sep=""), TestPeaks$Factors)
	TestPeaks$Zero_2 <- regexpr(paste(as.character(Factor),"_D0_E2",sep=""), TestPeaks$Factors)
	TestPeaks$Four_1 <- regexpr(paste(as.character(Factor),"_4h_E1",sep=""), TestPeaks$Factors)
	TestPeaks$Four_2 <- regexpr(paste(as.character(Factor),"_4h_E2",sep=""), TestPeaks$Factors)
	Tmp <- TestPeaks[ (TestPeaks$Zero_1 > -1 & TestPeaks$Zero_2 > -1) | (TestPeaks$Four_1 > -1 & TestPeaks$Four_2 > -1),]
	Tmp <- merge(Tmp,Counts[,c(1,i+1)], by="PeakID")
	Tmp <- merge(Tmp[,c(1,10)], Overlap2[,c(2,4,5)], by="PeakID")
	Tmp <- Tmp[ Tmp$Response != 0,]
	Direction[i,1] <- Factor
	Direction[i,2] <- nrow(Tmp[ (Tmp[,2] >= log2(1.5) & Tmp[,4] == 1) | (Tmp[,2] <= -log2(1.5) & Tmp[,4] == -1),])/nrow(Tmp)
	Direction[i,3] <- nrow(Tmp[ (Tmp[,2] >= log2(1.5) & Tmp[,4] == -1) | (Tmp[,2] <= -log2(1.5) & Tmp[,4] == 1),])/nrow(Tmp)
	Direction[i,4] <- nrow(Tmp[ abs(Tmp[,2]) < log2(1.5),])/nrow(Tmp)
	}

# Plot it
barplot(t(as.matrix(Direction[c(7,6,5,4,1,3,2,8,9),c(2,3,4)])), las=2, col=c("darkgreen","darkred","grey"), names=Direction[c(7,6,5,4,1,3,2,8,9),1])
```

[Back to start](../README.md)<br>
[Back to overview of Figure 6](../Links/Figure6.md)