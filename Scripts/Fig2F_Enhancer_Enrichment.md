```R
## Load the necessary packages
library(DESeq2)
library(edgeR)
library(org.Mm.eg.db)
library(GenomicRanges)
library(rtracklayer)

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

# Merge gene information and collaps to keep only unique interaction-gene pairs
FilteredInteractions_Genes <- merge(FilteredInteractions, Genes[,c("Fragment","Symbol")], by.x="baitID", by.y="Fragment")
FilteredInteractions_Genes$Collaps <- paste(FilteredInteractions_Genes$ID, FilteredInteractions_Genes$Symbol, sep="-")
FilteredInteractions_Genes <- FilteredInteractions_Genes[ duplicated(FilteredInteractions_Genes$Collaps)==F,]
FilteredInteractions_Genes <- FilteredInteractions_Genes[,c(1:36)]

## Analyze enhancer data
Peaks <- read.delim("Data/ChIPseq/Peak_lists/Activators.bed")
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

# Save the result
Fraction <- nrow(Overlap[ duplicated(Overlap$ID)==F,])/nrow(FilteredInteractions[ duplicated(FilteredInteractions$ID)==F,])

## Shuffle the interactions and overlap 1000 times
FilteredInteractions$Size <- FilteredInteractions$oeEnd - FilteredInteractions$oeStart
FilteredInteractions$baitCenter <- ((FilteredInteractions$baitEnd - FilteredInteractions$baitStart)/2)+FilteredInteractions$baitStart
SummaryData <- data.frame(matrix(ncol=1,nrow=1000))
EnhancerRanges <- GRanges(seqnames = Peaks$Chr, IRanges(start = Peaks$Start, end = Peaks$End), strand = rep("+",nrow(Peaks)), mcols = Peaks$PeakID)

for (i in 1:1000) {
  BaitList <- data.frame(baitID = FilteredInteractions[ duplicated(FilteredInteractions$baitID)==F,"baitID"])
  BaitList$Random <- sample(BaitList$baitID, size = nrow(BaitList), replace = T)
  BaitList <- merge(BaitList, FilteredInteractions[,c("dist","Size","baitID","baitCenter","baitChr","ID")], by.x="Random", by.y="baitID")
  BaitList$Mid <- BaitList$baitCenter + BaitList$dist
  BaitList$oeStart <- BaitList$Mid - (BaitList$Size/2)
  BaitList$oeEnd <- BaitList$Mid + (BaitList$Size/2)
  colnames(BaitList)[6] <- "oeChr"
  BaitList$oeID <- seq(1,nrow(BaitList), by=1)
  otherEndRanges <- GRanges(seqnames = BaitList$oeChr, IRanges(start = BaitList$oeStart, end = BaitList$oeEnd), strand = rep("+",nrow(BaitList)), mcols = BaitList$oeID)
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
  Overlap <- merge(Overlap, BaitList[,c("oeID","ID")], by="oeID")
  SummaryData[i,1] <- nrow(Overlap[ duplicated(Overlap$ID)==F,])/nrow(FilteredInteractions[ duplicated(FilteredInteractions$ID)==F,])
}

# Plot the results
barCenters <- barplot(c(Fraction,mean(SummaryData[,1])), las=1, ylim=c(0,0.5), names = c("PIR","Randomized"), ylab = "Overlap with enhancers", col=c("green3","grey"))
arrows(barCenters[2], mean(SummaryData[,1]) - sd(SummaryData[,1]) * 2, barCenters[2],mean(SummaryData[,1]) + sd(SummaryData[,1]) * 2, lwd = 1.5, angle = 90, code = 3, length = 0.05)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 2](../Links/Figure2.md)
