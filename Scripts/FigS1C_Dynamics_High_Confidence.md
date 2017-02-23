```R
## Load the necessary packages
library(edgeR)
library(GenomicRanges)

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

## Integrate the data
# Keep only interactions overlapping an enhancer
Combined <- Combined[ Combined$oeID %in% Overlap$oeID,]

# Count the number of interactions in different ranges of maximal log2 fold change and place the results into a new matrix
Counts <- matrix(ncol=1, nrow=7)
Counts[1] <- nrow(Combined[ Combined$logFC_max <= -log2(3),])
Counts[2] <- nrow(Combined[ Combined$logFC_max > -log2(3) & Combined$logFC_max <= -log2(2),])
Counts[3] <- nrow(Combined[ Combined$logFC_max > -log2(2) & Combined$logFC_max <= -log2(1.5),])
Counts[4] <- nrow(Combined[ Combined$logFC_max > -log2(1.5) & Combined$logFC_max < log2(1.5),])
Counts[5] <- nrow(Combined[ Combined$logFC_max >= log2(1.5) & Combined$logFC_max < log2(2),])
Counts[6] <- nrow(Combined[ Combined$logFC_max >= log2(2) & Combined$logFC_max < log2(3),])
Counts[7] <- nrow(Combined[ Combined$logFC_max >= log2(3),])

# Calculate the fraction of interactions in each ranges
Counts <- Counts / nrow(Combined)

# Plot a barplot
barplot(Counts, beside = T, las = 1, names=c("(-Inf,-3]","(-3,-2]","(-2,-1.5]","(-1.5,-1.5)","[1.5,2)","[2,3)","[3,Inf)"),
	col=c("#8B0000", "#C40000", "#FF0000","#DCDCDC","#00FF00", "#00B100", "#006400"), ylab="Fraction of interactions")
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)
