```R
## Load the necessary packages
library(GenomicRanges)
library(edgeR)

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

## Process TADs
# Import the data
Zero_TAD <- read.delim("Data/Interactions/TADs/D0_40kb.domains.bed", header=F, skip=1)
Four_TAD <- read.delim("Data/Interactions/TADs/4h_40kb.domains.bed", header=F, skip=1)
Two_TAD <- read.delim("Data/Interactions/TADs/D2_40kb.domains.bed", header=F, skip=1)

# Set column names
colnames(Zero_TAD) <- c("Chr","Start","End","PeakID","Score","Strand")
colnames(Four_TAD) <- c("Chr","Start","End","PeakID","Score","Strand")
colnames(Two_TAD) <- c("Chr","Start","End","PeakID","Score","Strand")

# Make new unique numeric PeakIDs
Zero_TAD$PeakID <- seq(1, nrow(Zero_TAD),by=1)
Four_TAD$PeakID <- seq(1, nrow(Four_TAD),by=1)
Two_TAD$PeakID <- seq(1, nrow(Two_TAD),by=1)

## Overlap the interactions with D0 TADs
# otherEnds
TADRanges <- GRanges(seqnames = Zero_TAD$Chr, IRanges(start = Zero_TAD$Start, end = Zero_TAD$End), strand = rep("+",nrow(Zero_TAD)), mcols = Zero_TAD$PeakID)
otherEnds <- FilteredInteractions[ duplicated(FilteredInteractions$oeID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$oeChr, IRanges(start = otherEnds$oeStart, end = otherEnds$oeEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$oeID)

# Get the overlap and width of overlap
Overlap <- suppressWarnings(findOverlaps(otherEndRanges, TADRanges))
Width <- as.data.frame(ranges(Overlap, ranges(otherEndRanges), ranges(TADRanges)))
Overlap <- as.data.frame(Overlap)

# Convert to data frames
otherEnds <- as.data.frame(otherEndRanges[Overlap[,1]])
TAD <- as.data.frame(TADRanges[Overlap[,2]])

# Set column names
colnames(otherEnds)[6] <- "oeID"
colnames(TAD)[6] <- "PeakID"

# Make an ID for merging
otherEnds$Merge <- seq(1, nrow(otherEnds), by=1)
TAD$Merge <- seq(1, nrow(TAD), by=1)
Width$Merge <- seq(1, nrow(Width), by=1)

# Merge the overlap results
Overlap <- merge(otherEnds[,c(6,7)], TAD[,c(6,7)])
Overlap <- merge(Overlap, Width[,c(3,4)], by="Merge")

# Collaps overlap to only the largest overlap
Overlap <- Overlap[ order(Overlap$oeID, -Overlap$width),]
Overlap <- Overlap[ duplicated(Overlap$oeID)==F,]
Overlap <- Overlap[,c(2,3)]
colnames(Overlap) <- c("oeID","TAD_Zero_oe")

# Add the information to the interaction data frame
FilteredInteractions <- merge(FilteredInteractions, Overlap, by="oeID", all.x=T)

## Baits
otherEnds <- FilteredInteractions[ duplicated(FilteredInteractions$baitID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$baitChr, IRanges(start = otherEnds$baitStart, end = otherEnds$baitEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$baitID)

# Get the overlap and width of overlap
Overlap <- suppressWarnings(findOverlaps(otherEndRanges, TADRanges))
Width <- as.data.frame(ranges(Overlap, ranges(otherEndRanges), ranges(TADRanges)))
Overlap <- as.data.frame(Overlap)

# Convert to data frames
otherEnds <- as.data.frame(otherEndRanges[Overlap[,1]])
TAD <- as.data.frame(TADRanges[Overlap[,2]])

# Set column names
colnames(otherEnds)[6] <- "oeID"
colnames(TAD)[6] <- "PeakID"

# Make an ID for merging
otherEnds$Merge <- seq(1, nrow(otherEnds), by=1)
TAD$Merge <- seq(1, nrow(TAD), by=1)
Width$Merge <- seq(1, nrow(Width), by=1)

# Merge the overlap results
Overlap <- merge(otherEnds[,c(6,7)], TAD[,c(6,7)])
Overlap <- merge(Overlap, Width[,c(3,4)], by="Merge")

# Collaps overlap to only the largest overlap
Overlap <- Overlap[ order(Overlap$oeID, -Overlap$width),]
Overlap <- Overlap[ duplicated(Overlap$oeID)==F,]
Overlap <- Overlap[,c(2,3)]
colnames(Overlap) <- c("baitID","TAD_Zero_bait")

# Add the information to the interaction data frame
FilteredInteractions <- merge(FilteredInteractions, Overlap, by="baitID", all.x=T)

## Overlap the interactions with 4h TADs
## otherEnds
TADRanges <- GRanges(seqnames = Four_TAD$Chr, IRanges(start = Four_TAD$Start, end = Four_TAD$End), strand = rep("+",nrow(Four_TAD)), mcols = Four_TAD$PeakID)
otherEnds <- FilteredInteractions[ duplicated(FilteredInteractions$oeID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$oeChr, IRanges(start = otherEnds$oeStart, end = otherEnds$oeEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$oeID)

# Get the overlap and width of overlap
Overlap <- suppressWarnings(findOverlaps(otherEndRanges, TADRanges))
Width <- as.data.frame(ranges(Overlap, ranges(otherEndRanges), ranges(TADRanges)))
Overlap <- as.data.frame(Overlap)

# Convert to data frames
otherEnds <- as.data.frame(otherEndRanges[Overlap[,1]])
TAD <- as.data.frame(TADRanges[Overlap[,2]])

# Set column names
colnames(otherEnds)[6] <- "oeID"
colnames(TAD)[6] <- "PeakID"

# Make an ID for merging
otherEnds$Merge <- seq(1, nrow(otherEnds), by=1)
TAD$Merge <- seq(1, nrow(TAD), by=1)
Width$Merge <- seq(1, nrow(Width), by=1)

# Merge the overlap results
Overlap <- merge(otherEnds[,c(6,7)], TAD[,c(6,7)])
Overlap <- merge(Overlap, Width[,c(3,4)], by="Merge")

# Collaps overlap to only the largest overlap
Overlap <- Overlap[ order(Overlap$oeID, -Overlap$width),]
Overlap <- Overlap[ duplicated(Overlap$oeID)==F,]
Overlap <- Overlap[,c(2,3)]
colnames(Overlap) <- c("oeID","TAD_Four_oe")

# Add the information to the interaction data frame
FilteredInteractions <- merge(FilteredInteractions, Overlap, by="oeID", all.x=T)

## Baits
otherEnds <- FilteredInteractions[ duplicated(FilteredInteractions$baitID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$baitChr, IRanges(start = otherEnds$baitStart, end = otherEnds$baitEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$baitID)

# Get the overlap and width of overlap
Overlap <- suppressWarnings(findOverlaps(otherEndRanges, TADRanges))
Width <- as.data.frame(ranges(Overlap, ranges(otherEndRanges), ranges(TADRanges)))
Overlap <- as.data.frame(Overlap)

# Convert to data frames
otherEnds <- as.data.frame(otherEndRanges[Overlap[,1]])
TAD <- as.data.frame(TADRanges[Overlap[,2]])

# Set column names
colnames(otherEnds)[6] <- "oeID"
colnames(TAD)[6] <- "PeakID"

# Make an ID for merging
otherEnds$Merge <- seq(1, nrow(otherEnds), by=1)
TAD$Merge <- seq(1, nrow(TAD), by=1)
Width$Merge <- seq(1, nrow(Width), by=1)

# Merge the overlap results
Overlap <- merge(otherEnds[,c(6,7)], TAD[,c(6,7)])
Overlap <- merge(Overlap, Width[,c(3,4)], by="Merge")

# Collaps overlap to only the largest overlap
Overlap <- Overlap[ order(Overlap$oeID, -Overlap$width),]
Overlap <- Overlap[ duplicated(Overlap$oeID)==F,]
Overlap <- Overlap[,c(2,3)]
colnames(Overlap) <- c("baitID","TAD_Four_bait")

# Add the information to the interaction data frame
FilteredInteractions <- merge(FilteredInteractions, Overlap, by="baitID", all.x=T)

## Overlap the interactions with D2 TADs
# otherEnds
TADRanges <- GRanges(seqnames = Two_TAD$Chr, IRanges(start = Two_TAD$Start, end = Two_TAD$End), strand = rep("+",nrow(Two_TAD)), mcols = Two_TAD$PeakID)
otherEnds <- FilteredInteractions[ duplicated(FilteredInteractions$oeID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$oeChr, IRanges(start = otherEnds$oeStart, end = otherEnds$oeEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$oeID)

# Get the overlap and width of overlap
Overlap <- suppressWarnings(findOverlaps(otherEndRanges, TADRanges))
Width <- as.data.frame(ranges(Overlap, ranges(otherEndRanges), ranges(TADRanges)))
Overlap <- as.data.frame(Overlap)

# Convert to data frames
otherEnds <- as.data.frame(otherEndRanges[Overlap[,1]])
TAD <- as.data.frame(TADRanges[Overlap[,2]])

# Set column names
colnames(otherEnds)[6] <- "oeID"
colnames(TAD)[6] <- "PeakID"

# Make an ID for merging
otherEnds$Merge <- seq(1, nrow(otherEnds), by=1)
TAD$Merge <- seq(1, nrow(TAD), by=1)
Width$Merge <- seq(1, nrow(Width), by=1)

# Merge the overlap results
Overlap <- merge(otherEnds[,c(6,7)], TAD[,c(6,7)])
Overlap <- merge(Overlap, Width[,c(3,4)], by="Merge")

# Collaps overlap to only the largest overlap
Overlap <- Overlap[ order(Overlap$oeID, -Overlap$width),]
Overlap <- Overlap[ duplicated(Overlap$oeID)==F,]
Overlap <- Overlap[,c(2,3)]
colnames(Overlap) <- c("oeID","TAD_Two_oe")

# Add the information to the interaction data frame
FilteredInteractions <- merge(FilteredInteractions, Overlap, by="oeID", all.x=T)

## Baits
otherEnds <- FilteredInteractions[ duplicated(FilteredInteractions$baitID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$baitChr, IRanges(start = otherEnds$baitStart, end = otherEnds$baitEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$baitID)

# Get the overlap and width of overlap
Overlap <- suppressWarnings(findOverlaps(otherEndRanges, TADRanges))
Width <- as.data.frame(ranges(Overlap, ranges(otherEndRanges), ranges(TADRanges)))
Overlap <- as.data.frame(Overlap)

# Convert to data frames
otherEnds <- as.data.frame(otherEndRanges[Overlap[,1]])
TAD <- as.data.frame(TADRanges[Overlap[,2]])

# Set column names
colnames(otherEnds)[6] <- "oeID"
colnames(TAD)[6] <- "PeakID"

# Make an ID for merging
otherEnds$Merge <- seq(1, nrow(otherEnds), by=1)
TAD$Merge <- seq(1, nrow(TAD), by=1)
Width$Merge <- seq(1, nrow(Width), by=1)

# Merge the overlap results
Overlap <- merge(otherEnds[,c(6,7)], TAD[,c(6,7)])
Overlap <- merge(Overlap, Width[,c(3,4)], by="Merge")

# Collaps overlap to only the largest overlap
Overlap <- Overlap[ order(Overlap$oeID, -Overlap$width),]
Overlap <- Overlap[ duplicated(Overlap$oeID)==F,]
Overlap <- Overlap[,c(2,3)]
colnames(Overlap) <- c("baitID","TAD_Two_bait")

# Add the information to the interaction data frame
FilteredInteractions <- merge(FilteredInteractions, Overlap, by="baitID", all.x=T)

## Calculate the fraction of interactions with both ends within any TAD
# Get interactions with any anchor inside any TAD at each end and timepoint
D0_Bait	<- FilteredInteractions[ !is.na(FilteredInteractions$TAD_Zero_bait),]
D0_otherEnd <- FilteredInteractions[ !is.na(FilteredInteractions$TAD_Zero_oe),]
H4_Bait	<- FilteredInteractions[ !is.na(FilteredInteractions$TAD_Four_bait),]
H4_otherEnd <- FilteredInteractions[ !is.na(FilteredInteractions$TAD_Four_oe),]
D2_Bait	<- FilteredInteractions[ !is.na(FilteredInteractions$TAD_Two_bait),]
D2_otherEnd <- FilteredInteractions[ !is.na(FilteredInteractions$TAD_Two_oe),]

# Combine the results and filter the interactions
D0 <- data.frame(ID = c(as.character(D0_Bait$ID),as.character(D0_otherEnd$ID)))
H4 <- data.frame(ID = c(as.character(H4_Bait$ID),as.character(H4_otherEnd$ID)))
D2 <- data.frame(ID = c(as.character(D2_Bait$ID),as.character(D2_otherEnd$ID)))
D0_AnchoredTAD <- FilteredInteractions[ FilteredInteractions$ID %in% D0$ID,]
H4_AnchoredTAD <- FilteredInteractions[ FilteredInteractions$ID %in% H4$ID,]
D2_AnchoredTAD <- FilteredInteractions[ FilteredInteractions$ID %in% D2$ID,]

# Combine the IDs of interactions anchored TADs at any time point
Anchored <- data.frame(ID = c(as.character(D0_AnchoredTAD$ID),as.character(H4_AnchoredTAD$ID),as.character(D2_AnchoredTAD$ID)))

# Get the interactions within TADs at each time point
D0_WithinTAD <- D0_AnchoredTAD[ D0_AnchoredTAD$TAD_Zero_bait == D0_AnchoredTAD$TAD_Zero_oe,]
H4_WithinTAD <- H4_AnchoredTAD[ H4_AnchoredTAD$TAD_Four_bait == H4_AnchoredTAD$TAD_Four_oe,]
D2_WithinTAD <- D2_AnchoredTAD[ D2_AnchoredTAD$TAD_Two_bait == D2_AnchoredTAD$TAD_Two_oe,]

# Get rid of all NA's
D0_WithinTAD <- D0_WithinTAD[ !is.na(D0_WithinTAD$TAD_Zero_oe),]
H4_WithinTAD <- H4_WithinTAD[ !is.na(H4_WithinTAD$TAD_Zero_oe),]
D2_WithinTAD <- D0_WithinTAD[ !is.na(D2_WithinTAD$TAD_Zero_oe),]

# Combine the IDs of interactions within TADs at any time point
Within <- data.frame(ID = c(as.character(D0_WithinTAD$ID),as.character(H4_WithinTAD$ID),as.character(D2_WithinTAD$ID)))

# Make a matrix to capture the results
Result <- matrix(ncol=2, nrow=2)

# Calculate the proportion of interactions within TADs
Result[1,1] <- nrow(FilteredInteractions[ !(FilteredInteractions$ID %in% Anchored$ID),])/nrow(FilteredInteractions)
Result[2,1] <- 0
Result[1,2] <- nrow(FilteredInteractions[ FilteredInteractions$ID %in% Anchored$ID & FilteredInteractions$ID %in% Within$ID,])/nrow(FilteredInteractions)
Result[2,2] <- nrow(FilteredInteractions[ FilteredInteractions$ID %in% Anchored$ID & !(FilteredInteractions$ID %in% Within$ID),])/nrow(FilteredInteractions)

# Plot the results
barplot(Result, ylim=c(0,1), las=1, ylab="Fraction of interactions", col=c("grey","blue3"), names=c("Not anchored","Anchored"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)