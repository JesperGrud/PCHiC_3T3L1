```R
## Load the necessary packages
library(GenomicRanges)
library(edgeR)

## Identify TADs and boundary regions
# Import the data
Zero <- read.table("Data/Interactions/TADs/D0.DI", quote="\"", comment.char="", skip=1)

## Call TADs on sDI
# Convert sDI to z-scores
Mean <- mean(Zero$V4)
SD <- sd(Zero$V4)
Zero$V4 <- (Zero$V4 - Mean)/SD
rm(Mean)
rm(SD)

All <- data.frame()
Bounds <- data.frame()
Chromosomes <- c(paste("chr",seq(1,19,by=1),sep=""),"chrX")
for (Chr in Chromosomes) {
  # Subset to a single chromosome
  Subset <- Zero[ Zero$V1 == Chr,]
  
  # Identify transition points (between positive and negative bias)
  Subset$Transition <- 0
  for (i in 1:(nrow(Subset)-1)) {
    SignCurr <- sign(Subset[i,4])
    SignNext <- sign(Subset[(i+1),4])
    if (SignCurr != SignNext) { Subset[c(i),"Transition"] <- 1}
  }
  
  # Identify each region of bias
  Region <- 1
  Subset$Region <- 0
  for (i in 1:nrow(Subset)) {
    Subset[i,"Region"] <- Region
    if (Subset[i,"Transition"] == 1) { Region <- Region + 1}
  }
  
  # Calculate the length of each region
  Subset$Count <- 1
  Counts <- aggregate(Subset$Count, by=list(Subset$Region), FUN="sum")
  colnames(Counts) <- c("Region","Length")
  Subset <- merge(Subset, Counts, by="Region")
  
  # Identify regions longer or equal to 10 bins
  Subset$Long <- 0
  Subset[ Subset$Length >= 10,"Long"] <- 1
  
  # Mark each region for sufficient bias and set the sign
  Subset$Bias <- 0
  Subset$Sign <- 0
  for (i in unique(Subset$Region)) {
    Tmp <- Subset[ Subset$Region == i,]
    if (max(abs(Tmp$V4)) >= 0.5) {
      Subset[ Subset$Region == i, "Bias"] <- 1
      Subset[ Subset$Region == i, "Sign"] <- sign(max(Tmp$V4))
    }
  }
  
  # Remove regions that are too short (< 10 bins = 50kb)
  Subset[ Subset$Length <= 10, "Bias"] <- 0
  
  # Identify biased regions
  PositiveRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == 1, "Region"])
  NegativeRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == -1, "Region"])
  
  ## Identify type 1 boundaries (left > right)
  # For each negative region, get the closest positive region
  T1 <- data.frame()
  for (i in 1:length(NegativeRegions)) {
    if (length(PositiveRegions[ PositiveRegions >= NegativeRegions[i]]) > 0) {
      T1[i,1] <- NegativeRegions[i]
      T1[i,2] <- min(PositiveRegions[ PositiveRegions >= NegativeRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T1 <- T1[ order(T1[,2], -T1[,1]),]
  T1 <- T1[ duplicated(T1[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T1$Keep <- 0
  T1[ T1[,2] == (T1[,1] + 1),"Keep"] <- 1
  Proceeding <- T1[ T1$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT1 <- T1[ T1$Keep == 0,]
  for (i in 1:nrow(GappedT1)) {
    Test <- Subset[ Subset$Region < GappedT1[i,2] & Subset$Region > GappedT1[i,1],]
    if (nrow(Test) <= 8) { GappedT1[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT1[i,"Keep"] <- 0 }
    if (GappedT1[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT1[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT1[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT1[ GappedT1$Keep == 1,]) > 0) {
  GappedT1Keep <- GappedT1[ GappedT1$Keep == 1, c(4,5,3)]
  colnames(GappedT1Keep) <- colnames(Proceeding)
  T1 <- rbind(GappedT1Keep, Proceeding)
  }
  T1 <- T1[ order(T1[,1]),]
  
  # Identify type 2 boundaries (right > left)
  T2 <- data.frame()
  for (i in 1:length(PositiveRegions)) {
    if (length(NegativeRegions[ NegativeRegions >= PositiveRegions[i]]) > 0) {
      T2[i,1] <- PositiveRegions[i]
      T2[i,2] <- min(NegativeRegions[ NegativeRegions >= PositiveRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T2 <- T2[ order(T2[,2], -T2[,1]),]
  T2 <- T2[ duplicated(T2[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T2$Keep <- 0
  T2[ T2[,2] == (T2[,1] + 1),"Keep"] <- 1
  Proceeding <- T2[ T2$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT2 <- T2[ T2$Keep == 0,]
  for (i in 1:nrow(GappedT2)) {
    Test <- Subset[ Subset$Region < GappedT2[i,2] & Subset$Region > GappedT2[i,1],]
    if (nrow(Test) <= 8) { GappedT2[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT2[i,"Keep"] <- 0 }
    if (GappedT2[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT2[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT2[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT2[ GappedT2$Keep == 1,]) > 0) {
    GappedT2Keep <- GappedT2[ GappedT2$Keep == 1, c(4,5,3)]
    colnames(GappedT2Keep) <- colnames(Proceeding)
    T2 <- rbind(GappedT2Keep, Proceeding)
  }
  T2 <- T2[ order(T2[,1]),]
  
  # Get non-boundary biased regions
  PositiveNonbias <- data.frame(Region = PositiveRegions[ !(PositiveRegions %in% T1[,2]) & !(PositiveRegions %in% GappedT1[ GappedT1$Keep == 1,2])])
  NegativeNonbias <- data.frame(Region = NegativeRegions[ !(NegativeRegions %in% T2[,2]) & !(NegativeRegions %in% GappedT2[ GappedT2$Keep == 1,2])])
  
  # Stich together T1 regions
  TADs <- data.frame()
  for (i in 1:(nrow(T1)-1)) {
    TADs[i,1] <- T1[i,2]
    TADs[i,2] <- T1[(i+1),2]
  }
  
  # Remove TADs that doesnt cross a T2 boundary or multiple
  TADs$Keep <- 0
  TADs$T2 <- 0
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,2]
    Tmp <- data.frame(T2 = T2[ T2[,2] >= Min & T2[,2] <= Max,])
    if (nrow(Tmp) == 1) { 
      TADs[i,"Keep"] <- 1
      TADs[i,"T2"] <- Tmp[1,1]
    }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have negative bias between first T1 and T2 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,4]
    Tmp <- data.frame(NegativeNonbias = NegativeNonbias[ NegativeNonbias > Min & NegativeNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have positive bias between T2 and second T1 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,4]
    Max <- TADs[i,2]
    Tmp <- data.frame(PositiveNonbias = PositiveNonbias[ PositiveNonbias > Min & PositiveNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate TAD boundaries
  TADs$chr <- as.character(Chr)
  for (i in 1:nrow(TADs)) {
    Tmp <- Subset[ Subset$Region == TADs[i,1],]
    TADs[i,"start"] <- Tmp[ which.max(Tmp[,5]),3]
    Tmp <- Subset[ Subset$Region == (TADs[i,2]-1),]
    TADs[i,"end"] <- Tmp[ which.min(Tmp[,5]),3]+50000
  }
  
  # Filter TADs in low coverage areas (+/- 20 bins)
  Subset$Forward <- Subset[2,3] - Subset[1,3]
  for (nr in 2:nrow(Subset)) { Subset[nr,"Forward"] <- Subset[nr,3] - Subset[(nr-1),3] }
  
  for (m in 1:nrow(TADs)) {
    Start <- as.numeric(TADs[m,6])-50000
    End <- as.numeric(TADs[m,7])+50000
    Tmp <- Subset
    Tmp <- Tmp[ Tmp[,3] >= Start & Tmp[,3] < End,]
    if (max(Tmp$Forward) > 5000) { TADs[m,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate T1 boundaries
  T1$chr <- as.character(Chr)
  for (i in 1:nrow(T1)) {
    Tmp <- Subset[ Subset$Region == T1[i,2],]
    T1[i,"boundary"] <- min(Tmp[,3])
  }
  
  # Finalize TADs
  TADs <- TADs[,c(5,6,7)]
  
  # Combine all TADs
  All <- rbind(All, TADs)
  
  # Save T1 Boundaries
  T1 <- T1[,c(4,5)]
  Bounds <- rbind(Bounds, T1)
}

# Rename the TAD calls
Zero_TADs <- All
Zero_Bounds <- Bounds

## Clean up
rm(list=setdiff(ls(),c("Zero_TADs")))

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

# Keep only D0 interactions
FilteredInteractions <- FilteredInteractions[ FilteredInteractions$Zero_Score >= 5,]

## Integrate interactions and TADs
# Make IDs
Zero_TADs$TADID <- seq(1, nrow(Zero_TADs), by=1)

## Overlap the interactions with D0 TADs
# otherEnds
TADRanges <- GRanges(seqnames = Zero_TADs$chr, IRanges(start = Zero_TADs$start, end = Zero_TADs$end), strand = rep("+",nrow(Zero_TADs)), mcols = Zero_TADs$TADID)
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
colnames(TAD)[6] <- "TADID"

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
colnames(TAD)[6] <- "TADID"

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

## Filter interactions
FilteredInteractions$TAD_Zero_oe <- as.character(FilteredInteractions$TAD_Zero_oe)
FilteredInteractions$TAD_Zero_bait <- as.character(FilteredInteractions$TAD_Zero_bait)

# Keep only interactions where the bait is inside a TAD
FilteredInteractions <- FilteredInteractions[ !is.na(FilteredInteractions$TAD_Zero_bait),]

# IntraTAD
Within <- FilteredInteractions[ !is.na(FilteredInteractions$TAD_Zero_oe) & !is.na(FilteredInteractions$TAD_Zero_bait),]
Within <- Within[ Within$TAD_Zero_oe == Within$TAD_Zero_bait,]
A <- nrow(Within)/nrow(FilteredInteractions)

# TransTAD
Outside <- FilteredInteractions[ !(FilteredInteractions$ID %in% Within$ID),]
B <- nrow(Outside)/nrow(FilteredInteractions)

## Process randomized control
# Get the relative position of the baitCenter
FilteredInteractions$baitCenter <- (FilteredInteractions$baitEnd + FilteredInteractions$baitStart)/2
FilteredInteractions <- merge(FilteredInteractions, Zero_TADs[,c(2,3,4)], by.x="TAD_Zero_bait",by.y="TADID")
FilteredInteractions$Size <- FilteredInteractions$end - FilteredInteractions$start
FilteredInteractions$Distance <- FilteredInteractions$baitCenter - FilteredInteractions$start
FilteredInteractions$RelativeDistance <- FilteredInteractions$Distance / FilteredInteractions$Size

# Calculate sizes for reestablishing the interaction network
FilteredInteractions$baitSize <- (FilteredInteractions$baitEnd - FilteredInteractions$baitStart)
FilteredInteractions$oeSize <- (FilteredInteractions$oeEnd - FilteredInteractions$oeStart)

# Collapsed by TADs
TADlist <- FilteredInteractions[ duplicated(FilteredInteractions$TAD_Zero_bait)==F, ]

# Make a data.frame to capture the results
Result <- data.frame(matrix(ncol=2, nrow= 100))

# Shuffle the TAD for each bait
for (i in 1:100) {
BaitList <- data.frame(baitID = FilteredInteractions[ duplicated(FilteredInteractions$baitID)==F, "baitID"])
BaitList$Shuffled <- sample(TADlist$TAD_Zero_bait, nrow(BaitList), replace=T)
BaitList <- merge(BaitList, TADlist[,c("TAD_Zero_bait","start","Size","baitChr")], by.y="TAD_Zero_bait", by.x="Shuffled")
BaitList <- BaitList[ duplicated(BaitList$baitID)==F,]
BaitList <- merge(BaitList, FilteredInteractions[,c("RelativeDistance","baitSize","baitID")], by="baitID")
BaitList <- BaitList[ duplicated(BaitList$baitID)==F,]
BaitList$baitCenter <- (BaitList$RelativeDistance * BaitList$Size) + BaitList$start
BaitList$baitStart <- BaitList$baitCenter - (BaitList$baitSize / 2)
BaitList$baitEnd <- BaitList$baitCenter + (BaitList$baitSize / 2)
BaitList$oeChr <- BaitList$baitChr
BaitList <- merge(BaitList, FilteredInteractions[,c("baitID","ID","dist","oeSize")], by="baitID")
BaitList$oeMid <- BaitList$baitCenter + BaitList$dist
BaitList$oeStart <- BaitList$oeMid - (BaitList$oeSize / 2)
BaitList$oeEnd <- BaitList$oeMid + (BaitList$oeSize / 2)

# Overlap otherEnds
otherEndRanges <- GRanges(seqnames = BaitList$oeChr, IRanges(start = BaitList$oeStart, end = BaitList$oeEnd), strand = rep("+",nrow(BaitList)), mcols = BaitList$ID)
Overlap <- suppressWarnings(findOverlaps(otherEndRanges, TADRanges))
Width <- as.data.frame(ranges(Overlap, ranges(otherEndRanges), ranges(TADRanges)))
Overlap <- as.data.frame(Overlap)

# Convert to data frames
otherEnds <- as.data.frame(otherEndRanges[Overlap[,1]])
TAD <- as.data.frame(TADRanges[Overlap[,2]])

# Set column names
colnames(otherEnds)[6] <- "ID"
colnames(TAD)[6] <- "TADID"

# Make an ID for merging
otherEnds$Merge <- seq(1, nrow(otherEnds), by=1)
TAD$Merge <- seq(1, nrow(TAD), by=1)
Width$Merge <- seq(1, nrow(Width), by=1)

# Merge the overlap results
Overlap <- merge(otherEnds[,c(6,7)], TAD[,c(6,7)], by="Merge")
Overlap <- merge(Overlap, Width[,c(3,4)], by="Merge")

# Collaps overlap to only the largest overlap
Overlap <- Overlap[ order(Overlap$ID, -Overlap$width),]
Overlap <- Overlap[ duplicated(Overlap$ID)==F,]
Overlap <- Overlap[,c(2,3)]
colnames(Overlap) <- c("ID","TAD_Zero_oe")

# Add the information to the interaction data frame
BaitList <- merge(BaitList, Overlap, by="ID", all.x=T)

## Baits
Tmp <- BaitList[ duplicated(BaitList$baitID)==F,]
otherEndRanges <- GRanges(seqnames = Tmp$baitChr, IRanges(start = Tmp$baitStart, end = Tmp$baitEnd), strand = rep("+",nrow(Tmp)), mcols = Tmp$baitID)

# Get the overlap and width of overlap
Overlap <- suppressWarnings(findOverlaps(otherEndRanges, TADRanges))
Width <- as.data.frame(ranges(Overlap, ranges(otherEndRanges), ranges(TADRanges)))
Overlap <- as.data.frame(Overlap)

# Convert to data frames
otherEnds <- as.data.frame(otherEndRanges[Overlap[,1]])
TAD <- as.data.frame(TADRanges[Overlap[,2]])

# Set column names
colnames(otherEnds)[6] <- "baitID"
colnames(TAD)[6] <- "TADID"

# Make an ID for merging
otherEnds$Merge <- seq(1, nrow(otherEnds), by=1)
TAD$Merge <- seq(1, nrow(TAD), by=1)
Width$Merge <- seq(1, nrow(Width), by=1)

# Merge the overlap results
Overlap <- merge(otherEnds[,c(6,7)], TAD[,c(6,7)], by="Merge")
Overlap <- merge(Overlap, Width[,c(3,4)], by="Merge")

# Collaps overlap to only the largest overlap
Overlap <- Overlap[ order(Overlap$baitID, -Overlap$width),]
Overlap <- Overlap[ duplicated(Overlap$baitID)==F,]
Overlap <- Overlap[,c(2,3)]
colnames(Overlap) <- c("baitID","TAD_Zero_bait")

# Add the information to the interaction data frame
BaitList <- merge(BaitList, Overlap, by="baitID", all.x=T)

## Filter interactions
BaitList$TAD_Zero_oe <- as.character(BaitList$TAD_Zero_oe)
BaitList$TAD_Zero_bait <- as.character(BaitList$TAD_Zero_bait)
BaitList <- BaitList[ !is.na(BaitList$TAD_Zero_bait),]
BaitList <- BaitList[ duplicated(BaitList$ID)==F,]

# IntraTAD
Within <- BaitList[ !is.na(BaitList$TAD_Zero_oe) & !is.na(BaitList$TAD_Zero_bait),]
Within <- Within[ Within$TAD_Zero_oe == Within$TAD_Zero_bait,]
Result[i,1] <- nrow(Within) / nrow(BaitList)

# TransTAD
Outside <- BaitList[ !(BaitList$ID %in% Within$ID),]
Result[i,2] <- nrow(Outside) / nrow(BaitList)
}

## Finalize the results
# Make a matrix with the results
Matrix <- matrix(ncol=2, nrow=2)
Matrix[1,1] <- A
Matrix[2,1] <- B
Matrix[1,2] <- mean(Result[,1])
Matrix[2,2] <- mean(Result[,2])

# Plot the barplot
barplot(Matrix, col=c("red3","green3"), las=1, ylab="Fraction of interactions", names=c("TADs","Rand TADs"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)
