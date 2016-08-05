```R
## Load packages
library(GenomicRanges)
library(regioneR)

## Process TADs
# Import TADs
Zero_TAD <- read.delim("Data/Interactions/TADs/D0_40kb.domains.bed", header=F, skip=1)
Four_TAD <- read.delim("Data/Interactions/TADs/4h_40kb.domains.bed", header=F, skip=1)
Two_TAD <- read.delim("Data/Interactions/TADs/D2_40kb.domains.bed", header=F, skip=1)
ES_TAD <- read.delim("Data/Interactions/TADs/ESC_40kb.domains.bed", header=F, skip=1)

# Calculate total area of TADs in each condition
Zero_TAD$Size <- Zero_TAD[,3] - Zero_TAD[,2]
Four_TAD$Size <- Four_TAD[,3] - Four_TAD[,2]
Two_TAD$Size <- Two_TAD[,3] - Two_TAD[,2]
ES_TAD$Size <- ES_TAD[,3] - ES_TAD[,2]

# Define unique IDs
Zero_TAD[,4] <- paste("Zero_",seq(1,nrow(Zero_TAD),by=1),sep="")
Four_TAD[,4] <- paste("Four_",seq(1,nrow(Four_TAD),by=1),sep="")
Two_TAD[,4] <- paste("Two_",seq(1,nrow(Two_TAD),by=1),sep="")
ES_TAD[,4] <- paste("ES_",seq(1,nrow(ES_TAD),by=1),sep="")

# Define GRanges objects
ZeroRanges <- GRanges(seqnames = Zero_TAD$V1, IRanges(start = Zero_TAD$V2, end = Zero_TAD$V3), strand = rep("+",nrow(Zero_TAD)), mcols = Zero_TAD$V4)
FourRanges <- GRanges(seqnames = Four_TAD$V1, IRanges(start = Four_TAD$V2, end = Four_TAD$V3), strand = rep("+",nrow(Four_TAD)), mcols = Four_TAD$V4)
TwoRanges <- GRanges(seqnames = Two_TAD$V1, IRanges(start = Two_TAD$V2, end = Two_TAD$V3), strand = rep("+",nrow(Two_TAD)), mcols = Two_TAD$V4)
ESRanges <- GRanges(seqnames = ES_TAD$V1, IRanges(start = ES_TAD$V2, end = ES_TAD$V3), strand = rep("+",nrow(ES_TAD)), mcols = ES_TAD$V4)

# Overlap with D0
FourOverlap <- suppressWarnings(findOverlaps(ZeroRanges, FourRanges))
FourWidth <- as.data.frame(ranges(FourOverlap, ranges(ZeroRanges), ranges(FourRanges)))
TwoOverlap <- suppressWarnings(findOverlaps(ZeroRanges, TwoRanges))
TwoWidth <- as.data.frame(ranges(TwoOverlap, ranges(ZeroRanges), ranges(TwoRanges)))
ESOverlap <- suppressWarnings(findOverlaps(ZeroRanges, ESRanges))
ESWidth <- as.data.frame(ranges(ESOverlap, ranges(ZeroRanges), ranges(ESRanges)))

# Merge width of overlap and TAD IDs
Four <- as.data.frame(cbind(as.data.frame(ZeroRanges[as.matrix(FourOverlap)[,1]])[,6],as.data.frame(FourRanges[as.matrix(FourOverlap)[,2]])[,6],FourWidth$width))
Two <- as.data.frame(cbind(as.data.frame(ZeroRanges[as.matrix(TwoOverlap)[,1]])[,6],as.data.frame(TwoRanges[as.matrix(TwoOverlap)[,2]])[,6],TwoWidth$width))
ES <- as.data.frame(cbind(as.data.frame(ZeroRanges[as.matrix(ESOverlap)[,1]])[,6],as.data.frame(ESRanges[as.matrix(ESOverlap)[,2]])[,6],ESWidth$width))

# Convert width to numeric
Four$V3 <- as.numeric(as.character(Four$V3))
Two$V3 <- as.numeric(as.character(Two$V3))
ES$V3 <- as.numeric(as.character(ES$V3))

# Keep only the largest overlap
Four <- Four[ order(Four$V1, -Four$V3),]
Two <- Two[ order(Two$V1, -Two$V3),]
ES <- ES[ order(ES$V1, -ES$V3),]
Four <- Four[ duplicated(Four$V1)==F,]
Two <- Two[ duplicated(Two$V1)==F,]
ES <- ES[ duplicated(ES$V1)==F,]

# Calculate fractional overlap
Four_Fraction <- sum(Four$V3) / (as.numeric(sum(Zero_TAD$Size)) + as.numeric(sum(Four_TAD$Size)) - as.numeric(sum(Four$V3)))
Two_Fraction <- sum(Two$V3) / (as.numeric(sum(Zero_TAD$Size)) + as.numeric(sum(Two_TAD$Size)) - as.numeric(sum(Two$V3)))
ES_Fraction <- sum(ES$V3) / (as.numeric(sum(Zero_TAD$Size)) + as.numeric(sum(ES_TAD$Size)) - as.numeric(sum(ES$V3)))

# Calculate overlap with randomized control
ZeroRanges <- GRanges(seqnames = Zero_TAD$V1, IRanges(start = Zero_TAD$V2, end = Zero_TAD$V3), strand = rep("+",nrow(Zero_TAD)), mcols = Zero_TAD$V4)
mm9 <- suppressWarnings(fetchExtendedChromInfoFromUCSC("mm9")[c(1:22),c(1,2)])
RandomResult <- data.frame(matrix(ncol=1,nrow=100))
for (i in 1:100) {
	RandomRanges <- randomizeRegions(ZeroRanges, mm9)
	RandomOverlap <- suppressWarnings(findOverlaps(ZeroRanges, RandomRanges))
	RandomWidth <- as.data.frame(ranges(RandomOverlap, ranges(ZeroRanges), ranges(RandomRanges)))
	Random <- as.data.frame(cbind(as.data.frame(ZeroRanges[as.matrix(RandomOverlap)[,1]])[,6],RandomWidth$width))
	Random$V2 <- as.numeric(as.character(Random$V2))
	Random <- Random[ order(Random$V1, -Random$V2),]
	Random <- Random[ duplicated(Random$V1)==F,]
	RandomResult[i,1] <- sum(Random$V2) / (as.numeric(sum(Zero_TAD$Size)) + as.numeric(sum(Zero_TAD$Size)) - as.numeric(sum(Random$V2)))
}

## Plot the overlap
barplot(c(Four_Fraction, Two_Fraction, ES_Fraction, mean(RandomResult[,1])), las=1, ylim=c(0,1), col=c("darkred","green","chocolate4","black"), names=c("4h","D2","ESC","Random"), ylab="Fractional overlap with D0")
```

[Back to start](../README.md)<br>
[Back to overview of Figure 1](../Links/Figure1.md)
