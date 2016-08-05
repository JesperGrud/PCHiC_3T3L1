```R
## Load the necessary packages
library(edgeR)
library(GenomicRanges)
library(DESeq2)
library(circlize)

## Process enhancer data
# Import the data
Counts <- read.delim("Data/ChIPseq/Counts/MED1.count")

# Keep only distal sites (>= 2kb away from closest TSS)
Counts <- Counts[ abs(Counts$DistanceToNearestGene) >= 2000,]

# Setup to run DESeq2
rownames(Counts) <- Counts$PeakID
Design <- data.frame(Condition = c("D0","h4","D2","D4","D7","D0","h4","D2","D4","D7"), Replicate = c(rep("a",5),rep("b",5)))
rownames(Design) <- colnames(Counts)[c(6:15)]

# Initialize DESeq2
DDS <- DESeqDataSetFromMatrix(countData = as.matrix(Counts[,c(6:15)]), colData = Design, design = ~Condition)

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
Result=cbind(Counts[,c(1:4)],assay(Normalized),D0,H4,D2,D4,D7)

# Identify differentially occupied enhancers (FDR-adjusted p-values <= 0.05)
Result$Differential <- 0
Result[ apply(Result[,c(15:34)],1,FUN="min") <= 0.05,"Differential"] <- 1

# Calculate average occupancy in replicates
Result$Count_D0 <- rowMeans(Result[,c("MED1_D0_exp1","MED1_D0_exp2")])
Result$Count_H4 <- rowMeans(Result[,c("MED1_4h_exp1","MED1_4h_exp2")])
Result$Count_D2 <- rowMeans(Result[,c("MED1_D2_exp1","MED1_D2_exp2")])
Result$Count_D4 <- rowMeans(Result[,c("MED1_D4_exp1","MED1_D4_exp2")])
Result$Count_D7 <- rowMeans(Result[,c("MED1_D7_exp1","MED1_D7_exp2")])

# Center the averaged counts
Result[,c(36:40)] <- t(scale(t(Result[,c(36:40)]), center=T, scale=F))

# Partition data into 5 cluster using k-means for the differential super-enhancers (on both replicates)
rownames(Result) <- Result$PeakID
set.seed(13)
Clusters <- kmeans(t(scale(t(as.matrix(Result[ Result$Differential == 1,c(5:14)])), center=T, scale=F)), 7, iter.max=500)
Clustering <- data.frame(Cluster = Clusters$cluster)
Clustering$PeakID <- rownames(Clustering)

# Rename clusters for intuitive cluster profiles
Clustering[ Clustering$Cluster == 1, "Cluster"] <- 8
Clustering[ Clustering$Cluster == 2, "Cluster"] <- 9
Clustering[ Clustering$Cluster == 3, "Cluster"] <- 10
Clustering[ Clustering$Cluster == 4, "Cluster"] <- 11
Clustering[ Clustering$Cluster == 6, "Cluster"] <- 12
Clustering[ Clustering$Cluster == 7, "Cluster"] <- 13
Clustering[ Clustering$Cluster == 8, "Cluster"] <- 6
Clustering[ Clustering$Cluster == 9, "Cluster"] <- 3
Clustering[ Clustering$Cluster == 10, "Cluster"] <- 7
Clustering[ Clustering$Cluster == 11, "Cluster"] <- 1
Clustering[ Clustering$Cluster == 12, "Cluster"] <- 4
Clustering[ Clustering$Cluster == 13, "Cluster"] <- 2

## Process interactions
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

## Overlap the data  and the interactions
# Get peaks from each group
C1 <- Counts[ Counts$PeakID %in% Clustering[ Clustering$Cluster == 1, "PeakID"],c(1:4)]
C2 <- Counts[ Counts$PeakID %in% Clustering[ Clustering$Cluster == 2, "PeakID"],c(1:4)]
C3 <- Counts[ Counts$PeakID %in% Clustering[ Clustering$Cluster == 3, "PeakID"],c(1:4)]
C4 <- Counts[ Counts$PeakID %in% Clustering[ Clustering$Cluster == 4, "PeakID"],c(1:4)]
C5 <- Counts[ Counts$PeakID %in% Clustering[ Clustering$Cluster == 5, "PeakID"],c(1:4)]
C6 <- Counts[ Counts$PeakID %in% Clustering[ Clustering$Cluster == 6, "PeakID"],c(1:4)]
C7 <- Counts[ Counts$PeakID %in% Clustering[ Clustering$Cluster == 7, "PeakID"],c(1:4)]

## Overlap SE constituents and otherEnds
# Make a GRanges object from otherEnds
otherEnds <- FilteredInteractions[ duplicated(FilteredInteractions$oeID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$oeChr, IRanges(start = otherEnds$oeStart, end = otherEnds$oeEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$oeID)
	
# Make a list with all constituents
Constituents <- list(C1, C2, C3, C4, C5, C6, C7)

# Define bin size for calculating enrichments
Binsize <- 20

# Make a list for capturing the results
ResultList <- list()

# Loop through all sets and plot
for (consti in 1:7) {
	# Grab the data from the list
	tmp <- Constituents[[consti]]
	# Create GRanges object for enhancers and otherEnds
	ConstRanges <- GRanges(seqnames = tmp$Chr, IRanges(start = tmp$Start, end = tmp$End), strand = rep("+",nrow(tmp)), mcols = tmp$PeakID)
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
	Result <- data.frame(matrix(ncol=5, nrow= length(seq(-180,160,by=Binsize))))
	# Set column names and define upper and lower bounds
	colnames(Result) <- c("LowerBound","UpperBound","Enrichment","P_Less","P_Greater")
	Result[,"LowerBound"] <- seq(-180,180-Binsize,by=Binsize)
	Result[,"UpperBound"] <- seq(-180+Binsize,180,by=Binsize)
	# Calculate enrichment and P-values in a loop
	for (i in 1:nrow(Result)) {
		# Devide interaction into a group within the hue range and a group outside
		Within <- FilteredInteractions[ FilteredInteractions$H >= as.numeric(Result[i,"LowerBound"]) & FilteredInteractions$H <= as.numeric(Result[i,"UpperBound"]),]
		Outside <- FilteredInteractions[ FilteredInteractions$H < Result[i,"LowerBound"] | FilteredInteractions$H > Result[i,"UpperBound"],]
		# Get the interactions within the hue range overlapping and not overlapping constituents
		WithinOverlap <- Within[ Within$oeID %in% ConstOverlap$oeID,]
		WithinNoOverlap <- Within[ !(Within$oeID %in% ConstOverlap$oeID),]
		# Get the interactions outside the hue range overlapping and not overlapping constituents
		OutsideOverlap <- Outside[ Outside$oeID %in% ConstOverlap$oeID,]
		OutsideNoOverlap <- Outside[ !(Outside$oeID %in% ConstOverlap$oeID),]
		# Calculate enrichment and p-values
		Result[i,"Enrichment"] <- log2((nrow(WithinOverlap) / ( nrow(WithinOverlap) + nrow(WithinNoOverlap) ) ) / (nrow(OutsideOverlap) / ( nrow(OutsideOverlap) + nrow(OutsideNoOverlap) ) ))
		Result[i,"P_Less"] <- phyper(nrow(WithinOverlap),nrow(OutsideOverlap), nrow(OutsideNoOverlap), nrow(WithinOverlap)+nrow(WithinNoOverlap), lower.tail=TRUE)
		Result[i,"P_Greater"] <- phyper(nrow(WithinOverlap),nrow(OutsideOverlap), nrow(OutsideNoOverlap), nrow(WithinOverlap)+nrow(WithinNoOverlap), lower.tail=FALSE)
		}
	# Place the result into a new list
	ResultList[[consti]] <- Result
	}
	
## Plot the enrichments in a loop
par(mfcol=c(2,4))
Limits <- c(1.25,0.75,1,0.5,0.75,0.75,1.25)
for (q in 1:7) {
	# Grab the data
	tmp <- ResultList[[q]]
	# Set min and max
	max <- Limits[q]
	min <- -1 * max
	# Initialize the plot
	circos.clear()
	circos.par(gap.degree=0, start.degree=-90)
	circos.initialize("A", xlim=c(-179,179))
	circos.trackPlotRegion("A", ylim=c(min,max), track.height=0.5, bg.border=NA)
	# Plot it
	for (i in 1:18) {circos.rect(tmp[i,1], 0, tmp[i,2], tmp[i,3], col=ifelse(tmp[i,4]<=0.05, "darkred", ifelse(tmp[i,5]<=0.05,"darkgreen","grey"))) }
	circos.lines(c(-180,180),c(0,0))
	circos.lines(c(-180,180),c(max,max),lty=2)
	circos.lines(c(-180,180),c(min,min),lty=2)
	}
```

[Back to start](../README.md)<br>
[Back to overview of Figure S3](../Links/FigureS3.md)