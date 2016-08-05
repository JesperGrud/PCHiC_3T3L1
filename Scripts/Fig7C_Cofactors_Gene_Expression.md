```R
## Load the necessary libraries
library(GenomicRanges)
library(edgeR)
library(DESeq2)
library(org.Mm.eg.db)

## Process cofactor ChIP-seq data
# Import read counts for the 6 cofactors
Cofactors <- read.delim("Data/ChIPseq/Counts/Cofactors.count")

# Take means of replicates
Cofactors$P300_4h <- rowMeans(Cofactors[,c(8,9)])
Cofactors$P300_D0 <- rowMeans(Cofactors[,c(10,11)])
Cofactors$MED1_4h <- rowMeans(Cofactors[,c(12,13)])
Cofactors$MED1_D0 <- rowMeans(Cofactors[,c(14,15)])
Cofactors$SMC1_4h <- rowMeans(Cofactors[,c(16,17)])
Cofactors$SMC1_D0 <- rowMeans(Cofactors[,c(18,19)])
Cofactors$HDAC2_4h <- rowMeans(Cofactors[,c(20,21)])
Cofactors$HDAC2_D0 <- rowMeans(Cofactors[,c(22,23)])
Cofactors$HDAC3_4h <- rowMeans(Cofactors[,c(24,25)])
Cofactors$HDAC3_D0 <- rowMeans(Cofactors[,c(26,27)])
Cofactors$NCoR_4h <- rowMeans(Cofactors[,c(28,29)])
Cofactors$NCoR_D0 <- rowMeans(Cofactors[,c(30,31)])

# Filter away weak sites (< 20 tags / 10 mill in all samples)
Cofactors$Max <- apply(Cofactors[,c(8:31)],1,FUN="max")
Cofactors <- Cofactors[ Cofactors$Max >= 20,]

# Filter away clonal sites (>= 100 tags / 10 mill in all samples)
Cofactors$Min <- apply(Cofactors[,c(8:31)],1,FUN="min")
Cofactors <- Cofactors[ Cofactors$Min <= 100,]

# Filter columns and set names
Cofactors <- Cofactors[,c(1:4,32:43)]
colnames(Cofactors)[1] <- "PeakID"

# Calculate log2 fold changes
Cofactors$logFC_MED1 <- log2((Cofactors$MED1_4h+0.25) / (Cofactors$MED1_D0+0.25))
Cofactors$logFC_P300 <- log2((Cofactors$P300_4h+0.25) / (Cofactors$P300_D0+0.25))
Cofactors$logFC_SMC1 <- log2((Cofactors$SMC1_4h+0.25) / (Cofactors$SMC1_D0+0.25))
Cofactors$logFC_HDAC2 <- log2((Cofactors$HDAC2_4h+0.25) / (Cofactors$HDAC2_D0+0.25))
Cofactors$logFC_HDAC3 <- log2((Cofactors$HDAC3_4h+0.25) / (Cofactors$HDAC3_D0+0.25))
Cofactors$logFC_NCoR <- log2((Cofactors$NCoR_4h+0.25) / (Cofactors$NCoR_D0+0.25))

# Calculate average log2 fold induction of coactivators and corepressors
Cofactors$Average_Coactivator <- rowMeans(Cofactors[,c("logFC_MED1","logFC_SMC1","logFC_P300")])
Cofactors$Average_Corepressor <- rowMeans(Cofactors[,c("logFC_NCoR","logFC_HDAC3","logFC_HDAC2")])

# Z-transform the log2 fold changes
for (i in 17:22) {
  Mean <- mean(Cofactors[,i])
  SD <- sd(Cofactors[,i])
  Cofactors[,i] <- (Cofactors[,i]-Mean)/SD
}

# Clean up the workspace
rm(list=c("Mean","SD","i"))

# Calculate composite scores
Cofactors$Coactivators <- apply(Cofactors[,c(17,18,19)],1,FUN="sum")
Cofactors$Corepressors <- apply(Cofactors[,c(20,21,22)],1,FUN="sum")

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

# Clean up the workspace
rm(list=setdiff(ls(),c("FilteredInteractions_Genes","Cofactors")))

## Overlap the otherEnds of interactions with cofactor peaks
# Create GRanges object for enhancers and otherEnds
EnhancerRanges <- GRanges(seqnames = Cofactors$Chr, IRanges(start = Cofactors$Start, end = Cofactors$End), strand = rep("+",nrow(Cofactors)), mcols = Cofactors$PeakID)
otherEnds <- FilteredInteractions_Genes[ duplicated(FilteredInteractions_Genes$oeID)==F,]
otherEndRanges <- GRanges(seqnames = otherEnds$oeChr, IRanges(start = otherEnds$oeStart, end = otherEnds$oeEnd), strand = rep("+",nrow(otherEnds)), mcols = otherEnds$oeID)

# Get the overlap and width of overlap
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

# Combine the interaction data and cofactor data
Overlap <- merge(Overlap, Cofactors[,c("Corepressors","Coactivators","PeakID")], by="PeakID")
Overlap <- merge(Overlap, FilteredInteractions_Genes[,c("oeID","dist","Symbol","logFC_Four_Zero","ID")], by="oeID")

## Process RNA-seq data
# Import RNA-seq data
RNA <- read.delim("Data/RNAseq/RNA.txt", header=T)

# Perform differential expression analysis using DESeq2
MetaData <- data.frame(Condition = c("D0","H4","D1","D2","D7","D0","H4","D1","D2","D7"), Replicate = c(rep("a",5), rep("b",5)))
rownames(MetaData) <- colnames(RNA)[c(6:15)]
rownames(RNA) <- RNA$RefSeq
DGE <- DESeqDataSetFromMatrix(as.matrix(RNA[,c(6:15)]), colData = MetaData, design = ~ Condition)
DGE <- estimateSizeFactors(DGE)
DGE <- estimateDispersions(DGE)
DGE <- nbinomWaldTest(DGE, betaPrior=F)

# Extract p-values and do correction for multiple testing
Pvalues <- as.matrix(mcols(DGE)[,grep("WaldPvalue_Condition", colnames(mcols(DGE)))])
Padj <- as.data.frame(apply(as.matrix(Pvalues),2,function(X) {p.adjust(X, method='BH')}))
colnames(Padj) <- c("Padj_One","Padj_Two","Padj_Seven","Padj_Four")
Padj$RefSeq <- names(DGE)

# Get variance stabilized data 
VSD <- as.data.frame(assay(varianceStabilizingTransformation(DGE, blind=F)))
VSD$RefSeq <- rownames(VSD)

# Merge all of the data and remove NAs (All zero counts)
RNA <- merge(RNA[,c(1:5)], VSD, by="RefSeq", sort=F)
RNA <- merge(RNA, Padj, by="RefSeq", sort=F)
RNA <- na.omit(RNA)

# Calculate mean expression across replicates
RNA$D0 <- rowMeans(2^RNA[,c("D0_E1","D0_E2")])
RNA$H4 <- rowMeans(2^RNA[,c("H4_E1","H4_E2")])
RNA$D1 <- rowMeans(2^RNA[,c("D1_E1","D1_E2")])
RNA$D2 <- rowMeans(2^RNA[,c("D2_E1","D2_E2")])
RNA$D7 <- rowMeans(2^RNA[,c("D7_E1","D7_E2")])

# Do clustering only for genes significant between D0 and either H4 or D2
rownames(RNA) <- RNA$RefSeq
ClusterData <- as.matrix(RNA[ RNA$Padj_Two < 0.01 | RNA$Padj_Four < 0.01,c(6,7,9,11,12,14)])
set.seed(13)
Clusters <- kmeans(t(scale(t(ClusterData), center=T, scale=F)), 7, iter.max=500)
Clusters <- data.frame(Cluster = Clusters$cluster)
Clusters$RefSeq <- rownames(Clusters)
RNA <- merge(RNA, Clusters, by="RefSeq", all.x=T)
RNA[ is.na(RNA$Cluster), "Cluster"] <- 0

# Import transcript length information and merge it
Lengths <- read.delim("Data/RNAseq/Annotation/Gene.lengths", header=T)
RNA <- merge(RNA, Lengths, by="RefSeq")

# Extract mapping to official symbols
Convert <- suppressMessages(select(org.Mm.eg.db, as.character(RNA$RefSeq), columns = c("SYMBOL","REFSEQ"), "REFSEQ"))
colnames(Convert) <- c("RefSeq","Symbol")
RNA <- merge(RNA, Convert, by="RefSeq")

# Calculate log2 fold change for genes after 4 hours
RNA$logFC_Four <- log2(RNA$H4 / RNA$D0)

## Loop through bins of coactivator and corepressor composite scores and calculate the average log2 fold change of interactions
# Make a data frame to capture the results
Heatmap_Genes <- data.frame(matrix(ncol=length(seq(-12.4,12.4, by=0.4))-1,nrow=length(seq(-12.4,12.4, by=0.4))-1))

# Loop through the data and count the number of sites
for (i in 1:(length(seq(-12.4,12.4, by=0.4))-1)) {
  Coact_Low <- seq(-12.4,12.4, by=0.4)[i]
  Coact_High <- seq(-12.4,12.4, by=0.4)[(i+1)]
  for (q in 1:(length(seq(-12.4,12.4, by=0.4))-1) ) {
    Corep_Low <- seq(-12.4,12.4, by=0.4)[q]
    Corep_High <- seq(-12.4,12.4, by=0.4)[(q+1)]
    Tmp <- Overlap[ Overlap$Coactivators >= Coact_Low & Overlap$Coactivators < Coact_High & Overlap$Corepressors >= Corep_Low & Overlap$Corepressors < Corep_High,]
    if (nrow(Tmp) == 0) { Heatmap_Genes[i,q] <- NA } else { 
		Tmp <- merge(Tmp, RNA[,c("Symbol","logFC_Four")], by="Symbol")
		Heatmap_Genes[i,q] <- mean(Tmp$logFC_Four) 
		}
  }
}

# Softcap the log2 fold changes at 2 
Heatmap_Genes[ Heatmap_Genes > 2 & !is.na(Heatmap_Genes)] <- 2
Heatmap_Genes[ Heatmap_Genes < -2 & !is.na(Heatmap_Genes)] <- -2

# Setup for plotting the heatmap
Labels <- signif(seq(-12.4,12.4, by=0.4),3)
Breaks <- seq(-2,2,length.out = 111)
Colors <- colorRampPalette(c("darkred","red","white","green","darkgreen"))(110)

# Plot the heatmap
par(fig=c(0,0.8,0,1))
plot(Labels, Labels, type="n", axes=T, xlab="Coactivators", ylab="Corepressors",,xaxs="i", yaxs="i", las=1)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
image(Labels, Labels, as.matrix(Heatmap_Genes), breaks=Breaks,col=Colors, add=T)

# Plot the legend
par(fig=c(0.8,1,0,1), new=t)
plot(1,1,t="n",ylim=c(min(Breaks),max(Breaks)), xlim=c(0,1), xaxt="n", yaxt="s", xlab="", ylab="Mean log2 fold change",xaxs="i", yaxs="i", las = 1, frame.plot=F) 
for (i in 1:(length(Breaks)-1)) {
	rect(0,Breaks[[i]],1,Breaks[[(i+1)]], col=Colors[[i]], border=NA)
}
box()
```

[Back to start](../README.md)<br>
[Back to overview of Figure 7](../Links/Figure7.md)