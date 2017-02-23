```R
## Load packages
library(DESeq2)
library(edgeR)
library(org.Mm.eg.db)
library(GenomicRanges)
library(rtracklayer)

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

# Add gene information 
FilteredInteractions_Genes <- merge(FilteredInteractions, Genes[,c("Fragment","Symbol")], by.x="baitID", by.y="Fragment")
FilteredInteractions_Genes$Collaps <- paste(FilteredInteractions_Genes$ID, FilteredInteractions_Genes$Symbol, sep="-")
FilteredInteractions_Genes <- FilteredInteractions_Genes[ duplicated(FilteredInteractions_Genes$Collaps)==F,]
FilteredInteractions_Genes <- FilteredInteractions_Genes[,c(1:41)]

## Process eQTL data
# Import the data
eQTL <- read.delim("Data/Other/eQTL.txt", header=T)
eQTL$UniqueID <- seq(1, nrow(eQTL), by=1)

## Lift over SNP positions from mm10 to mm9
# Liftover
ch = import.chain("Data/Other/Annotation/mm10ToMm9.over.chain")
eQTLRange <- GRanges(seqnames = paste("chr",eQTL$snp_chr,sep=""), IRanges(start = eQTL$snp_bp_mm10, end = eQTL$snp_bp_mm10), strand = rep("+",nrow(eQTL)), mcols=eQTL$UniqueID)
LO <- liftOver(eQTLRange, ch)
LO <- as.data.frame(unlist(LO))
colnames(LO)[c(2,6)] <- c("snp_bp_mm9","UniqueID")
LO <- LO[ order(LO$UniqueID),]
# Combine the results
Unchained <- eQTL[ !(eQTL$UniqueID %in% LO$UniqueID),]
Unchained$snp_bp_mm9 <- Unchained$snp_bp_mm10
Chained <- eQTL[ (eQTL$UniqueID %in% LO$UniqueID),]
Chained <- Chained[ order(Chained$UniqueID),]
Chained <- cbind(Chained, LO$snp_bp_mm9)
colnames(Chained)[20] <- "snp_bp_mm9"
eQTL <- rbind(Chained, Unchained)

## Lift over TSS positions from mm10 to mm9
# Liftover
eQTLRange <- GRanges(seqnames = paste("chr",eQTL$snp_chr,sep=""), IRanges(start = eQTL$gene_start_bp, end = eQTL$gene_start_bp), strand = rep("+",nrow(eQTL)), mcols=eQTL$UniqueID)
LO <- liftOver(eQTLRange, ch)
LO <- as.data.frame(unlist(LO))
colnames(LO)[c(2,6)] <- c("gene_start_bp_mm9","UniqueID")
LO <- LO[ order(LO$UniqueID),]
# Combine the results
Unchained <- eQTL[ !(eQTL$UniqueID %in% LO$UniqueID),]
Unchained$gene_start_bp_mm9 <- Unchained$gene_start_bp
Chained <- eQTL[ (eQTL$UniqueID %in% LO$UniqueID),]
Chained <- Chained[ order(Chained$UniqueID),]
Chained <- cbind(Chained, LO$gene_start_bp_mm9)
colnames(Chained)[21] <- "gene_start_bp_mm9"
eQTL <- rbind(Chained, Unchained)

# Remove eQTLs spanning more than 1Mb
eQTL$snp_to_gene_start <- abs(eQTL$gene_start_bp_mm9 - eQTL$snp_bp_mm9)
eQTL <- eQTL[ eQTL$snp_to_gene_start <= 1000000, ]

# Keep only the lead eQTL for each gene
eQTL <- eQTL[ order(eQTL$gene_symbol, eQTL$pvalue),]
eQTL <- eQTL[ duplicated(eQTL$gene_symbol)==F,]

# Keep only interactions mapping to an RefSeq gene
FilteredInteractions_Genes <- FilteredInteractions_Genes[ FilteredInteractions_Genes$Symbol %in% RNA$Symbol,]
FilteredInteractions_Genes$UniqueID <- seq(1, nrow(FilteredInteractions_Genes), by = 1)

# Overlap PIR and eQTLs
eQTLRange <- GRanges(seqnames = paste("chr",eQTL$snp_chr,sep=""), IRanges(start = eQTL$snp_bp_mm9, end = eQTL$snp_bp_mm9), strand = rep("+",nrow(eQTL)), mcols=eQTL$UniqueID)
Tmp <- FilteredInteractions[ duplicated(FilteredInteractions$oeID)==F,]
otherEndRanges <- GRanges(seqnames = Tmp$oeChr, IRanges(start = Tmp$oeStart, end = Tmp$oeEnd), strand = rep("+",nrow(Tmp)), mcols = Tmp$oeID)
Overlap <- as.data.frame(suppressWarnings(findOverlaps(otherEndRanges, eQTLRange)))
colnames(Overlap) <- c("oeRow","eQTLRow")
eQTL$eQTLRow <- seq(1,nrow(eQTL), by=1)
Tmp$oeRow <- seq(1,nrow(Tmp), by=1)
Overlap <- merge(Overlap, Tmp[,c("oeRow","oeID")])
Overlap <- merge(Overlap, FilteredInteractions_Genes[,c("ID","oeID","Symbol","dist","UniqueID")], by="oeID")
Overlap <- merge(Overlap, eQTL[,c("eQTLRow","rsID","pvalue","gene_symbol","snp_to_gene_start")], by="eQTLRow")
Overlap$Symbol <- as.character(Overlap$Symbol)
Overlap$gene_symbol <- as.character(Overlap$gene_symbol)

# Filter SNPs: One per gene and one per otherEnd
SNPs <- Overlap[ duplicated(Overlap$rsID)==F,]
OnTarget <- Overlap[ Overlap$Symbol == Overlap$gene_symbol,]
OnTarget <- OnTarget[ duplicated(OnTarget$rsID)==F,]
All <- nrow(OnTarget)/nrow(SNPs)

## Make a randomized control
# Move otherEnds to another random gene
FilteredInteractions_Genes$Size <- FilteredInteractions_Genes$oeEnd - FilteredInteractions_Genes$oeStart
SummaryData <- data.frame(matrix(ncol=1,nrow=1000))

for (i in 1:1000) {
	# Make a list of all genes and randomize it
	GeneList <- data.frame(Symbol = FilteredInteractions_Genes[ duplicated(FilteredInteractions_Genes$Symbol)==F,c("Symbol")])
	GeneList$Random <- sample(RNA$Symbol, size = nrow(GeneList), replace = T)
	
	# Merge on information about the position of the random gene, and the interactions from the source gene
	GeneList <- merge(GeneList, RNA[,c("Symbol","strand","start","end","chr")], by.x="Random", by.y="Symbol")
	GeneList <- merge(GeneList, FilteredInteractions_Genes[,c("dist","Symbol","ID","Size")], by="Symbol")

	# Move the interactions from the source gene to same relative positions of the random gene
	Plus <- GeneList[ GeneList$strand == "+",]
	Plus$Mid <- Plus$start - Plus$dist
	Plus$oeStart <- Plus$Mid - (Plus$Size/2)
	Plus$oeEnd <- Plus$Mid + (Plus$Size/2)
	Minus <- GeneList[ GeneList$strand == "-",]
	Minus$Mid <- Minus$start - Minus$dist
	Minus$oeStart <- Minus$Mid - (Minus$Size/2)
	Minus$oeEnd <- Minus$Mid + (Minus$Size/2)
	GeneList <- rbind(Plus, Minus)
	colnames(GeneList)[6] <- "oeChr"
	GeneList$oeID <- seq(1, nrow(GeneList), by=1)

	# Overlap the new positions of the PIR with the SNPs
	otherEndRanges <- GRanges(seqnames = GeneList$oeChr, IRanges(start = GeneList$oeStart, end = GeneList$oeEnd), strand = rep("+",nrow(GeneList)), mcols = GeneList$oeID)
	Overlap <- as.data.frame(suppressWarnings(findOverlaps(otherEndRanges, eQTLRange)))
	colnames(Overlap) <- c("oeRow","eQTLRow")
	eQTL$eQTLRow <- seq(1,nrow(eQTL), by=1)
	GeneList$oeRow <- seq(1,nrow(GeneList), by=1)
	Overlap <- merge(Overlap, GeneList[,c("oeRow","Random","dist","ID")])
	Overlap <- merge(Overlap, eQTL[,c("eQTLRow","rsID","pvalue","gene_symbol")])
	Overlap$Random <- as.character(Overlap$Random)
	Overlap$gene_symbol <- as.character(Overlap$gene_symbol)
	SNPs <- Overlap[ duplicated(Overlap$rsID)==F,]
	
	# Find the SNPs that are on target with the random gene
	OnTarget <- Overlap[ Overlap$Random == Overlap$gene_symbol,]
	OnTarget <- OnTarget[ duplicated(OnTarget$rsID)==F,]
	
	# Save the result
	SummaryData[i,1] <- nrow(OnTarget)/nrow(SNPs)
}

# Plot the result
barCenters <- barplot(c(All,mean(SummaryData[,1])), las=1, ylim=c(0,0.12), names = c("PIR","Randomized"), col=c("green3","grey"), ylab="Fraction of eQTLs on  target")
arrows(barCenters[2], mean(SummaryData[,1]) - sd(SummaryData[,1]) * 2, barCenters[2],mean(SummaryData[,1]) + sd(SummaryData[,1]) * 2, lwd = 1.5, angle = 90, code = 3, length = 0.05)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 2](../Links/Figure2.md)
