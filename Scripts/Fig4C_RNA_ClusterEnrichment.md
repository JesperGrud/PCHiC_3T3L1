```R
# Load libraries
library(DESeq2)
library(org.Mm.eg.db)
library(edgeR)
library(circlize)
library(e1071)
library(fdrtools)

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
Pvalues <- as.matrix(as.data.frame(mcols(DGE)[,grep("WaldPvalue_Condition", colnames(mcols(DGE)))]))
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

# Import transcript length information and merge it
Lengths <- read.delim("Data/RNAseq/Annotation/Gene.lengths", header=T)
RNA <- merge(RNA, Lengths, by="RefSeq")

# Normalize RNA-seq counts to txLength
RNA$D0 <- RNA$D0 / (RNA$txLength/1000)
RNA$H4 <- RNA$H4 / (RNA$txLength/1000)
RNA$D1 <- RNA$D1 / (RNA$txLength/1000)
RNA$D2 <- RNA$D2 / (RNA$txLength/1000)
RNA$D7 <- RNA$D7 / (RNA$txLength/1000)

# Calculate maximal expression
RNA$Max <- apply(RNA[,c("D0","H4","D1","D2","D7")],1,FUN="max")

# Calculate log2 fold changes
RNA$logFC_H4 <- RNA$H4 / RNA$D0
RNA$logFC_H4 <- log2(RNA$logFC_H4)
RNA$logFC_D1 <- RNA$D1 / RNA$D0
RNA$logFC_D1 <- log2(RNA$logFC_D1)
RNA$logFC_D2 <- RNA$D2 / RNA$D0
RNA$logFC_D2 <- log2(RNA$logFC_D2)
RNA$logFC_D7 <- RNA$D7 / RNA$D0
RNA$logFC_D7 <- log2(RNA$logFC_D7)

# Define genes to cluster
RNA$Include <- 0
RNA[ (abs(RNA$logFC_H4) > 0 & RNA$Padj_Four <= 0.01) | (abs(RNA$logFC_D2) > 0 & RNA$Padj_Two <= 0.01), "Include"] <- 1

# Do clustering only for genes significant between D0 and either H4 or D2
rownames(RNA) <- RNA$RefSeq
ClusterData <- as.data.frame(RNA[ ,c("D0","H4","D2")])
ClusterData <- log2(ClusterData)
ClusterData <- as.data.frame(t(scale(t(ClusterData))))
ClusterData <- ClusterData[ rownames(ClusterData) %in% RNA[ RNA$Include == 1,"RefSeq"],]

# Set fuzzy-like parameter
Mest <- 2

# Do the clustering
NoClusters <- 8
set.seed(15)
Clusters <- cmeans(ClusterData, NoClusters, m = Mest, method = "cmeans")
Centers <- Clusters$centers

# Calculate membership for each replicate and the mean within each cluster
Rep1 <- as.matrix(RNA[ RNA$Include == 1 ,c("D0_E1","H4_E1","D2_E1")])
Rep1 <- t(scale(t(Rep1)))
dm <- sapply(seq_len(nrow(Rep1)),function(i) apply(Centers, 1, function(v) sqrt(sum((Rep1[i, ]-v)^2))))
Rep1Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Rep1Membership <- as.data.frame(Rep1Membership)
rownames(Rep1Membership) <- RNA[ RNA$Include == 1,c("RefSeq")]

Rep2 <- as.matrix(RNA[ RNA$Include == 1,c("D0_E2","H4_E2","D2_E2")])
Rep2 <- t(scale(t(Rep2)))
dm <- sapply(seq_len(nrow(Rep2)),function(i) apply(Centers, 1, function(v) sqrt(sum((Rep2[i, ]-v)^2))))
Rep2Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Rep2Membership <- as.data.frame(Rep2Membership)
rownames(Rep2Membership) <- RNA[ RNA$Include == 1,c("RefSeq")]

All <- as.matrix(RNA[ RNA$Include == 1,c("D0","H4","D2")])
All <- log2(All)
All <- t(scale(t(All), center = T, scale = F))
dm <- sapply(seq_len(nrow(All)),function(i) apply(Centers, 1, function(v) sqrt(sum((All[i, ]-v)^2))))
Membership <- t(apply(dm, 2, function(x) {
  tmp <- 1/((x/sum(x))^(2/(Mest-1)))
  tmp/sum(tmp)
}))
Membership <- as.data.frame(Membership)
rownames(Membership) <- RNA[ RNA$Include == 1,c("RefSeq")]

# Define clusters based on membership using replicates
Rep1Membership$Max <- apply(Rep1Membership,1,FUN="max")
Rep1Membership$Include <- 0
Rep1Membership[ Rep1Membership$Max >= 0.2,"Include"] <- 1
Rep1Membership$Which <- 0
for (i in 1:nrow(Rep1Membership)) { Rep1Membership[i,"Which"] <- which.max(Rep1Membership[i,c(1:(NoClusters))]) }
Rep1Membership$Which <- Rep1Membership$Which * Rep1Membership$Include
Rep1Membership$RefSeq <- rownames(Rep1Membership)

Rep2Membership$Max <- apply(Rep2Membership,1,FUN="max")
Rep2Membership$Include <- 0
Rep2Membership[ Rep2Membership$Max >= 0.2,"Include"] <- 1
Rep2Membership$Which <- 0
for (i in 1:nrow(Rep2Membership)) { Rep2Membership[i,"Which"] <- which.max(Rep2Membership[i,c(1:(NoClusters))]) }
Rep2Membership$Which <- Rep2Membership$Which * Rep2Membership$Include
Rep2Membership$RefSeq <- rownames(Rep2Membership)

Membership$Max <- apply(Membership,1,FUN="max")
Membership$Include <- 0
Membership[ Membership$Max >= 0.3,"Include"] <- 1
Membership$Which <- 0
for (i in 1:nrow(Membership)) { Membership[i,"Which"] <- which.max(Membership[i,c(1:(NoClusters))]) }
Membership$Which <- Membership$Which * Membership$Include
Membership$RefSeq <- rownames(Membership)

Clusters <- merge(Membership[,c("Max","RefSeq","Which")], Rep1Membership[,c("RefSeq","Which")], by="RefSeq")
Clusters <- merge(Clusters, Rep2Membership[,c("RefSeq","Which")], by="RefSeq")
Clusters <- Clusters[ Clusters$Which.x != 0,]
Clusters <- Clusters[ Clusters$Which.x == Clusters$Which.y,]
Clusters <- Clusters[ Clusters$Which.x == Clusters$Which,]
colnames(Clusters)[3] <- "Cluster"
colnames(Clusters)[2] <- "Membership"
RNA <- merge(RNA[,c(1:29)], Clusters[,c(1,2,3)], by="RefSeq", all.x=T)
RNA[ is.na(RNA$Cluster), "Cluster"] <- 0
RNA[ is.na(RNA$Membership), "Membership"] <- 0

# Filter the results by expression
RNA[ RNA$Include == 0 & RNA$Cluster < (NoClusters), "Cluster"] <- 0

# Remove too small clusters
RNA[ RNA$Cluster == 5, "Cluster"] <- 0

# Rename clusters for intuitive grouping during plots
RNA[ RNA$Cluster == 4, "Cluster"] <- 10
RNA[ RNA$Cluster == 8, "Cluster"] <- 11
RNA[ RNA$Cluster == 3, "Cluster"] <- 12
RNA[ RNA$Cluster == 7, "Cluster"] <- 13
RNA[ RNA$Cluster == 2, "Cluster"] <- 14
RNA[ RNA$Cluster == 1, "Cluster"] <- 15
RNA[ RNA$Cluster == 6, "Cluster"] <- 16
RNA[ RNA$Cluster == 10, "Cluster"] <- 1
RNA[ RNA$Cluster == 11, "Cluster"] <- 2
RNA[ RNA$Cluster == 12, "Cluster"] <- 3
RNA[ RNA$Cluster == 13, "Cluster"] <- 4
RNA[ RNA$Cluster == 14, "Cluster"] <- 5
RNA[ RNA$Cluster == 15, "Cluster"] <- 6
RNA[ RNA$Cluster == 16, "Cluster"] <- 7

# Calculate the hue for each gene
RNA$V <- apply(RNA[,c("D0","H4","D2")],1,FUN="max")

# Calculate the S values
RNA$Minimum <- apply(RNA[,c("D0","H4","D2")],1,FUN="min")
RNA$S <- 1 - (RNA$Minimum / RNA$V)

# Calculate the H values 
RNA$H <- RNA$D0 + RNA$H4 - RNA$D2 - RNA$V
RNA$H <- RNA$H / (RNA$V * RNA$S)
RNA$H <- RNA$H + 2
RNA$H <- RNA$H * sign(RNA$H4 - RNA$D0)
RNA$H <- RNA$H * 60

# Extract mapping to official symbols
Convert <- suppressMessages(select(org.Mm.eg.db, as.character(RNA$RefSeq), columns = c("SYMBOL","REFSEQ"), "REFSEQ"))
colnames(Convert) <- c("RefSeq","Symbol")
RNA <- merge(RNA, Convert, by="RefSeq")

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

# Merge gene information and collaps to keep only unique interaction-gene pairs
FilteredInteractions_Genes <- merge(FilteredInteractions, Genes[,c("Fragment","Symbol")], by.x="baitID", by.y="Fragment")
FilteredInteractions_Genes$Collaps <- paste(FilteredInteractions_Genes$ID, FilteredInteractions_Genes$Symbol, sep="-")
FilteredInteractions_Genes <- FilteredInteractions_Genes[ duplicated(FilteredInteractions_Genes$Collaps)==F,]
FilteredInteractions_Genes <- FilteredInteractions_Genes[,c(1:41)]

# Keep only interactions within 100kb
FilteredInteractions_Genes <- FilteredInteractions_Genes[ abs(FilteredInteractions_Genes$dist) <= 100000,]

## Perform enrichment analysis
Binsize <- 40
WindowSize <- 1
Length <- length(seq(-180,180-WindowSize,by=WindowSize))

# Make data frames to capture the results
Enrichment <- data.frame(matrix(ncol=Length, nrow=6))
PBinom <- data.frame(matrix(ncol=Length, nrow=6))

# Bin the hues (not crossing -180|180)
HueLower <- seq(-180,180-Binsize,by=WindowSize)
HueUpper <- seq(-180+Binsize,180,by=WindowSize)

# Loop through all clusters
for (cluster in 1:7) {
	# Loop through hue ranges within each cluster that dont cross -180|180
	for (hue in 1:Length) {
		# Define the hue range
		Lower <- HueLower[hue]
		Upper <- HueUpper[hue]
		# Get interactions within and outside of hue range
		Within <- FilteredInteractions_Genes[ FilteredInteractions_Genes$H > Lower & FilteredInteractions_Genes$H < Upper,]
		Outside <- FilteredInteractions_Genes[ FilteredInteractions_Genes$H < Lower | FilteredInteractions_Genes$H > Upper,]
		# Get interactions inside and outside gene clusters
		WithinOutside <- Within[ Within$Symbol %in% RNA[ abs(RNA$logFC_H4) <= log2(1.5) & abs(RNA$logFC_D2) <= log2(1.5) & RNA$Padj_Two >= 0.05 & RNA$Padj_Four >= 0.05,"Symbol"],]
		OutsideOutside <- Outside[ Outside$Symbol %in% RNA[ abs(RNA$logFC_H4) <= log2(1.5) & abs(RNA$logFC_D2) <= log2(1.5) & RNA$Padj_Two >= 0.05 & RNA$Padj_Four >= 0.05,"Symbol"],]
		Tmp <- RNA[ RNA$Cluster == cluster,]
		WithinWithin <- Within[ Within$Symbol %in% Tmp[ , "Symbol"],]
		OutsideWithin <- Outside[ Outside$Symbol %in% Tmp[ , "Symbol"],]
		
		# Impose rules about the number of interactions
		# Save enrichment and pvalue
		Enrichment[(cluster),hue] <- log2(((nrow(WithinWithin)+1)/(nrow(WithinWithin)+nrow(WithinOutside)+1))/((nrow(OutsideWithin)+1)/(nrow(OutsideWithin)+nrow(OutsideOutside)+1)))
		PBinom[cluster,hue] <- binom.test(x = nrow(WithinWithin), n = nrow(WithinWithin)+nrow(WithinOutside), p = nrow(OutsideWithin)/(nrow(OutsideWithin)+nrow(OutsideOutside)))$p.value
	}
  # Process hues that cross the boundaries
  for (q in 1:(Binsize-WindowSize)/WindowSize) {
  hue <- length(HueLower) + q
  Lower <- max(HueLower) + (q * WindowSize)
	Upper <- -180 + (q * WindowSize)
	Within <- FilteredInteractions_Genes[ FilteredInteractions_Genes$H > Lower | FilteredInteractions_Genes$H < Upper,]
	Outside <- FilteredInteractions_Genes[ FilteredInteractions_Genes$H < Lower & FilteredInteractions_Genes$H > Upper,]
	# Get interactions inside and outside gene clusters
	WithinOutside <- Within[ Within$Symbol %in% RNA[ abs(RNA$logFC_H4) <= log2(1.5) & abs(RNA$logFC_D2) <= log2(1.5) & RNA$Padj_Two >= 0.05 & RNA$Padj_Four >= 0.05,"Symbol"],]
		OutsideOutside <- Outside[ Outside$Symbol %in% RNA[ abs(RNA$logFC_H4) <= log2(1.5) & abs(RNA$logFC_D2) <= log2(1.5) & RNA$Padj_Two >= 0.05 & RNA$Padj_Four >= 0.05,"Symbol"],]
		Tmp <- RNA[ RNA$Cluster == cluster,]
		WithinWithin <- Within[ Within$Symbol %in% Tmp[ , "Symbol"],]
		
		OutsideWithin <- Outside[ Outside$Symbol %in% Tmp[ , "Symbol"],]
		# Impose rules about the number of interactions
		# Save enrichment and pvalue
		Enrichment[(cluster),hue] <- log2(((nrow(WithinWithin)+1)/(nrow(WithinWithin)+nrow(WithinOutside)+1))/((nrow(OutsideWithin)+1)/(nrow(OutsideWithin)+nrow(OutsideOutside)+1)))
				PBinom[cluster,hue] <- binom.test(x = nrow(WithinWithin), n = nrow(WithinWithin)+nrow(WithinOutside), p = nrow(OutsideWithin)/(nrow(OutsideWithin)+nrow(OutsideOutside)))$p.value
	 }
}

# Perform FDR correction
for (i in 1:7) {
PBinom[i,] <- fdrtool(as.vector(as.matrix(PBinom[i,])), statistic = "pvalue", plot=F, verbose=F)$qval
}

# Plot the enrichments for cluster 1 to 3
par(mfcol=c(1,3))
Limits <- c(1,1,1)
for (cluster in 1:3) {
# Define axis limits
Min=-1 * Limits[cluster]
Max=Limits[cluster]
	# Plot the areas with significant enrichments
	circos.clear()
	circos.par(gap.degree=0, start.degree=-90)
	circos.initialize("A", xlim=c(-179,179))
  # Plot the raw enrichments that dont cross
	circos.trackPlotRegion("A", ylim=c(Min,Max), track.height=0.6, bg.border=NA)	
	for (i in 1:length(HueLower)) {circos.rect(HueLower[i]+19.5, 0, HueUpper[i]-19.5, Enrichment[cluster,i], col = ifelse(PBinom[cluster,i] <= 0.1, ifelse(Enrichment[cluster,i] < 0, "darkred","darkgreen"),"grey"), border = NA) }
	# Add the ones that cross in a loop except the one on both sides
	 # Process hues that cross the boundaries
  for (q in 1:(Binsize-WindowSize)/WindowSize) {
  Lower <- max(HueLower) + (q * WindowSize)+19.5
	Upper <- max(HueLower) + (q * WindowSize)+20.5
  circos.rect(Lower, 0, Upper, Enrichment[cluster,(length(HueLower)+q)], col = ifelse(PBinom[cluster,(length(HueLower)+q)] <= 0.1,ifelse(Enrichment[cluster,(length(HueLower)+q)] < 0, "darkred","darkgreen"),"grey"), border = NA)
	}
  # Add lines
	circos.lines(c(-180,180),c(0,0))
	circos.lines(c(-180,180),c(Max,Max),lty=2)
	circos.lines(c(-180,180),c(Min,Min),lty=2)
	if ( cluster < 7 ){
	  circos.lines(c(max(RNA[ RNA$Cluster == cluster, "H"])+10,max(RNA[ RNA$Cluster == cluster, "H"])+10),c(Min,Max), col="blue", lwd=2)
	  circos.lines(c(min(RNA[ RNA$Cluster == cluster, "H"])-10,min(RNA[ RNA$Cluster == cluster, "H"])-10),c(Min,Max), col="blue", lwd=2)
	} else {
	  circos.lines(c(min(RNA[ RNA$Cluster == cluster & RNA$H > 0, "H"])-10,min(RNA[ RNA$Cluster == cluster & RNA$H > 0, "H"]))-10,c(Min,Max), col="blue", lwd=2)
    circos.lines(c(max(RNA[ RNA$Cluster == cluster & RNA$H < 0, "H"])+10,max(RNA[ RNA$Cluster == cluster & RNA$H < 0, "H"])+10),c(Min,Max), col="blue", lwd=2)
  }
	# Add cluster indicator
	circos.trackPlotRegion("A", ylim=c(0,1),track.height=0.2, bg.border=NA)
	circos.text(0,-1.5,paste("Cluster",cluster,sep=" "))
}
```

[Back to start](../README.md)<br>
[Back to overview of Figure 4](../Links/Figure4.md)
