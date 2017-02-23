```R
## Load packages
library(DESeq2)
library(edgeR)
library(org.Mm.eg.db)

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

# Define groups of dynamic genes and interactions
RepressedGenes <- RNA[ RNA$Padj_Four <= 0.01 & RNA$logFC_Four <= -log2(2),]
InducedGenes <- RNA[ RNA$Padj_Four <= 0.01 & RNA$logFC_Four >= log2(2),]
RepressedInt <- FilteredInteractions_Genes[ FilteredInteractions_Genes$logFC_Four_Zero <= -log2(1.5),]
InducedInt <- FilteredInteractions_Genes[ FilteredInteractions_Genes$logFC_Four_Zero >= log2(1.5),]

# Get the overlap between the different groups and all genes / interactions
Ind_Ind <- InducedInt[ InducedInt$Symbol %in% InducedGenes$Symbol,]
Ind_Rep <- InducedInt[ InducedInt$Symbol %in% RepressedGenes$Symbol,]
Ind_All <- InducedInt[ InducedInt$Symbol %in% RNA$Symbol,]
Rep_Ind <- RepressedInt[ RepressedInt$Symbol %in% InducedGenes$Symbol,]
Rep_Rep <- RepressedInt[ RepressedInt$Symbol %in% RepressedGenes$Symbol,]
Rep_All <- RepressedInt[ RepressedInt$Symbol %in% RNA$Symbol,]
All_Ind <- FilteredInteractions_Genes[ FilteredInteractions_Genes$Symbol %in% InducedGenes$Symbol,]
All_Rep <- FilteredInteractions_Genes[ FilteredInteractions_Genes$Symbol %in% RepressedGenes$Symbol,]
All_All <- FilteredInteractions_Genes[ FilteredInteractions_Genes$Symbol %in% RNA$Symbol,]
  
# Count the number of interactions within different distances
Matrix <- matrix(ncol=22, nrow=10)
Matrix[,1] <- c(50000,100000,150000,200000,250000,300000,350000,400000,450000,500000)
for(i in 1:10) { 
  Matrix[i,2] <- nrow(Ind_Ind[ abs(Ind_Ind$dist) <= Matrix[i,1],]) 
  Matrix[i,3] <- nrow(Ind_Rep[ abs(Ind_Rep$dist) <= Matrix[i,1],]) 
  Matrix[i,4] <- nrow(Ind_All[ abs(Ind_All$dist) <= Matrix[i,1],]) 
  Matrix[i,5] <- nrow(Rep_Ind[ abs(Rep_Ind$dist) <= Matrix[i,1],]) 
  Matrix[i,6] <- nrow(Rep_Rep[ abs(Rep_Rep$dist) <= Matrix[i,1],]) 
  Matrix[i,7] <- nrow(Rep_All[ abs(Rep_All$dist) <= Matrix[i,1],]) 
  Matrix[i,8] <- nrow(All_Ind[ abs(All_Ind$dist) <= Matrix[i,1],]) 
  Matrix[i,9] <- nrow(All_Rep[ abs(All_Rep$dist) <= Matrix[i,1],]) 
  Matrix[i,10] <- nrow(All_All[ abs(All_All$dist) <= Matrix[i,1],]) 
}

# Calculate the P-value (both gain and loss)
for(i in 1:10) {
	Matrix[i,11] <- prop.test(c(Matrix[i,2],Matrix[i,8]),c(Matrix[i,4],Matrix[i,10]), alternative = "greater")$p.value
	Matrix[i,12] <- prop.test(c(Matrix[i,2],Matrix[i,8]),c(Matrix[i,4],Matrix[i,10]), alternative = "less")$p.value
	Matrix[i,13] <- prop.test(c(Matrix[i,3],Matrix[i,9]),c(Matrix[i,4],Matrix[i,10]), alternative = "greater")$p.value
	Matrix[i,14] <- prop.test(c(Matrix[i,3],Matrix[i,9]),c(Matrix[i,4],Matrix[i,10]), alternative = "less")$p.value
	Matrix[i,15] <- prop.test(c(Matrix[i,5],Matrix[i,8]),c(Matrix[i,7],Matrix[i,10]), alternative = "greater")$p.value
	Matrix[i,16] <- prop.test(c(Matrix[i,5],Matrix[i,8]),c(Matrix[i,7],Matrix[i,10]), alternative = "less")$p.value
	Matrix[i,17] <- prop.test(c(Matrix[i,6],Matrix[i,9]),c(Matrix[i,7],Matrix[i,10]), alternative = "greater")$p.value
	Matrix[i,18] <- prop.test(c(Matrix[i,6],Matrix[i,9]),c(Matrix[i,7],Matrix[i,10]), alternative = "less")$p.value
	}
		
# Normalize to the total number of interactions within the distance threshold
Matrix[,2] <- Matrix[,2] / Matrix[,4]
Matrix[,3] <- Matrix[,3] / Matrix[,4]
Matrix[,5] <- Matrix[,5] / Matrix[,7]
Matrix[,6] <- Matrix[,6] / Matrix[,7]
Matrix[,8] <- Matrix[,8] / Matrix[,10]
Matrix[,9] <- Matrix[,9] / Matrix[,10]

# Calculate log2 enrichment compared to all genes
Matrix[,19] <- log2(Matrix[,2] / Matrix[,8])
Matrix[,20] <- log2(Matrix[,3] / Matrix[,9])
Matrix[,21] <- log2(Matrix[,5] / Matrix[,8])
Matrix[,22] <- log2(Matrix[,6] / Matrix[,9])

## Plot the result
# Setup for plotting
par(mfcol=c(1,4))

# For each comparision, make a color vector based on P-values and plot the result
col_vector <- ifelse(Matrix[1,11] <= 0.05, "green3", ifelse(Matrix[1,12] <= 0.05,"red3","grey"))
for (i in 2:10) { col_vector <- c(col_vector,ifelse(Matrix[i,11] <= 0.05, "green3", ifelse(Matrix[i,12] <= 0.05,"red3","grey"))) }
barplot(Matrix[,19], ylim=c(0,2.5), las=1, main = "Induced genes > Induced interactions", col=col_vector)

col_vector <- ifelse(Matrix[1,15] <= 0.05, "green3", ifelse(Matrix[1,16] <= 0.05,"red3","grey"))
for (i in 2:10) { col_vector <- c(col_vector,ifelse(Matrix[i,15] <= 0.05, "green3", ifelse(Matrix[i,16] <= 0.05,"red3","grey"))) }
barplot(Matrix[,21], ylim=c(-1.5,0), las=1, main = "Repressed genes > Induced interactions", col=col_vector)

col_vector <- ifelse(Matrix[1,13] <= 0.05, "green3", ifelse(Matrix[1,14] <= 0.05,"red3","grey"))
for (i in 2:10) { col_vector <- c(col_vector,ifelse(Matrix[i,13] <= 0.05, "green3", ifelse(Matrix[i,14] <= 0.05,"red3","grey"))) }
barplot(Matrix[,20], ylim=c(-3,0), las=1, main = "Induced genes > Repressed interactions", col=col_vector)

col_vector <- ifelse(Matrix[1,17] <= 0.05, "green3", ifelse(Matrix[1,18] <= 0.05,"red3","grey"))
for (i in 2:10) { col_vector <- c(col_vector,ifelse(Matrix[i,17] <= 0.05, "green3", ifelse(Matrix[i,18] <= 0.05,"red3","grey"))) }
barplot(Matrix[,22], ylim=c(0,1.5), las=1, main = "Repressed genes > Repressed interactions", col=col_vector)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 4](../Links/Figure4.md)