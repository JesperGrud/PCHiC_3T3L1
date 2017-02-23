```R
## Load the necessary packages
library(DESeq2)
library(org.Mm.eg.db)

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

# Extract mapping to official symbols and entrez ids
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

# Keep only interactions between two protein coding genes that are measured in RNA-seq 
Genes <- read.delim("Data/Interactions/PCHiC/Annotation/Genes_Fragments.txt", header=T)

# Filter gene list to contain only protein coding genes measured in the RNA-seq data
Genes <- Genes[ Genes$Type == "protein_coding",]
Genes <- Genes[ Genes$Symbol %in% RNA$Symbol,]

# Filter away fragments with multiple genes
Genes$Count <- 1
Count <- aggregate(Genes$Count, by=list(Genes$Fragment), FUN="sum")
Genes <- Genes[ Genes$Fragment %in% Count[ Count$x == 1, "Group.1"],]

# Filter interactions to include only interactions between the filtered genes
FilteredInteractions <- FilteredInteractions[ FilteredInteractions$baitID %in% Genes$Fragment,] 
FilteredInteractions <- FilteredInteractions[ FilteredInteractions$oeID %in% Genes$Fragment,] 

# Merge symbols
colnames(Genes)[c(2,3)] <- c("oeID","oeSymbol")
FilteredInteractions <- merge(FilteredInteractions, Genes[,c(2,3)], by="oeID")
colnames(Genes)[c(2,3)] <- c("baitID","baitSymbol")
FilteredInteractions <- merge(FilteredInteractions, Genes[,c(2,3)], by="baitID")

# Remove duplicated connections and reduce the size of the list
FilteredInteractions <- FilteredInteractions[,c(29,30)]
FilteredInteractions$Collaps <- paste(FilteredInteractions$oeSymbol, FilteredInteractions$baitSymbol, sep="-")
FilteredInteractions <- FilteredInteractions[ duplicated(FilteredInteractions$Collaps)==F,c(1,2)]

# Define groups (induced = 1, repressed = -1) for both bait and otherEnd gene
FilteredInteractions$oeGroup <- 0
FilteredInteractions[ FilteredInteractions$oeSymbol %in% RNA[ RNA$Padj_Four <= 0.05 & log2(RNA$H4/RNA$D0) > 0 ,"Symbol"],"oeGroup"] <- 1
FilteredInteractions[ FilteredInteractions$oeSymbol %in% RNA[ RNA$Padj_Four <= 0.05 & log2(RNA$H4/RNA$D0) < 0 ,"Symbol"],"oeGroup"] <- -1
FilteredInteractions$baitGroup <- 0
FilteredInteractions[ FilteredInteractions$baitSymbol %in% RNA[ RNA$Padj_Four <= 0.05 & log2(RNA$H4/RNA$D0) > 0 ,"Symbol"],"baitGroup"] <- 1
FilteredInteractions[ FilteredInteractions$baitSymbol %in% RNA[ RNA$Padj_Four <= 0.05 & log2(RNA$H4/RNA$D0) < 0 ,"Symbol"],"baitGroup"] <- -1

# Overlap
Equi <- (nrow(FilteredInteractions[ FilteredInteractions$baitGroup == 1 & FilteredInteractions$oeGroup == 1,])+nrow(FilteredInteractions[ FilteredInteractions$baitGroup == -1 & FilteredInteractions$oeGroup == -1,]))/nrow(FilteredInteractions[ FilteredInteractions$oeGroup != 0 | FilteredInteractions$baitGroup != 0,])
Anti <- (nrow(FilteredInteractions[ (FilteredInteractions$baitGroup == 1 & FilteredInteractions$oeGroup == -1),])+nrow(FilteredInteractions[ (FilteredInteractions$baitGroup == -1 & FilteredInteractions$oeGroup == 1),]))/nrow(FilteredInteractions[ FilteredInteractions$oeGroup != 0 | FilteredInteractions$baitGroup != 0,])

# Generate random 4k random pairs 1000 times
Result <- as.data.frame(matrix(ncol=2, nrow=1000))
for (q in 1:1000) {
  Random <- as.data.frame(matrix(ncol=2, nrow=4000))
  for (i in 1:4000) { Random[i,1] <- RNA[ sample(1:nrow(RNA),1),"Symbol"] }
  for (i in 1:4000) { Random[i,2] <- RNA[ sample(1:nrow(RNA),1),"Symbol"] }
  Random <- Random[ Random$V1 != Random$V2, ]
  Random$oeGroup <- 0
  Random[ Random$V1 %in% RNA[ RNA$Padj_Four <= 0.05 & log2(RNA$H4/RNA$D0) > 0 ,"Symbol"],"oeGroup"] <- 1
  Random[ Random$V1 %in% RNA[ RNA$Padj_Four <= 0.05 & log2(RNA$H4/RNA$D0) < 0 ,"Symbol"],"oeGroup"] <- -1
  Random$baitGroup <- 0
  Random[ Random$V2 %in% RNA[ RNA$Padj_Four <= 0.05 & log2(RNA$H4/RNA$D0) > 0 ,"Symbol"],"baitGroup"] <- 1
  Random[ Random$V2 %in% RNA[ RNA$Padj_Four <= 0.05 & log2(RNA$H4/RNA$D0) < 0 ,"Symbol"],"baitGroup"] <- -1
  Result[q,1] <- (nrow(Random[ Random$baitGroup == 1 & Random$oeGroup == 1,])+nrow(Random[ Random$baitGroup == -1 & Random$oeGroup == -1,]))/nrow(Random[ Random$oeGroup != 0 | Random$baitGroup != 0,])
  Result[q,2] <- (nrow(Random[ (Random$baitGroup == 1 & Random$oeGroup == -1),])+nrow(Random[ (Random$baitGroup == -1 & Random$oeGroup == 1),]))/nrow(Random[ Random$oeGroup != 0 | Random$baitGroup != 0,])
  print(q)
}

# Plot it
B <- barplot(c(Equi, mean(Result[,1]), Anti,mean(Result[,2])), las=1, ylim=c(0,0.15), col=c("green","grey","red","grey"))
arrows(as.vector(B), c(0,mean(Result[,1])-sd(Result[,1]),0,mean(Result[,2])-sd(Result[,2])), as.vector(B), c(0,mean(Result[,1])+sd(Result[,1]),0,mean(Result[,2])+sd(Result[,2])), lwd = 1.5, angle = 90, code = 3, length = 0.05)
```

[Back to start](../README.md)<br>
[Back to overview of Figure S2](../Links/FigureS2.md)