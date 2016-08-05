```R
# Load libraries
library(DESeq2)
library(org.Mm.eg.db)
library(clusterProfiler)

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

# Rename clusters for intuitive grouping during plots
RNA[ RNA$Cluster == 2, "Cluster"] <- 8
RNA[ RNA$Cluster == 4, "Cluster"] <- 9
RNA[ RNA$Cluster == 6, "Cluster"] <- 10
RNA[ RNA$Cluster == 7, "Cluster"] <- 11
RNA[ RNA$Cluster == 8, "Cluster"] <- 4
RNA[ RNA$Cluster == 9, "Cluster"] <- 2
RNA[ RNA$Cluster == 10, "Cluster"] <- 7
RNA[ RNA$Cluster == 11, "Cluster"] <- 6

# Import transcript length information and merge it
Lengths <- read.delim("Data/RNAseq/Annotation/Gene.lengths", header=T)
RNA <- merge(RNA, Lengths, by="RefSeq")

# Extract mapping to official symbols and entrez ids
Convert <- suppressMessages(select(org.Mm.eg.db, as.character(RNA$RefSeq), columns = c("SYMBOL","REFSEQ"), "REFSEQ"))
colnames(Convert) <- c("RefSeq","Symbol")
RNA <- merge(RNA, Convert, by="RefSeq")

## Perform GO analysis
# Get Entrez IDs for each cluster and place into a list
Clusters <- list()
Clusters$C1 <- bitr(RNA[ RNA$Cluster == 1, "RefSeq"], "REFSEQ", "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
Clusters$C2 <- bitr(RNA[ RNA$Cluster == 2, "RefSeq"], "REFSEQ", "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
Clusters$C3 <- bitr(RNA[ RNA$Cluster == 3, "RefSeq"], "REFSEQ", "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
Clusters$C4 <- bitr(RNA[ RNA$Cluster == 4, "RefSeq"], "REFSEQ", "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
Clusters$C5 <- bitr(RNA[ RNA$Cluster == 5, "RefSeq"], "REFSEQ", "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
Clusters$C6 <- bitr(RNA[ RNA$Cluster == 6, "RefSeq"], "REFSEQ", "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
Clusters$C7 <- bitr(RNA[ RNA$Cluster == 7, "RefSeq"], "REFSEQ", "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID

# Define the background
Background <- bitr(RNA[ , "RefSeq"], "REFSEQ", "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID

# Run GO analysis
CGO <- compareCluster(Clusters, fun = "enrichGO", OrgDb = org.Mm.eg.db, universe = Background, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, pvalueCutoff = 0.05, ont = "BP")

# Simplify the results
CGO_Simple <- simplify(CGO)

# Plot the results
dotplot(CGO_Simple)
```

[Back to start](../README.md)<br>
[Back to overview of Figure S2](../Links/FigureS2.md)