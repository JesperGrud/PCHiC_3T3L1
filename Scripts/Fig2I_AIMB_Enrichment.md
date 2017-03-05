```R
## Load packages
library(edgeR)
library(org.Mm.eg.db)

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

## Import and process allele-specificity results
Genes <- read.delim("Data/Other/CHT_GeneExpression.result.txt", header=T)
Genes <- Genes[ duplicated(paste(Genes$Chr, Genes$Pos,sep="-"))==F,]
PCHiC <- read.delim("Data/Other/CHT_PCHiC.result.txt", header=F)

## Cleanup
rm(list=setdiff(ls(),c("FilteredInteractions_Genes","Genes","PCHiC")))

## Calculate heterozygocity
Genes$Input_Genotype <- abs((Genes$Input_Genotype - 1))
Genes$RNA_Genotype <- abs((Genes$RNA_Genotype - 1))
PCHiC$V10 <- abs((PCHiC$V10 - 1))
PCHiC$V9 <- abs((PCHiC$V9 - 1))


## Keep only tests of interest for genes
Genes <- Genes[ order(Genes$RefSeq, Genes$Pval),]
Genes <- Genes[ duplicated(Genes$RefSeq)==F,]
Genes <- Genes[ apply(Genes[,c(5:8)], 1, FUN="var") >= 30,]

## Add genes information
Convert <- suppressMessages(select(org.Mm.eg.db, as.character(Genes$RefSeq), columns = c("SYMBOL","REFSEQ"), "REFSEQ"))
colnames(Convert) <- c("RefSeq","Symbol")
Genes <- merge(Genes, Convert, by="RefSeq")
Genes <- unique(Genes)
rm(Convert)

## Add interaction information
PCHiC <- merge(PCHiC, FilteredInteractions_Genes[,c("ID","dist","Symbol")], by.x="V6",by.y="ID")

## Keep only tests of interest for interactions
PCHiC <- PCHiC[ apply(PCHiC[,c(5:8)], 1, FUN="var") >= 30,]

## FDR correct
PCHiC$V3 <- fdrtool(PCHiC$V3, statistic = "pvalue")$qval
Genes$Pval <- fdrtool(Genes$Pval, statistic = "pvalue")$qval

## Group the genes by AIMB
GenesImbalanced <- Genes[ Genes$Pval <= 0.1 ,]
PCHICImbalanced <- PCHiC[ PCHiC$V3 <= 0.1 ,]

## Collaps interactions to with multiple SNPs single interactions (keep most significant for imbalanced)
PCHICImbalanced$snpID <- paste(PCHICImbalanced[,1], PCHICImbalanced$Symbol, sep="-")
PCHICImbalanced <- PCHICImbalanced[ order(PCHICImbalanced$snpID, PCHICImbalanced$V3),]
PCHICImbalanced <- PCHICImbalanced[ duplicated(PCHICImbalanced$snpID)==F,]

## Collaps SNPs that connect multiple times to same gene
PCHICImbalanced$snpID <- paste(PCHICImbalanced[,3], PCHICImbalanced$Symbol, sep="-")
PCHICImbalanced <- PCHICImbalanced[ order(PCHICImbalanced$snpID, PCHICImbalanced$V3),]
PCHICImbalanced <- PCHICImbalanced[ duplicated(PCHICImbalanced$snpID)==F,]

## Random sampling for all interactions and balanced non-significant interactions (100x)
Results <- data.frame(matrix(ncol=1, nrow=100))
for (q in 1:100) {
	# Make a data.frame to save distance matched interactions
	RandomAll <- data.frame(matrix(ncol=2,nrow(PCHICImbalanced)))
	
	# For each imbalanced interaction, grab a random distance matched one (+/- 1%)
	for (i in 1:nrow(PCHICImbalanced)) {
		Tmp <- FilteredInteractions_Genes[ FilteredInteractions_Genes$Zero_Score >= 5,]
		Tmp <- Tmp[ abs(Tmp$dist) <= abs(PCHICImbalanced[i,"dist"])+(abs(PCHICImbalanced[i,"dist"])*0.01),]
		Tmp <- Tmp[ abs(Tmp$dist) >= abs(PCHICImbalanced[i,"dist"])-(abs(PCHICImbalanced[i,"dist"])*0.01),]
		Tmp <- Tmp[ sign(Tmp$dist) == sign(PCHICImbalanced[i,"dist"]),]
		Tmp <- Tmp[ !(Tmp$ID %in% PCHICImbalanced[i,"V6"]), ]
		Tmp <- Tmp[ sample(1:nrow(Tmp),1),]
		RandomAll[i,1] <- as.character(Tmp[1,"Symbol"])
		RandomAll[i,2] <- Tmp[1,"dist"]
	}
	
	# Save the fraction of connected to imbalanced genes
	Results[q,1] <- nrow(RandomAll[ RandomAll[,1] %in% GenesImbalanced$Symbol,])/nrow(RandomAll)
	print(q)
}

# Calcualte the results from the imbalanced interactions
AIMB <- nrow(PCHICImbalanced[ PCHICImbalanced$Symbol %in% GenesImbalanced$Symbol,])/nrow(PCHICImbalanced)

# Plot the results
B <- barplot(c(AIMB, mean(Results[,1])), beside = T, las=1, ylim=c(0,0.15), col=c("green3","grey"), ylab="Fraction connected to imbalanced gene")
arrows(B, c(0,mean(Results[,1])+sd(Results[,1])), B,c(0,mean(Results[,1])-sd(Results[,1])), lwd = 1.5, angle = 90, code = 3, length = 0.05)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 2](../Links/Figure2.md)
