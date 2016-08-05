```R
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

# Define variable for counting each type for region
FilteredInteractions$Count <- 1

# Count the number of interactions from each non-baited fragment
Nonbaited <- aggregate(x = FilteredInteractions$Count, by=list(FilteredInteractions$oeID), FUN="sum")

# Calculate the percent of nonbaited fragments with 1, 2, 3, 4 and 5 more or connections
PercentNonbaited <- as.matrix(c(
	(nrow(Nonbaited[ Nonbaited$x == 1,])/nrow(Nonbaited))*100,
	(nrow(Nonbaited[ Nonbaited$x == 2,])/nrow(Nonbaited))*100,
	(nrow(Nonbaited[ Nonbaited$x == 3,])/nrow(Nonbaited))*100,
	(nrow(Nonbaited[ Nonbaited$x == 4,])/nrow(Nonbaited))*100,
	(nrow(Nonbaited[ Nonbaited$x >= 5,])/nrow(Nonbaited))*100))

# Plot a barplot for the non-baited fragments
barplot(PercentNonbaited, beside = T, las = 1, names=c("1","2","5","4",">= 5"), 
col = colorRampPalette(c("black","blue"))(5), ylab = "% non-baited fragments", xlab="# connections")
```

[Back to start](../README.md)<br>
[Back to overview of Figure 1](../Links/Figure1.md)