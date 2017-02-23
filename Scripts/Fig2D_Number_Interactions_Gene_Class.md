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

# Import gene information and split interactions into bait-to-gene and bait-to-genome
Genes <- read.delim("Data/Interactions/PCHiC/Annotation/Genes_Fragments.txt", header=T)
BtGene <- FilteredInteractions[ FilteredInteractions$oeID %in% Genes$Fragment,] 
BtGenome <- FilteredInteractions[ !(FilteredInteractions$oeID %in% Genes$Fragment),] 

# Count the number of different types of interactions 
CountPP <- nrow(BtGene) / nrow(FilteredInteractions) * 100
CountPCG <- nrow(BtGenome[ BtGenome$baitID %in% Genes[ Genes$Type == "protein_coding","Fragment"],]) / nrow(FilteredInteractions) * 100
CountNCG <- nrow(BtGenome[ BtGenome$baitID %in% Genes[ Genes$Type != "protein_coding","Fragment"],]) / nrow(FilteredInteractions) * 100

# Setup the data in a matrix for a stacked barplot
Count <- matrix(c(CountPP,0,CountPCG,CountNCG), ncol=2)

# Make the plot
par(mfcol=c(1,1))
barplot(Count, col=c("green","red"), las = 1, ylim=c(0,100),
	ylab="% interactions", names=c("P-P","P-G"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 2](../Links/Figure2.md)