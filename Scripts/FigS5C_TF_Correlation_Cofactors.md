```R
# Load the required packages
library(pheatmap)

# Set factor to process
Factors <- c("CEBPb","CEBPd","KLF4","Jun","RXR")

# Define result matrix
CorMatrix <- data.frame(matrix(ncol=5, nrow=6))

# Loop-de-loop
for (i in 1:5) {

# Grab the name of the factor
Factor <- Factors[i]

# Load data
Point <- read.delim(paste("Data/ChIPseq/Counts/",Factor,".point", sep=""))
colnames(Point)[1] <- "PeakID"
colnames(Point)[8] <- paste(Factor,"_4h",sep="")
colnames(Point)[9] <- paste(Factor,"_D0",sep="")
Broad <- read.delim(paste("Data/ChIPseq/Counts/",Factor,".broad", sep=""))
colnames(Broad)[1] <- "PeakID"

# Process cofactor data
Point$HDAC2_4h <- rowMeans(Point[,c(10,11)])
Point$HDAC2_D0 <- rowMeans(Point[,c(12,13)])
Point$HDAC3_4h <- rowMeans(Point[,c(14,15)])
Point$HDAC3_D0 <- rowMeans(Point[,c(16,17)])
Point$NCoR_4h <- rowMeans(Point[,c(18,19)])
Point$NCoR_D0 <- rowMeans(Point[,c(20,21)])
Point$MED1_4h <- rowMeans(Point[,c(22,23)])
Point$MED1_D0 <- rowMeans(Point[,c(24,25)])
Point$P300_4h <- rowMeans(Point[,c(26,27)])
Point$P300_D0 <- rowMeans(Point[,c(28,29)])
Broad$H3K27ac_4h <- rowMeans(Broad[,c(8,9)])
Broad$H3K27ac_D0 <- rowMeans(Broad[,c(10,11)])

# Merge the data and clean up
Counts <- merge(Point[,c(1,8,9,30:39)], Broad[,c(1,12,13)], by="PeakID")
rm(Point)
rm(Broad)

# Calculate log2 fold changes
Counts$logFC <- log2(((Counts[,2]+0.25) / (Counts[,3]+0.25)))
Counts$HDAC2_logFC <- log2(((Counts[,4]+0.25) / (Counts[,5]+0.25)))
Counts$HDAC3_logFC <- log2(((Counts[,6]+0.25) / (Counts[,7]+0.25)))
Counts$NCoR_logFC <- log2(((Counts[,8]+0.25) / (Counts[,9]+0.25)))
Counts$MED1_logFC <- log2(((Counts[,10]+0.25) / (Counts[,11]+0.25)))
Counts$P300_logFC <- log2(((Counts[,12]+0.25) / (Counts[,13]+0.25)))
Counts$H3K27ac_logFC <- log2(((Counts[,14]+0.25) / (Counts[,15]+0.25)))

# Get the max occupancy for each factor across time points
CorMatrix[1,i] <- cor(Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,4] >= 20 | Counts[,5] >= 20),"logFC"], Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,4] >= 20 | Counts[,5] >= 20),"HDAC2_logFC"])
CorMatrix[2,i] <- cor(Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,6] >= 20 | Counts[,7] >= 20),"logFC"], Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,6] >= 20 | Counts[,7] >= 20),"HDAC3_logFC"])
CorMatrix[3,i] <- cor(Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,8] >= 20 | Counts[,9] >= 20),"logFC"], Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,8] >= 20 | Counts[,9] >= 20),"NCoR_logFC"])
CorMatrix[4,i] <- cor(Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,10] >= 20 | Counts[,11] >= 20),"logFC"], Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,10] >= 20 | Counts[,11] >= 20),"MED1_logFC"])
CorMatrix[5,i] <- cor(Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,12] >= 20 | Counts[,13] >= 20),"logFC"], Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,12] >= 20 | Counts[,13] >= 20),"P300_logFC"])
CorMatrix[6,i] <- cor(Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,14] >= 20 | Counts[,15] >= 20),"logFC"], Counts[ (Counts[,2] >= 20 | Counts[,3] >= 20) & (Counts[,14] >= 20 | Counts[,15] >= 20),"H3K27ac_logFC"])

# Set the colnames
colnames(CorMatrix)[i] <- Factor
}

# Set the rownames and plot the heatmap
rownames(CorMatrix) <- c("HDAC2","HDAC3","NCoR","MED1","P300","H3K27ac")

# Plot it
barplot(as.matrix(CorMatrix), beside=T, legend.text = T, col = c("darkred","red","pink","darkgreen","green","lightgreen"), las=1, ylim=c(0,1), ylab="Pearson's correlation coefficient", xlab="Transcription factor")
```

[Back to start](../README.md)<br>
[Back to overview of Figure S5](../Links/FigureS5.md)