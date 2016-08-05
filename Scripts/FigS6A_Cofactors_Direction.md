```R
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

# Define groups based on coactivator change and remove sites that are mixed
InducedCoactivators <- Cofactors[ Cofactors$logFC_MED1 >= log2(1.5) & Cofactors$logFC_P300 >= log2(1.5) & Cofactors$logFC_SMC1 >= log2(1.5),]
ConstantCoactivators <- Cofactors[ abs(Cofactors$logFC_MED1) <= log2(1.5) & abs(Cofactors$logFC_P300) <= log2(1.5) & abs(Cofactors$logFC_SMC1) <= log2(1.5),]
RepressedCoactivators <- Cofactors[ Cofactors$logFC_MED1 <= -log2(1.5) & Cofactors$logFC_P300 <= -log2(1.5) & Cofactors$logFC_SMC1 <= -log2(1.5),]

# Define groups based on coactivator change and remove sites that are mixed
InducedCorepressors <- Cofactors[ Cofactors$logFC_HDAC2 >= log2(1.5) & Cofactors$logFC_HDAC3 >= log2(1.5) & Cofactors$logFC_NCoR >= log2(1.5),]
ConstantCorepressors <- Cofactors[ abs(Cofactors$logFC_HDAC2) <= log2(1.5) & abs(Cofactors$logFC_HDAC3) <= log2(1.5) & abs(Cofactors$logFC_NCoR) <= log2(1.5),]
RepressedCorepressors <- Cofactors[ Cofactors$logFC_HDAC2 <= -log2(1.5) & Cofactors$logFC_HDAC3 <= -log2(1.5) & Cofactors$logFC_NCoR <= -log2(1.5),]

# Make a matrix to capture the fractions
Result <- matrix(ncol=3, nrow=3)

# Count the number of sites
Result[1,1] <- nrow(InducedCoactivators[ InducedCoactivators$PeakID %in% InducedCorepressors$PeakID,])
Result[2,1] <- nrow(InducedCoactivators[ InducedCoactivators$PeakID %in% ConstantCorepressors$PeakID,])
Result[3,1] <- nrow(InducedCoactivators[ InducedCoactivators$PeakID %in% RepressedCorepressors$PeakID,])
Result[1,2] <- nrow(ConstantCoactivators[ ConstantCoactivators$PeakID %in% InducedCorepressors$PeakID,])
Result[2,2] <- nrow(ConstantCoactivators[ ConstantCoactivators$PeakID %in% ConstantCorepressors$PeakID,])
Result[3,2] <- nrow(ConstantCoactivators[ ConstantCoactivators$PeakID %in% RepressedCorepressors$PeakID,])
Result[1,3] <- nrow(RepressedCoactivators[ RepressedCoactivators$PeakID %in% InducedCorepressors$PeakID,])
Result[2,3] <- nrow(RepressedCoactivators[ RepressedCoactivators$PeakID %in% ConstantCorepressors$PeakID,])
Result[3,3] <- nrow(RepressedCoactivators[ RepressedCoactivators$PeakID %in% RepressedCorepressors$PeakID,])

# Calculate the fraction of the total
Result[,1] <- Result[,1] / sum(Result[,1])
Result[,2] <- Result[,2] / sum(Result[,2])
Result[,3] <- Result[,3] / sum(Result[,3])

## Plot the results
par(mfcol=c(1,3))
barplot(Result[,1], ylim=c(0,1), col=c("green","grey","red"), las=1, names=c("Gain","Const","Loss"), main="Gain of coactivators")
barplot(Result[,2], ylim=c(0,1), col=c("green","grey","red"), las=1, names=c("Gain","Const","Loss"), main="Constant coactivators")
barplot(Result[,3], ylim=c(0,1), col=c("green","grey","red"), las=1, names=c("Gain","Const","Loss"), main="Loss of coactivators")
```

[Back to start](../README.md)<br>
[Back to overview of Figure S6](../Links/FigureS6.md)