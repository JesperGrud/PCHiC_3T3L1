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

# Z-transform the log2 fold changes
for (i in 17:22) {
  Mean <- mean(Cofactors[,i])
  SD <- sd(Cofactors[,i])
  Cofactors[,i] <- (Cofactors[,i]-Mean)/SD
}

# Clean up the workspace
rm(list=c("Mean","SD","i"))

# Calculate composite scores
Cofactors$Coactivators <- apply(Cofactors[,c(17,18,19)],1,FUN="sum")
Cofactors$Corepressors <- apply(Cofactors[,c(20,21,22)],1,FUN="sum")

# Reset the graphics device
graphics.off()

# Plot the smoothend scatter plot of composite scores
smoothScatter(Cofactors$Coactivators, Cofactors$Corepressors, pch=16, cex=0, xlab="Coactivators", ylab="Corepressors", las=1, transformation = function(x) x^0.25)

# Add a diagnoal and regression line
abline(lm(Cofactors$Corepressors ~ Cofactors$Coactivators), col="red")
abline(0,1,lty=2)

# Add regression information
Regression <- lm(Cofactors$Corepressors ~ Cofactors$Coactivators)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 7](../Links/Figure7.md)