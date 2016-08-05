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

## Loop through bins of coactivator and corepressor composite scores and calculate the average log2 fold of coactivators
# Make a data frame to capture the results
Heatmap_Coactivators <- data.frame(matrix(ncol=length(seq(-12.4,12.4, by=0.4))-1,nrow=length(seq(-12.4,12.4, by=0.4))-1))

# Loop through the data and count the number of sites
for (i in 1:(length(seq(-12.4,12.4, by=0.4))-1)) {
  Coact_Low <- seq(-12.4,12.4, by=0.4)[i]
  Coact_High <- seq(-12.4,12.4, by=0.4)[(i+1)]
  for (q in 1:(length(seq(-12.4,12.4, by=0.4))-1) ) {
    Corep_Low <- seq(-12.4,12.4, by=0.4)[q]
    Corep_High <- seq(-12.4,12.4, by=0.4)[(q+1)]
    Tmp <- Cofactors[ Cofactors$Coactivators >= Coact_Low & Cofactors$Coactivators < Coact_High & Cofactors$Corepressors >= Corep_Low & Cofactors$Corepressors < Corep_High,]
    if (nrow(Tmp) == 0) { Heatmap_Coactivators[i,q] <- NA } else { Heatmap_Coactivators[i,q] <- mean(Tmp$Average_Coactivator) }
  }
}

## Loop through bins of coactivator and corepressor composite scores and calculate the average log2 fold of corepressors
# Make a data frame to capture the results
Heatmap_Corepressors <- data.frame(matrix(ncol=length(seq(-12.4,12.4, by=0.4))-1,nrow=length(seq(-12.4,12.4, by=0.4))-1))

# Loop through the data and count the number of sites
for (i in 1:(length(seq(-12.4,12.4, by=0.4))-1)) {
  Coact_Low <- seq(-12.4,12.4, by=0.4)[i]
  Coact_High <- seq(-12.4,12.4, by=0.4)[(i+1)]
  for (q in 1:(length(seq(-12.4,12.4, by=0.4))-1) ) {
    Corep_Low <- seq(-12.4,12.4, by=0.4)[q]
    Corep_High <- seq(-12.4,12.4, by=0.4)[(q+1)]
    Tmp <- Cofactors[ Cofactors$Coactivators >= Coact_Low & Cofactors$Coactivators < Coact_High & Cofactors$Corepressors >= Corep_Low & Cofactors$Corepressors < Corep_High,]
    if (nrow(Tmp) == 0) { Heatmap_Corepressors[i,q] <- NA } else { Heatmap_Corepressors[i,q] <- mean(Tmp$Average_Corepressor) }
  }
}

## Setup for plotting the heatmap
par(mfcol=c(1,3))
Labels <- signif(seq(-12.4,12.4, by=0.4),3)
Breaks <- seq(-6,6,length.out = 61)
Colors <- colorRampPalette(c("red","white","green"))(60)

# Plot the heatmap of coactivators
#par(fig=c(0,0.8,0,1))
plot(Labels, Labels, type="n", axes=T, xlab="Coactivators", ylab="Corepressors",,xaxs="i", yaxs="i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
image(Labels, Labels, as.matrix(Heatmap_Coactivators), breaks=Breaks,col=Colors, add=T)

# Plot the heatmap of coactivators
plot(Labels, Labels, type="n", axes=T, xlab="Coactivators", ylab="Corepressors",,xaxs="i", yaxs="i")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
image(Labels, Labels, as.matrix(Heatmap_Corepressors), breaks=Breaks,col=Colors, add=T)

# Plot the legend
#par(fig=c(0.8,1,0,1), new=t)
plot(1,1,t="n",ylim=c(-6,6), xlim=c(0,1), xaxt="n", yaxt="s", xlab="", ylab="Mean cofactor binding intensity",xaxs="i", yaxs="i", las = 1, frame.plot=F) 
for (i in 1:(length(Breaks)-1)) {
	rect(0,Breaks[[i]],1,Breaks[[(i+1)]], col=Colors[[i]], border=NA)
}
box()
```

[Back to start](../README.md)<br>
[Back to overview of Figure S5](../Links/FigureS5.md)
