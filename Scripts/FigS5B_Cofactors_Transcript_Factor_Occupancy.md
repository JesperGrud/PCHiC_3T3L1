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

## Process TF counts
# Import data and set colnames
Counts <- read.delim("Data/ChIPseq/Counts/Cofactors_TFs.count", header=T)

# Average replicates for CTCF
Counts$CTCF_4h <- rowMeans(Counts[,c("CTCF_4h_1","CTCF_4h_2")])
Counts$CTCF_D0 <- rowMeans(Counts[,c("CTCF_D0_1","CTCF_D0_2")])

# Calculate log2 fold changes using 0.25 as a pseudocount
Counts$CEBPb <- log2((Counts$CEBPb_4h + 0.25) / (Counts$CEBPb_D0 + 0.25))
Counts$CEBPd <- log2((Counts$CEBPd_4h + 0.25) / (Counts$CEBPd_D0 + 0.25))
Counts$Jun <- log2((Counts$Jun_4h + 0.25) / (Counts$Jun_D0 + 0.25))
Counts$KLF4 <- log2((Counts$KLF_4h + 0.25) / (Counts$KLF_D0 + 0.25))
Counts$RXR <- log2((Counts$RXR_4h + 0.25) / (Counts$RXR_D0 + 0.25))
Counts$CTCF <- log2((Counts$CTCF_4h + 0.25) / (Counts$CTCF_D0 + 0.25))

## Loop through bins of coactivator and corepressor composite scores and calculate the average log2 fold of TFs
# Make a data frame for each TF to capture the results
Heatmap_CEBPb <- data.frame(matrix(ncol=length(seq(-12.4,12.4, by=0.4))-1,nrow=length(seq(-12.4,12.4, by=0.4))-1))
Heatmap_CEBPd <- data.frame(matrix(ncol=length(seq(-12.4,12.4, by=0.4))-1,nrow=length(seq(-12.4,12.4, by=0.4))-1))
Heatmap_Jun <- data.frame(matrix(ncol=length(seq(-12.4,12.4, by=0.4))-1,nrow=length(seq(-12.4,12.4, by=0.4))-1))
Heatmap_KLF4 <- data.frame(matrix(ncol=length(seq(-12.4,12.4, by=0.4))-1,nrow=length(seq(-12.4,12.4, by=0.4))-1))
Heatmap_RXR <- data.frame(matrix(ncol=length(seq(-12.4,12.4, by=0.4))-1,nrow=length(seq(-12.4,12.4, by=0.4))-1))
Heatmap_CTCF <- data.frame(matrix(ncol=length(seq(-12.4,12.4, by=0.4))-1,nrow=length(seq(-12.4,12.4, by=0.4))-1))

# Loop through the data and count the number of sites
for (i in 1:(length(seq(-12.4,12.4, by=0.4))-1)) {
  Coact_Low <- seq(-12.4,12.4, by=0.4)[i]
  Coact_High <- seq(-12.4,12.4, by=0.4)[(i+1)]
  for (q in 1:(length(seq(-12.4,12.4, by=0.4))-1) ) {
    Corep_Low <- seq(-12.4,12.4, by=0.4)[q]
    Corep_High <- seq(-12.4,12.4, by=0.4)[(q+1)]
    Tmp <- Cofactors[ Cofactors$Coactivators >= Coact_Low & Cofactors$Coactivators < Coact_High & Cofactors$Corepressors >= Corep_Low & Cofactors$Corepressors < Corep_High,]
	Tmp2 <- Counts[ Counts$PeakID %in% Tmp$PeakID,]
    if (nrow(Tmp) == 0) { 
		Heatmap_CEBPb[i,q] <- NA 
		Heatmap_CEBPd[i,q] <- NA 
		Heatmap_Jun[i,q] <- NA 
		Heatmap_KLF4[i,q] <- NA 
		Heatmap_RXR[i,q] <- NA 
		Heatmap_CTCF[i,q] <- NA 
		} else { 
		Heatmap_CEBPb[i,q] <- mean(Tmp2$CEBPb)
		Heatmap_CEBPd[i,q] <- mean(Tmp2$CEBPd)
		Heatmap_Jun[i,q] <- mean(Tmp2$Jun) 
		Heatmap_KLF4[i,q] <- mean(Tmp2$KLF4) 
		Heatmap_RXR[i,q] <- mean(Tmp2$RXR) 
		Heatmap_CTCF[i,q] <- mean(Tmp2$CTCF) 
		}
  }
}

# Setup for plotting the heatmap
Labels <- signif(seq(-12.4,12.4, by=0.4),3)
Breaks <- seq(-6,6,length.out = 61)
Colors <- colorRampPalette(c("red","white","green"))(60)

# Set the area
par(mfcol=c(2,4))

# Plot the heatmaps
plot(Labels, Labels, type="n", axes=T, xlab="Coactivators", ylab="Corepressors",,xaxs="i", yaxs="i", main="CEBPb")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
image(Labels, Labels, as.matrix(Heatmap_CEBPb), breaks=Breaks,col=Colors, add=T)

plot(Labels, Labels, type="n", axes=T, xlab="Coactivators", ylab="Corepressors",,xaxs="i", yaxs="i", main="CEBPd")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
image(Labels, Labels, as.matrix(Heatmap_CEBPd), breaks=Breaks,col=Colors, add=T)

plot(Labels, Labels, type="n", axes=T, xlab="Coactivators", ylab="Corepressors",,xaxs="i", yaxs="i", main="cJun")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
image(Labels, Labels, as.matrix(Heatmap_Jun), breaks=Breaks,col=Colors, add=T)

plot(Labels, Labels, type="n", axes=T, xlab="Coactivators", ylab="Corepressors",,xaxs="i", yaxs="i", main="KLF4")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
image(Labels, Labels, as.matrix(Heatmap_KLF4), breaks=Breaks,col=Colors, add=T)

plot(Labels, Labels, type="n", axes=T, xlab="Coactivators", ylab="Corepressors",,xaxs="i", yaxs="i", main="RXR")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
image(Labels, Labels, as.matrix(Heatmap_RXR), breaks=Breaks,col=Colors, add=T)

plot(Labels, Labels, type="n", axes=T, xlab="Coactivators", ylab="Corepressors",,xaxs="i", yaxs="i", main="CTCF")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
image(Labels, Labels, as.matrix(Heatmap_CTCF), breaks=Breaks,col=Colors, add=T)

# Plot the legend
plot(1,1,t="n",ylim=c(min(Breaks),max(Breaks)), xlim=c(0,1), xaxt="n", yaxt="s", xlab="", ylab="Mean corepressor binding intensity",xaxs="i", yaxs="i", las = 1, frame.plot=F) 
for (i in 1:(length(Breaks)-1)) {
	rect(0,Breaks[[i]],1,Breaks[[(i+1)]], col=Colors[[i]], border=NA)
}
box()
```

[Back to start](../README.md)<br>
[Back to overview of Figure S5](../Links/FigureS5.md)