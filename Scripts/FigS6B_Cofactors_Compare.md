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

# Reset the graphics device
graphics.off()

# Set up the plot area
par(mfcol=c(5,5))

# Plot the log2 fold changes of all vs all with regression information
x = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$P300_4h >= 20 | Cofactors$P300_4h >= 20),"logFC_MED1"]
y = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$P300_4h >= 20 | Cofactors$P300_4h >= 20),"logFC_P300"]
plot(x, y, xlab="MED1", ylab="P300", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_4h >= 20),"logFC_MED1"]
y = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_4h >= 20),"logFC_SMC1"]
plot(x, y, xlab="MED1", ylab="SMC1", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_4h >= 20),"logFC_MED1"]
y = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_4h >= 20),"logFC_HDAC2"]
plot(x, y, xlab="MED1", ylab="HDAC2", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_4h >= 20),"logFC_MED1"]
y = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_4h >= 20),"logFC_HDAC3"]
plot(x, y, xlab="MED1", ylab="HDAC3", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_MED1"]
y = Cofactors[(Cofactors$MED1_4h >= 20 | Cofactors$MED1_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_NCoR"]
plot(x, y, xlab="MED1", ylab="NCoR", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

plot(1, type="n", axes=F, xlab="", ylab="")

x = Cofactors[(Cofactors$P300_4h >= 20 | Cofactors$P300_D0 >= 20) & (Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_4h >= 20),"logFC_P300"]
y = Cofactors[(Cofactors$P300_4h >= 20 | Cofactors$P300_D0 >= 20) & (Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_4h >= 20),"logFC_SMC1"]
plot(x, y, xlab="P300", ylab="SMC1", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$P300_4h >= 20 | Cofactors$P300_D0 >= 20) & (Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_4h >= 20),"logFC_P300"]
y = Cofactors[(Cofactors$P300_4h >= 20 | Cofactors$P300_D0 >= 20) & (Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_4h >= 20),"logFC_HDAC2"]
plot(x, y, xlab="P300", ylab="HDAC2", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$P300_4h >= 20 | Cofactors$P300_D0 >= 20) & (Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_4h >= 20),"logFC_P300"]
y = Cofactors[(Cofactors$P300_4h >= 20 | Cofactors$P300_D0 >= 20) & (Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_4h >= 20),"logFC_HDAC3"]
plot(x, y, xlab="P300", ylab="HDAC3", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$P300_4h >= 20 | Cofactors$P300_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_P300"]
y = Cofactors[(Cofactors$P300_4h >= 20 | Cofactors$P300_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_NCoR"]
plot(x, y, xlab="P300", ylab="NCoR", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")

x = Cofactors[(Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_D0 >= 20) & (Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_4h >= 20),"logFC_SMC1"]
y = Cofactors[(Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_D0 >= 20) & (Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_4h >= 20),"logFC_HDAC2"]
plot(x, y, xlab="SMC1", ylab="HDAC2", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_D0 >= 20) & (Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_4h >= 20),"logFC_SMC1"]
y = Cofactors[(Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_D0 >= 20) & (Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_4h >= 20),"logFC_HDAC3"]
plot(x, y, xlab="SMC1", ylab="HDAC3", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_SMC1"]
y = Cofactors[(Cofactors$SMC1_4h >= 20 | Cofactors$SMC1_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_NCoR"]
plot(x, y, xlab="SMC1", ylab="NCoR", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")

x = Cofactors[(Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_D0 >= 20) & (Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_4h >= 20),"logFC_HDAC2"]
y = Cofactors[(Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_D0 >= 20) & (Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_4h >= 20),"logFC_HDAC3"]
plot(x, y, xlab="HDAC2", ylab="HDAC3", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

x = Cofactors[(Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_HDAC2"]
y = Cofactors[(Cofactors$HDAC2_4h >= 20 | Cofactors$HDAC2_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_NCoR"]
plot(x, y, xlab="HDAC2", ylab="NCoR", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))

plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")

x = Cofactors[(Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_HDAC3"]
y = Cofactors[(Cofactors$HDAC3_4h >= 20 | Cofactors$HDAC3_D0 >= 20) & (Cofactors$NCoR_4h >= 20 | Cofactors$NCoR_4h >= 20),"logFC_NCoR"]
plot(x, y, xlab="HDAC3", ylab="NCoR", xlim=c(-10,10), ylim=c(-10,10), las=1, pch=16, cex=0.5)
abline(0,1,lty=2)
Regression <- lm(y ~ x)
Pearson <- signif(sqrt(summary(Regression)$r.square),2)
Alpha <- as.numeric(signif(Regression$coefficients[2],2))
text(-8.5,10,bquote(alpha == .(Alpha)))
text(-8.5,9,bquote(R == .(Pearson)))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S6](../Links/FigureS6.md)
