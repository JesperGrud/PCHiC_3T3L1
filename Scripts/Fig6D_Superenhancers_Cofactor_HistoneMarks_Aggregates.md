```R
## Import the data
Counts <- read.delim("Data/ChIPseq/Counts/Superenhancers_Cofactor_HistoneMarks_aggregate.count")
colnames(Counts)[1] <- "Distance"
Counts <- Counts[,c(1,grep("Coverage", colnames(Counts)))]

## Split out counts for each data set for averaging and place into a list
# Cofactors
MED1 <- Counts[,c(1,grep("MED1_4h", colnames(Counts)),grep("MED1_D0", colnames(Counts)),52)]
P300 <- Counts[,c(1,grep("P300_4h", colnames(Counts)),grep("P300_D0", colnames(Counts)),52)]
SMC1 <- Counts[,c(1,grep("SMC1_4h", colnames(Counts)),grep("SMC1_D0", colnames(Counts)),52)]
HDAC2 <- Counts[,c(1,grep("HDAC2_4h", colnames(Counts)),grep("HDAC2_D0", colnames(Counts)),52)]
HDAC3 <- Counts[,c(1,grep("HDAC3_4h", colnames(Counts)),grep("HDAC3_D0", colnames(Counts)),52)]
NCoR <- Counts[,c(1,grep("NCoR_4h", colnames(Counts)),grep("NCoR_D0", colnames(Counts)),52)]

# Histone marks
H3K27ac <- Counts[,c(1,grep("H3K27ac_4h", colnames(Counts)),grep("H3K27ac_D0", colnames(Counts)),52)]
H3K4me1 <- Counts[,c(1,grep("H3K4me1_4h", colnames(Counts)),grep("H3K4me1_D0", colnames(Counts)),52)]
H3K4me2 <- Counts[,c(1,grep("H3K4me2_4h", colnames(Counts)),grep("H3K4me2_D0", colnames(Counts)),52)]

## Average replicates
MED1$MED1_4h <- rowMeans(MED1[,c(2,3)])
MED1$MED1_D0 <- rowMeans(MED1[,c(4,5)])
MED1 <- MED1[,c(1,7,8,6)]

P300$P300_4h <- rowMeans(P300[,c(2,3)])
P300$P300_D0 <- rowMeans(P300[,c(4,5)])
P300 <- P300[,c(1,7,8,6)]

SMC1$SMC1_4h <- rowMeans(SMC1[,c(2,3)])
SMC1$SMC1_D0 <- rowMeans(SMC1[,c(4,5)])
SMC1 <- SMC1[,c(1,7,8,6)]

HDAC2$HDAC2_4h <- rowMeans(HDAC2[,c(2,3)])
HDAC2$HDAC2_D0 <- rowMeans(HDAC2[,c(4,5)])
HDAC2 <- HDAC2[,c(1,7,8,6)]

HDAC3$HDAC3_4h <- rowMeans(HDAC3[,c(2,3)])
HDAC3$HDAC3_D0 <- rowMeans(HDAC3[,c(4,5)])
HDAC3 <- HDAC3[,c(1,7,8,6)]

NCoR$NCoR_4h <- rowMeans(NCoR[,c(2,3)])
NCoR$NCoR_D0 <- rowMeans(NCoR[,c(4,5)])
NCoR <- NCoR[,c(1,7,8,6)]

H3K27ac$H3K27ac_4h <- rowMeans(H3K27ac[,c(2,3)])
H3K27ac$H3K27ac_D0 <- rowMeans(H3K27ac[,c(4,5)])
H3K27ac <- H3K27ac[,c(1,7,8,6)]

H3K4me1$H3K4me1_4h <- rowMeans(H3K4me1[,c(2,3)])
H3K4me1$H3K4me1_D0 <- rowMeans(H3K4me1[,c(4,5)])
H3K4me1 <- H3K4me1[,c(1,7,8,6)]

H3K4me2$H3K4me2_4h <- rowMeans(H3K4me2[,c(2,3)])
H3K4me2$H3K4me2_D0 <- rowMeans(H3K4me2[,c(4,5)])
H3K4me2 <- H3K4me2[,c(1,7,8,6)]

# Place everything into a list
Counts <- list()
Counts$MED1 <- MED1
Counts$HDAC2 <- HDAC2
Counts$H3K4me1 <- H3K4me1
Counts$P300 <- P300
Counts$HDAC3 <- HDAC3
Counts$H3K4me2 <- H3K4me2
Counts$SMC1 <- SMC1
Counts$NCoR <- NCoR
Counts$H3K27ac <- H3K27ac

## Plot the stuff
par(mfcol=c(3,3))
for (i in 1:9) {
plot(-10,-10,pch='',xlim=c(0,1), ylim=c(0,max(Counts[[i]][c(2,3)])), ylab="Normalized tag count", xlab="Fraction", main=names(Counts[i]), las=1, xaxt="n")
axis(1, at=seq(0,1,by=0.25), las=2)
lines(Counts[[i]][,1], Counts[[i]][,4], col="grey", lwd=2)
lines(Counts[[i]][,1], Counts[[i]][,2], col="darkblue", lwd=2)
lines(Counts[[i]][,1], Counts[[i]][,3], col="darkred", lwd=2)
}
```

[Back to start](../README.md)<br>
[Back to overview of Figure 6](../Links/Figure6.md)