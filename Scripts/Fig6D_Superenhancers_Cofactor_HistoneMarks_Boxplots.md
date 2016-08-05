```R
## Import the data
Counts <- read.delim("Data/ChIPseq/Counts/Superenhancers_Cofactor_HistoneMarks_window.count")
colnames(Counts)[1] <- "PeakID"

## Split out counts for each data set for averaging and place into a list
# Cofactors
MED1 <- Counts[,c(grep("MED1_4h", colnames(Counts)),grep("MED1_D0", colnames(Counts)),58)]
P300 <- Counts[,c(grep("P300_4h", colnames(Counts)),grep("P300_D0", colnames(Counts)),58)]
SMC1 <- Counts[,c(grep("SMC1_4h", colnames(Counts)),grep("SMC1_D0", colnames(Counts)),58)]
HDAC2 <- Counts[,c(grep("HDAC2_4h", colnames(Counts)),grep("HDAC2_D0", colnames(Counts)),58)]
HDAC3 <- Counts[,c(grep("HDAC3_4h", colnames(Counts)),grep("HDAC3_D0", colnames(Counts)),58)]
NCoR <- Counts[,c(grep("NCoR_4h", colnames(Counts)),grep("NCoR_D0", colnames(Counts)),58)]

# Histone marks
H3K27ac <- Counts[,c(grep("H3K27ac_4h", colnames(Counts)),grep("H3K27ac_D0", colnames(Counts)),58)]
H3K4me1 <- Counts[,c(grep("H3K4me1_4h", colnames(Counts)),grep("H3K4me1_D0", colnames(Counts)),58)]
H3K4me2 <- Counts[,c(grep("H3K4me2_4h", colnames(Counts)),grep("H3K4me2_D0", colnames(Counts)),58)]

## Average replicates
MED1$MED1_4h <- rowMeans(MED1[,c(1,2)])
MED1$MED1_D0 <- rowMeans(MED1[,c(3,4)])
MED1 <- MED1[,c(6,7,5)]

P300$P300_4h <- rowMeans(P300[,c(1,2)])
P300$P300_D0 <- rowMeans(P300[,c(3,4)])
P300 <- P300[,c(6,7,5)]

SMC1$SMC1_4h <- rowMeans(SMC1[,c(1,2)])
SMC1$SMC1_D0 <- rowMeans(SMC1[,c(3,4)])
SMC1 <- SMC1[,c(6,7,5)]

HDAC2$HDAC2_4h <- rowMeans(HDAC2[,c(1,2)])
HDAC2$HDAC2_D0 <- rowMeans(HDAC2[,c(3,4)])
HDAC2 <- HDAC2[,c(6,7,5)]

HDAC3$HDAC3_4h <- rowMeans(HDAC3[,c(1,2)])
HDAC3$HDAC3_D0 <- rowMeans(HDAC3[,c(3,4)])
HDAC3 <- HDAC3[,c(6,7,5)]

NCoR$NCoR_4h <- rowMeans(NCoR[,c(1,2)])
NCoR$NCoR_D0 <- rowMeans(NCoR[,c(3,4)])
NCoR <- NCoR[,c(6,7,5)]

H3K27ac$H3K27ac_4h <- rowMeans(H3K27ac[,c(1,2)])
H3K27ac$H3K27ac_D0 <- rowMeans(H3K27ac[,c(3,4)])
H3K27ac <- H3K27ac[,c(6,7,5)]

H3K4me1$H3K4me1_4h <- rowMeans(H3K4me1[,c(1,2)])
H3K4me1$H3K4me1_D0 <- rowMeans(H3K4me1[,c(3,4)])
H3K4me1 <- H3K4me1[,c(6,7,5)]

H3K4me2$H3K4me2_4h <- rowMeans(H3K4me2[,c(1,2)])
H3K4me2$H3K4me2_D0 <- rowMeans(H3K4me2[,c(3,4)])
H3K4me2 <- H3K4me2[,c(6,7,5)]

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

## Plot it all
par(mfcol=c(3,3))
for(i in 1:9) {
  boxplot(Counts[[i]][,2], Counts[[i]][,1], outline=F, col=c("darkblue", "darkred","grey"), ylab="ChIP-seq signal", names=c("D0", "4h"), xlab="Time point", main=names(Counts[i]), las=1)
}
```

[Back to start](../README.md)<br>
[Back to overview of Figure 6](../Links/Figure6.md)