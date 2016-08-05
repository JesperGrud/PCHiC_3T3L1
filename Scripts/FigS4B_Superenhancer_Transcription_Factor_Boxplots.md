```R
## Import the data
Counts <- read.delim("Data/ChIPseq/Counts/Superenhancers_Cofactor_HistoneMarks_window.count")
colnames(Counts)[1] <- "PeakID"

## Split out counts for each data set for averaging and place into a list
# TFs
CEBPb <- Counts[,c(grep("CEBPb_4h", colnames(Counts)),grep("CEBPb_D0", colnames(Counts)),58)]
CEBPd <- Counts[,c(grep("CEBPd_4h", colnames(Counts)),grep("CEBPd_D0", colnames(Counts)),58)]
Jun <- Counts[,c(grep("Jun_4h", colnames(Counts)),grep("Jun_D0", colnames(Counts)),58)]
KLF4 <- Counts[,c(grep("KLF4_4h", colnames(Counts)),grep("KLF4_D0", colnames(Counts)),58)]
RXR <- Counts[,c(grep("RXR_4h", colnames(Counts)),grep("RXR_D0", colnames(Counts)),58)]
CTCF <- Counts[,c(grep("CTCF_4h", colnames(Counts)),grep("CTCF_D0", colnames(Counts)),58)]

## Average replicates
CTCF$CTCF_4h <- rowMeans(CTCF[,c(1,2)])
CTCF$CTCF_D0 <- rowMeans(CTCF[,c(3,4)])
CTCF <- CTCF[,c(6,7,5)]

# Place everything into a list
Counts <- list()
Counts$CEBPb <- CEBPb
Counts$KLF4 <- KLF4
Counts$RXR <- RXR
Counts$CEBPd <- CEBPd
Counts$Jun <- Jun 
Counts$CTCF <- CTCF

## Plot it all
par(mfcol=c(3,2))
for(i in 1:6) {
  boxplot(Counts[[i]][,2], Counts[[i]][,1], outline=F, col=c("darkblue", "darkred","grey"), ylab="ChIP-seq signal", names=c("D0", "4h"), xlab="Time point", main=names(Counts[i]), las=1)
}
```

[Back to start](../README.md)<br>
[Back to overview of Figure S4](../Links/FigureS4.md)