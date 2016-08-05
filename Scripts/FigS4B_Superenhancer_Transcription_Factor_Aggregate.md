```R
## Import the data
Counts <- read.delim("Data/ChIPseq/Counts/Superenhancers_Cofactor_HistoneMarks_aggregate.count")
colnames(Counts)[1] <- "Distance"
Counts <- Counts[,c(1,grep("Coverage", colnames(Counts)))]

## Split out counts for each data set for averaging and place into a list
# TFs
CEBPb <- Counts[,c(1,grep("CEBPb_4h", colnames(Counts)),grep("CEBPb_D0", colnames(Counts)),52)]
CEBPd <- Counts[,c(1,grep("CEBPd_4h", colnames(Counts)),grep("CEBPd_D0", colnames(Counts)),52)]
Jun <- Counts[,c(1,grep("Jun_4h", colnames(Counts)),grep("Jun_D0", colnames(Counts)),52)]
KLF4 <- Counts[,c(1,grep("KLF4_4h", colnames(Counts)),grep("KLF4_D0", colnames(Counts)),52)]
RXR <- Counts[,c(1,grep("RXR_4h", colnames(Counts)),grep("RXR_D0", colnames(Counts)),52)]
CTCF <- Counts[,c(1,grep("CTCF_4h", colnames(Counts)),grep("CTCF_D0", colnames(Counts)),52)]

## Average replicates
CTCF$CTCF_4h <- rowMeans(CTCF[,c(2,3)])
CTCF$CTCF_D0 <- rowMeans(CTCF[,c(4,5)])
CTCF <- CTCF[,c(1,7,8,6)]

# Place everything into a list
Counts <- list()
Counts$CEBPb <- CEBPb
Counts$KLF4 <- KLF4
Counts$RXR <- RXR
Counts$CEBPd <- CEBPd
Counts$Jun <- Jun 
Counts$CTCF <- CTCF

## Plot the stuff
par(mfcol=c(3,2))
for (i in 1:6) {
plot(-10,-10,pch='',xlim=c(0,1), ylim=c(0,max(Counts[[i]][c(2,3)])), ylab="Normalized tag count", xlab="Fraction", main=names(Counts[i]), las=1, xaxt="n")
axis(1, at=seq(0,1,by=0.25), las=2)
lines(Counts[[i]][,1], Counts[[i]][,4], col="grey", lwd=2)
lines(Counts[[i]][,1], Counts[[i]][,2], col="darkblue", lwd=2)
lines(Counts[[i]][,1], Counts[[i]][,3], col="darkred", lwd=2)
}
```

[Back to start](../README.md)<br>
[Back to overview of Figure S4](../Links/FigureS4.md)