```R
## Process count data from sites with induced H3K27ac
# Import the data
Data <- read.delim("Data/ChIPseq/Counts/Histone_Marks_Induced_TFs_Cofactors.count")

# Take the relevant columns and set columns names
Data <- Data[,c(1,8,33:46)]
colnames(Data) <- c("PeakID","Input","CEBPb_D0","CEBPb_4h","CEBPd_D0","CEBPd_4h","KLF4_D0","KLF4_4h","RXR_D0","RXR_4h","Jun_D0","Jun_4h","CTCF_D0_E1","CTCF_D0_E2","CTCF_4h_E1","CTCF_4h_E2")

# Average replicates
Data$CTCF_D0 <- rowMeans(Data[,c("CTCF_D0_E1","CTCF_D0_E2")])
Data$CTCF_4h <- rowMeans(Data[,c("CTCF_4h_E1","CTCF_4h_E2")])

# Boxplot the results
par(mfcol=c(1,6))
boxplot(Data[ Data$CEBPb_D0 >= 20 | Data$CEBPb_4h >= 20 ,c("CEBPb_D0","CEBPb_4h")], outline=F, las=1, names=c("D0","4h"), main="CEBPb", ylab="Normalized tag count", col=c("blue3","red3"), ylim=c(0,150))
boxplot(Data[ Data$CEBPd_D0 >= 20 | Data$CEBPd_4h >= 20 ,c("CEBPd_D0","CEBPd_4h")], outline=F, las=1, names=c("D0","4h"), main="CEBPd", col=c("blue3","red3"), ylim=c(0,110))
boxplot(Data[ Data$KLF4_D0 >= 20 | Data$KLF4_4h >= 20 ,c("KLF4_D0","KLF4_4h")], outline=F, las=1, names=c("D0","4h"), main="KLF4", col=c("blue3","red3"), ylim=c(0,75))
boxplot(Data[ Data$RXR_D0 >= 20 | Data$RXR_4h >= 20 ,c("RXR_D0","RXR_4h")], outline=F, las=1, names=c("D0","4h"), main="RXR", col=c("blue3","red3"), ylim=c(0,90))
boxplot(Data[ Data$Jun_D0 >= 20 | Data$Jun_4h >= 20 ,c("Jun_D0","Jun_4h")], outline=F, las=1, names=c("D0","4h"), main="cJun", col=c("blue3","red3"), ylim=c(0,70))
boxplot(Data[ Data$CTCF_D0 >= 20 | Data$CTCF_4h >= 20 ,c("CTCF_D0","CTCF_4h")], outline=F, las=1, names=c("D0","4h"), main="CTCF", col=c("blue3","red3"), ylim=c(0,280))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S7](../Links/FigureS7.md)