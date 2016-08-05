```R
## Process count data from sites with induced H3K27ac
# Import the data
Data <- read.delim("Data/ChIPseq/Counts/Histone_Marks_Induced_TFs_Cofactors.count")

# Take the relevant columns and set columns names
Data <- Data[,c(1,8:32)]
colnames(Data) <- c("PeakID","Input","MED1_D0_E1","MED1_D0_E2","MED1_4h_E1","MED1_4h_E2","P300_D0_E1","P300_D0_E2","P300_4h_E1","P300_4h_E2","SMC1_D0_E1","SMC1_D0_E2","SMC1_4h_E1","SMC1_4h_E2","HDAC2_D0_E1","HDAC2_D0_E2","HDAC2_4h_E1","HDAC2_4h_E2","HDAC3_D0_E1","HDAC3_D0_E2","HDAC3_4h_E1","HDAC3_4h_E2","NCoR_D0_E1","NCoR_D0_E2","NCoR_4h_E1","NCoR_4h_E2")

# Average replicates
Data$MED1_D0 <- rowMeans(Data[,c("MED1_D0_E1","MED1_D0_E2")])
Data$MED1_4h <- rowMeans(Data[,c("MED1_4h_E1","MED1_4h_E2")])
Data$P300_D0 <- rowMeans(Data[,c("P300_D0_E1","P300_D0_E2")])
Data$P300_4h <- rowMeans(Data[,c("P300_4h_E1","P300_4h_E2")])
Data$SMC1_D0 <- rowMeans(Data[,c("SMC1_D0_E1","SMC1_D0_E2")])
Data$SMC1_4h <- rowMeans(Data[,c("SMC1_4h_E1","SMC1_4h_E2")])
Data$HDAC2_D0 <- rowMeans(Data[,c("HDAC2_D0_E1","HDAC2_D0_E2")])
Data$HDAC2_4h <- rowMeans(Data[,c("HDAC2_4h_E1","HDAC2_4h_E2")])
Data$HDAC3_D0 <- rowMeans(Data[,c("HDAC3_D0_E1","HDAC3_D0_E2")])
Data$HDAC3_4h <- rowMeans(Data[,c("HDAC3_4h_E1","HDAC3_4h_E2")])
Data$NCoR_D0 <- rowMeans(Data[,c("NCoR_D0_E1","NCoR_D0_E2")])
Data$NCoR_4h <- rowMeans(Data[,c("NCoR_4h_E1","NCoR_4h_E2")])

# Boxplot the results
par(mfcol=c(1,6))
boxplot(Data[ Data$MED1_D0 >= 20 | Data$MED1_4h >= 20 ,c("MED1_D0","MED1_4h")], outline=F, las=1, names=c("D0","4h"), main="MED1", ylab="Normalized tag count", col=c("blue3","red3"), ylim=c(0,160))
boxplot(Data[ Data$P300_D0 >= 20 | Data$P300_4h >= 20 ,c("P300_D0","P300_4h")], outline=F, las=1, names=c("D0","4h"), main="P300", col=c("blue3","red3"), ylim=c(0,50))
boxplot(Data[ Data$SMC1_D0 >= 20 | Data$SMC1_4h >= 20 ,c("SMC1_D0","SMC1_4h")], outline=F, las=1, names=c("D0","4h"), main="SMC1", col=c("blue3","red3"), ylim=c(0,110))
boxplot(Data[ Data$HDAC2_D0 >= 20 | Data$HDAC2_4h >= 20 ,c("HDAC2_D0","HDAC2_4h")], outline=F, las=1, names=c("D0","4h"), main="HDAC2", col=c("blue3","red3"), ylim=c(0,200))
boxplot(Data[ Data$HDAC3_D0 >= 20 | Data$HDAC3_4h >= 20 ,c("HDAC3_D0","HDAC3_4h")], outline=F, las=1, names=c("D0","4h"), main="HDAC3", col=c("blue3","red3"), ylim=c(0,40))
boxplot(Data[ Data$NCoR_D0 >= 20 | Data$NCoR_4h >= 20 ,c("NCoR_D0","NCoR_4h")], outline=F, las=1, names=c("D0","4h"), main="NCoR", col=c("blue3","red3"), ylim=c(0,60))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S7](../Links/FigureS7.md)