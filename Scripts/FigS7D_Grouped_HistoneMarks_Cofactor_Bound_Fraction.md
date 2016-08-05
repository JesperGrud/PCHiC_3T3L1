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

# Define time points and factors
Factors <- c("MED1","P300","SMC1","HDAC2","HDAC3","NCoR")
Time <- c("D0","4h")

# Make a matrix to capture the result
Result <- matrix(nrow=2, ncol=6)
colnames(Result) <- Factors
rownames(Result) <- Time

# Loop through factors and time points and calculate the bound fraction
for (factor in 1:6) {
	for (time in 1:2) {
		Result[time, factor]<- sum(Data[ ,colnames(Data) == paste(Factors[factor],Time[time], sep="_")] >= 20)/nrow(Data)
		}
	}
	
# Plot the results
par(mfcol=c(1,1))
barplot(Result, beside=T, ylim=c(0,0.6), las=2, col=c("blue3","red3"), ylab="Fraction bound by factor")
```

[Back to start](../README.md)<br>
[Back to overview of Figure S7](../Links/FigureS7.md)