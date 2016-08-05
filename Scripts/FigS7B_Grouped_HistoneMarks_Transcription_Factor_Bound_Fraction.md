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

# Define time points and factors
Factors <- c("CEBPb","CEBPd","KLF4","RXR","Jun","CTCF")
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
barplot(Result, beside=T, ylim=c(0,0.4), las=2, col=c("blue3","red3"), ylab="Fraction bound by factor")
```

[Back to start](../README.md)<br>
[Back to overview of Figure S7](../Links/FigureS7.md)