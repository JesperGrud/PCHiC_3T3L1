```R
## Load libraries
library(zoo)

## Process the count data 
# Import the count data
Counts <- read.delim("Data/Interactions/TADs/TAD_aggregate.count")

# Summarized plus and minus strand tags
Counts$CTCF_E1 <- rowSums(Counts[,c(3,4)])
Counts$CTCF_E2 <- rowSums(Counts[,c(6,7)])
Counts$MED1_E1 <- rowSums(Counts[,c(9,10)])
Counts$MED1_E2 <- rowSums(Counts[,c(12,13)])
Counts$SMC1_E1 <- rowSums(Counts[,c(15,16)])
Counts$SMC1_E2 <- rowSums(Counts[,c(18,19)])
Counts$Input <- rowSums(Counts[,c(21,22)])

# Average replicates
Counts$CTCF <- rowMeans(Counts[,c("CTCF_E1","CTCF_E2")])
Counts$MED1 <- rowMeans(Counts[,c("MED1_E1","MED1_E2")])
Counts$SMC1 <- rowMeans(Counts[,c("SMC1_E1","SMC1_E2")])

# Plot the rolling mean
plot(rollmean(Counts$CTCF, 50), type="l", lwd=2, xaxt="n", xlab="Fraction of TAD", ylab="Tag density", ylim=c(0.0038,0.008), col="deeppink", las=1)
lines(rollmean(Counts$MED1, 50), type="l", col="purple", lwd=2)
lines(rollmean(Counts$SMC1, 50), type="l", col="darkgreen", lwd=2)
lines(rollmean(Counts$Input, 50), type="l", col="grey", lwd=2)
axis(1, at=seq(0,10000,by=2000), labels=seq(0,1,by=0.2))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 1](../Links/Figure1.md)