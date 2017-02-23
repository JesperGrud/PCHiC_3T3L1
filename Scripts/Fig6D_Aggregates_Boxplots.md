```R
## Import all data
# Aggregate plots
Aggregates <- read.delim("Data/ChIPseq/Counts/Superenhancers_Cofactor_HistoneMarks_aggregate.count")
colnames(Aggregates)[1] <- "Bin"
Aggregates <- Aggregates[,c(1,grep("Coverage",colnames(Aggregates)))]

# Boxplots
Boxplots <- read.delim("Data/ChIPseq/Counts/Superenhancers_Cofactor_HistoneMarks_window.count")
colnames(Boxplots)[1] <- "PeakID"

# Plot everything side-by-side
par(fig=c(0,0.2,0.66,1))
plot(Aggregates[,1], rowMeans(Aggregates[,grep("MED1_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="MED1", ylim=c(0,20), las=1, xlab="")
lines(Aggregates[,1], rowMeans(Aggregates[,grep("MED1_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.2,0.33,0.66,1), new=T)
boxplot(rowMeans(Boxplots[,grep("MED1_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("MED1_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="MED1", ylim=c(0,2500), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.33,0.53,0.66,1), new=T)
plot(Aggregates[,1], rowMeans(Aggregates[,grep("P300_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="P300", ylim=c(0,5), las=1, xlab="")
lines(Aggregates[,1], rowMeans(Aggregates[,grep("P300_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.53,0.66,0.66,1), new=T)
boxplot(rowMeans(Boxplots[,grep("P300_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("P300_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="P300", ylim=c(0,750), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.66,0.86,0.66,1), new=T)
plot(Aggregates[,1], rowMeans(Aggregates[,grep("SMC1_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="SMC1", ylim=c(0,15), las=1, xlab="")
lines(Aggregates[,1], rowMeans(Aggregates[,grep("SMC1_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.86,1,0.66,1), new=T)
boxplot(rowMeans(Boxplots[,grep("SMC1_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("SMC1_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="SMC1", ylim=c(0,1600), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0,0.2,0.33,0.66), new=T)
plot(Aggregates[,1], rowMeans(Aggregates[,grep("HDAC2_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="HDAC2", ylim=c(0,20), las=1, xlab="")
lines(Aggregates[,1], rowMeans(Aggregates[,grep("HDAC2_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.2,0.33,0.33,0.66), new=T)
boxplot(rowMeans(Boxplots[,grep("HDAC2_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("HDAC2_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="HDAC2", ylim=c(0,3300), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.33,0.53,0.33,0.66), new=T)
plot(Aggregates[,1], rowMeans(Aggregates[,grep("HDAC3_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="HDAC3", ylim=c(0,5), las=1, xlab="")
lines(Aggregates[,1], rowMeans(Aggregates[,grep("HDAC3_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.53,0.66,0.33,0.66), new=T)
boxplot(rowMeans(Boxplots[,grep("HDAC3_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("HDAC3_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="HDAC3", ylim=c(0,700), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.66,0.86,0.33,0.66), new=T)
plot(Aggregates[,1], rowMeans(Aggregates[,grep("NCoR_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="NCoR", ylim=c(0,5), las=1, xlab="")
lines(Aggregates[,1], rowMeans(Aggregates[,grep("NCoR_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.86,1,0.33,0.66), new=T)
boxplot(rowMeans(Boxplots[,grep("NCoR_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("NCoR_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="NCoR", ylim=c(0,850), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0,0.2,0,0.33), new=T)
plot(Aggregates[,1], rowMeans(Aggregates[,grep("H3K4me1_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="H3K4me1", ylim=c(0,5), las=1, xlab="")
lines(Aggregates[,1], rowMeans(Aggregates[,grep("H3K4me1_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.2,0.33,0,0.33), new=T)
boxplot(rowMeans(Boxplots[,grep("H3K4me1_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("H3K4me1_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="H3K4me1", ylim=c(0,1500), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.33,0.53,0,0.33), new=T)
plot(Aggregates[,1], rowMeans(Aggregates[,grep("H3K4me2_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="H3K4me2", ylim=c(0,5), las=1, xlab="")
lines(Aggregates[,1], rowMeans(Aggregates[,grep("H3K4me2_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.53,0.66,0,0.33), new=T)
boxplot(rowMeans(Boxplots[,grep("H3K4me2_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("H3K4me2_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="H3K4me2", ylim=c(0,1500), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.66,0.86,0,0.33), new=T)
plot(Aggregates[,1], rowMeans(Aggregates[,grep("H3K27ac_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="H3K27ac", ylim=c(0,10), las=1, xlab="")
lines(Aggregates[,1], rowMeans(Aggregates[,grep("H3K27ac_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.86,1,0,0.33), new=T)
boxplot(rowMeans(Boxplots[,grep("H3K27ac_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("H3K27ac_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="H3K27ac", ylim=c(0,4000), las=1, outline=F, names=c("D0","4h"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 6](../Links/Figure6.md)