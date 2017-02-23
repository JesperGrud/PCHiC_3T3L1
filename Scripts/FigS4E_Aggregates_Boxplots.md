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
par(fig=c(0,0.2,0.5,1))
plot(Aggregates[,1], Aggregates[,grep("CEBPb_4h",colnames(Aggregates))], type="l", col ="red3", ylab="Normalized tag count", main="CEBPb", ylim=c(0,10), las=1)
lines(Aggregates[,1], Aggregates[,grep("CEBPb_D0",colnames(Aggregates))], col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.2,0.33,0.5,1), new=T)
boxplot(Boxplots[,grep("CEBPb_D0",colnames(Boxplots))],Boxplots[,grep("CEBPb_4h",colnames(Boxplots))], type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="CEBPb", ylim=c(0,1500), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.33,0.53,0.5,1), new=T)
plot(Aggregates[,1], Aggregates[,grep("CEBPd_4h",colnames(Aggregates))], type="l", col ="red3", ylab="Normalized tag count", main="CEBPd", ylim=c(0,10), las=1)
lines(Aggregates[,1], Aggregates[,grep("CEBPd_D0",colnames(Aggregates))], col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.53,0.66,0.5,1), new=T)
boxplot(Boxplots[,grep("CEBPd_D0",colnames(Boxplots))],Boxplots[,grep("CEBPd_4h",colnames(Boxplots))], type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="CEBPd", ylim=c(0,1250), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.66,0.86,0.5,1), new=T)
plot(Aggregates[,1], Aggregates[,grep("RXR_4h",colnames(Aggregates))], type="l", col ="red3", ylab="Normalized tag count", main="RXR", ylim=c(0,10), las=1)
lines(Aggregates[,1], Aggregates[,grep("RXR_D0",colnames(Aggregates))], col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.86,1,0.5,1), new=T)
boxplot(Boxplots[,grep("RXR_D0",colnames(Boxplots))],Boxplots[,grep("RXR_4h",colnames(Boxplots))], type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="RXR", ylim=c(0,1200), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0,0.2,0,0.5), new=T)
plot(Aggregates[,1], Aggregates[,grep("KLF4_4h",colnames(Aggregates))], type="l", col ="red3", ylab="Normalized tag count", main="KLF4", ylim=c(0,10), las=1)
lines(Aggregates[,1], Aggregates[,grep("KLF4_D0",colnames(Aggregates))], col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.2,0.33,0,0.5), new=T)
boxplot(Boxplots[,grep("KLF4_D0",colnames(Boxplots))],Boxplots[,grep("KLF4_4h",colnames(Boxplots))], type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="KLF4", ylim=c(0,1000), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.33,0.53,0,0.5), new=T)
plot(Aggregates[,1], Aggregates[,grep("cJun_4h",colnames(Aggregates))], type="l", col ="red3", ylab="Normalized tag count", main="cJun", ylim=c(0,10), las=1)
lines(Aggregates[,1], Aggregates[,grep("cJun_D0",colnames(Aggregates))], col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.53,0.66,0,0.5), new=T)
boxplot(Boxplots[,grep("cJun_D0",colnames(Boxplots))],Boxplots[,grep("cJun_4h",colnames(Boxplots))], type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="cJun", ylim=c(0,750), las=1, outline=F, names=c("D0","4h"))

par(fig=c(0.66,0.86,0,0.5), new=T)
plot(Aggregates[,1], rowMeans(Aggregates[,grep("CTCF_4h",colnames(Aggregates))]), type="l", col ="red3", ylab="Normalized tag count", main="CTCF", ylim=c(0,10), las=1)
lines(Aggregates[,1], rowMeans(Aggregates[,grep("CTCF_D0",colnames(Aggregates))]), col="blue3")
lines(Aggregates[,1], Aggregates[,52], col ="grey")
par(fig=c(0.86,1,0,0.5), new=T)
boxplot(rowMeans(Boxplots[,grep("CTCF_D0",colnames(Boxplots))]),rowMeans(Boxplots[,grep("CTCF_4h",colnames(Boxplots))]), type="l", col = c("blue3","red3"), ylab="Normalized tag count", main="CTCF", ylim=c(0,1000), las=1, outline=F, names=c("D0","4h"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S4](../Links/FigureS4.md)
