```R
# Import TADs
Zero_TAD <- read.delim("Data/Interactions/TADs/D0_40kb.domains.bed", header=F, skip=1)
Four_TAD <- read.delim("Data/Interactions/TADs/4h_40kb.domains.bed", header=F, skip=1)
Two_TAD <- read.delim("Data/Interactions/TADs/D2_40kb.domains.bed", header=F, skip=1)
ES_TAD <- read.delim("Data/Interactions/TADs/ESC_40kb.domains.bed", header=F, skip=1)

# Plot the number of TADs (NumberTADs.pdf)
barplot(c(nrow(Zero_TAD), nrow(Four_TAD), nrow(Two_TAD), nrow(ES_TAD)), las=1, ylab="Number of TADs",names=c("D0","4h","D2","ESC"), col=c("red3","green3","chocolate4","black"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 1](../Links/Figure1.md)
