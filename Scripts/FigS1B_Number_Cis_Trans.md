```R
## Import the interactions
Interactions <- read.delim("Data/Interactions/PCHiC/Interactions.txt", header = T)

## Remove interactions from or to blacklisted fragments
Blacklist <- read.delim("Data/Interactions/PCHiC/Annotation/Blacklisted_Fragments.txt", header=F)
FilteredInteractions <- Interactions[ !(Interactions$oeID %in% Blacklist[,1]),]
FilteredInteractions <- FilteredInteractions[ !(FilteredInteractions$baitID %in% Blacklist[,1]),]

## Calculate the percent of interactions in cis and in trans
Trans <- nrow(FilteredInteractions[ FilteredInteractions$oeChr != FilteredInteractions$baitChr,]) / nrow(FilteredInteractions)
Cis <- nrow(FilteredInteractions[ FilteredInteractions$oeChr == FilteredInteractions$baitChr,]) / nrow(FilteredInteractions)

## Plot a pie chart
par(mfcol=c(1,1))
pie(c(Cis,Trans), labels = c(paste("Cis"," - ",signif(Cis*100,3),"%",sep=""),
	paste("Trans"," - ",signif(Trans*100,3),"%",sep="")), col = c("blue","red"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)
