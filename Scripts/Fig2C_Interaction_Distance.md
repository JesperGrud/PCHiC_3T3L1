```R
# Import the interactions
Interactions <- read.delim("Data/Interactions/PCHiC/Interactions.txt", header = T)

# Remove interactions from or to blacklisted fragments
Blacklist <- read.delim("Data/Interactions/PCHiC/Annotation/Blacklisted_Fragments.txt", header=F)
FilteredInteractions <- Interactions[ !(Interactions$oeID %in% Blacklist[,1]),]
FilteredInteractions <- FilteredInteractions[ !(FilteredInteractions$baitID %in% Blacklist[,1]),]

# Filter away interactions in trans
FilteredInteractions <- FilteredInteractions[ FilteredInteractions$oeChr == FilteredInteractions$baitChr,]

# Save the total number of interactions
Total <- nrow(FilteredInteractions)

# Filter away interactions spanning more than 1 megabase
FilteredInteractions <- FilteredInteractions[ abs(FilteredInteractions$dist) <= 1000000,]

# Get the absolute distance
Distances <- abs(FilteredInteractions$dist)

# Get the number of interactions in 1bp-sized bins
H <- hist(Distances, breaks = 1000000, plot=F)

# Extract mids (in kb) and calculate % of interactions
Mids <- H$mids/1000
Interactions <- cumsum(H$counts)/Total*100

# Make the plot
par(mfcol=c(1,1))
plot(x = Mids, y =, Interactions,las=1, 
	ylab="% interactions", xlab="Distance (kb)", 
	type = "l", col = "darkblue", lwd = 5, ylim=c(0,100))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 2](../Links/Figure2.md)