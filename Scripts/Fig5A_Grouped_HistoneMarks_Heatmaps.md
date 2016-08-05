```R
## Load the necessary packages
library(pheatmap)

# Import the data
Absent <- read.delim("Data/ChIPseq/Counts/Histone_Marks_Absent_Heatmap.count")
Induced <- read.delim("Data/ChIPseq/Counts/Histone_Marks_Induced_Heatmap.count")
Constitutive <- read.delim("Data/ChIPseq/Counts/Histone_Marks_Constitutive_Heatmap.count")
Repressed <- read.delim("Data/ChIPseq/Counts/Histone_Marks_Repressed_Heatmap.count")

# Average replicates
Starts <- c(403,1205,2007,2809,3611,4413)
Ends <- c(803,1605,2407,3209,4011,4813)
Insert <- 5215

# Loop through the different conditions and take mean in each bin
for (cond in 1:6) {
	Start <- Starts[cond]
	End <- Ends[cond]
	while (Start <= End) {
		Absent[,Insert] <- rowMeans(Absent[,c(Start,(Start+401))])
		Induced[,Insert] <- rowMeans(Induced[,c(Start,(Start+401))])
		Constitutive[,Insert] <- rowMeans(Constitutive[,c(Start,(Start+401))])
		Repressed[,Insert] <- rowMeans(Repressed[,c(Start,(Start+401))])
		Insert <- Insert + 1
		Start <- Start + 1
		}
	}

# Calculate the mean signal in all samples
Absent$Mean <- apply(Absent[,c(2:402,5215:7620)],1,FUN="mean")
Induced$Mean <- apply(Induced[,c(2:402,5215:7620)],1,FUN="mean")
Constitutive$Mean <- apply(Constitutive[,c(2:402,5215:7620)],1,FUN="mean")
Repressed$Mean <- apply(Repressed[,c(2:402,5215:7620)],1,FUN="mean")

# Order by the mean
Absent <- Absent[ order(-Absent$Mean),]
Induced <- Induced[ order(-Induced$Mean),]
Constitutive <- Constitutive[ order(-Constitutive$Mean),]
Repressed <- Repressed[ order(-Repressed$Mean),]

# Plot the heatmaps
pheatmap(Absent[,c(2:402,6418:6818,6017:6417,7220:7620,6819:7219,5616:6016,5215:5615)], cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white","blue","black"))(21), breaks = c(exp(seq(0,5,by=0.25))), show_rownames = F, show_colnames = F, border_color = NA, legend = F, filename="Absent.jpg")
pheatmap(Induced[,c(2:402,6418:6818,6017:6417,7220:7620,6819:7219,5616:6016,5215:5615)], cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white","blue","black"))(21), breaks = c(exp(seq(0,5,by=0.25))), show_rownames = F, show_colnames = F, border_color = NA, legend = F, filename="Induced.jpg")
pheatmap(Constitutive[,c(2:402,6418:6818,6017:6417,7220:7620,6819:7219,5616:6016,5215:5615)], cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white","blue","black"))(21), breaks = c(exp(seq(0,5,by=0.25))), show_rownames = F, show_colnames = F, border_color = NA, legend = F, filename="Constitutive.jpg")
pheatmap(Repressed[,c(2:402,6418:6818,6017:6417,7220:7620,6819:7219,5616:6016,5215:5615)], cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white","blue","black"))(21), breaks = c(exp(seq(0,5,by=0.25))), show_rownames = F, show_colnames = F, border_color = NA, legend = F, filename="Repressed.jpg")
```

[Back to start](../README.md)<br>
[Back to overview of Figure 5](../Links/Figure5.md)