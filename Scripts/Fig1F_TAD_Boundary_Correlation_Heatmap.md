```R
## Identify TADs and boundary regions
# Source the necessary packages
library(pheatmap)

# Import the data
Zero <- read.table("Data/Interactions/TADs/D0.DI", quote="\"", comment.char="", skip=1)

## Call TADs on sDI
# Convert sDI to z-scores
Mean <- mean(Zero$V4)
SD <- sd(Zero$V4)
Zero$V4 <- (Zero$V4 - Mean)/SD
rm(Mean)
rm(SD)

All <- data.frame()
Bounds <- data.frame()
Chromosomes <- c(paste("chr",seq(1,19,by=1),sep=""),"chrX")
for (Chr in Chromosomes) {
  # Subset to a single chromosome
  Subset <- Zero[ Zero$V1 == Chr,]
  
  # Identify transition points (between positive and negative bias)
  Subset$Transition <- 0
  for (i in 1:(nrow(Subset)-1)) {
    SignCurr <- sign(Subset[i,4])
    SignNext <- sign(Subset[(i+1),4])
    if (SignCurr != SignNext) { Subset[c(i),"Transition"] <- 1}
  }
  
  # Identify each region of bias
  Region <- 1
  Subset$Region <- 0
  for (i in 1:nrow(Subset)) {
    Subset[i,"Region"] <- Region
    if (Subset[i,"Transition"] == 1) { Region <- Region + 1}
  }
  
  # Calculate the length of each region
  Subset$Count <- 1
  Counts <- aggregate(Subset$Count, by=list(Subset$Region), FUN="sum")
  colnames(Counts) <- c("Region","Length")
  Subset <- merge(Subset, Counts, by="Region")
  
  # Identify regions longer or equal to 10 bins
  Subset$Long <- 0
  Subset[ Subset$Length >= 10,"Long"] <- 1
  
  # Mark each region for sufficient bias and set the sign
  Subset$Bias <- 0
  Subset$Sign <- 0
  for (i in unique(Subset$Region)) {
    Tmp <- Subset[ Subset$Region == i,]
    if (max(abs(Tmp$V4)) >= 0.5) {
      Subset[ Subset$Region == i, "Bias"] <- 1
      Subset[ Subset$Region == i, "Sign"] <- sign(max(Tmp$V4))
    }
  }
  
  # Remove regions that are too short (< 10 bins = 50kb)
  Subset[ Subset$Length <= 10, "Bias"] <- 0
  
  # Identify biased regions
  PositiveRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == 1, "Region"])
  NegativeRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == -1, "Region"])
  
  ## Identify type 1 boundaries (left > right)
  # For each negative region, get the closest positive region
  T1 <- data.frame()
  for (i in 1:length(NegativeRegions)) {
    if (length(PositiveRegions[ PositiveRegions >= NegativeRegions[i]]) > 0) {
      T1[i,1] <- NegativeRegions[i]
      T1[i,2] <- min(PositiveRegions[ PositiveRegions >= NegativeRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T1 <- T1[ order(T1[,2], -T1[,1]),]
  T1 <- T1[ duplicated(T1[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T1$Keep <- 0
  T1[ T1[,2] == (T1[,1] + 1),"Keep"] <- 1
  Proceeding <- T1[ T1$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT1 <- T1[ T1$Keep == 0,]
  for (i in 1:nrow(GappedT1)) {
    Test <- Subset[ Subset$Region < GappedT1[i,2] & Subset$Region > GappedT1[i,1],]
    if (nrow(Test) <= 8) { GappedT1[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT1[i,"Keep"] <- 0 }
    if (GappedT1[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT1[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT1[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT1[ GappedT1$Keep == 1,]) > 0) {
  GappedT1Keep <- GappedT1[ GappedT1$Keep == 1, c(4,5,3)]
  colnames(GappedT1Keep) <- colnames(Proceeding)
  T1 <- rbind(GappedT1Keep, Proceeding)
  }
  T1 <- T1[ order(T1[,1]),]
  
  # Identify type 2 boundaries (right > left)
  T2 <- data.frame()
  for (i in 1:length(PositiveRegions)) {
    if (length(NegativeRegions[ NegativeRegions >= PositiveRegions[i]]) > 0) {
      T2[i,1] <- PositiveRegions[i]
      T2[i,2] <- min(NegativeRegions[ NegativeRegions >= PositiveRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T2 <- T2[ order(T2[,2], -T2[,1]),]
  T2 <- T2[ duplicated(T2[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T2$Keep <- 0
  T2[ T2[,2] == (T2[,1] + 1),"Keep"] <- 1
  Proceeding <- T2[ T2$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT2 <- T2[ T2$Keep == 0,]
  for (i in 1:nrow(GappedT2)) {
    Test <- Subset[ Subset$Region < GappedT2[i,2] & Subset$Region > GappedT2[i,1],]
    if (nrow(Test) <= 8) { GappedT2[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT2[i,"Keep"] <- 0 }
    if (GappedT2[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT2[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT2[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT2[ GappedT2$Keep == 1,]) > 0) {
    GappedT2Keep <- GappedT2[ GappedT2$Keep == 1, c(4,5,3)]
    colnames(GappedT2Keep) <- colnames(Proceeding)
    T2 <- rbind(GappedT2Keep, Proceeding)
  }
  T2 <- T2[ order(T2[,1]),]
  
  # Get non-boundary biased regions
  PositiveNonbias <- data.frame(Region = PositiveRegions[ !(PositiveRegions %in% T1[,2]) & !(PositiveRegions %in% GappedT1[ GappedT1$Keep == 1,2])])
  NegativeNonbias <- data.frame(Region = NegativeRegions[ !(NegativeRegions %in% T2[,2]) & !(NegativeRegions %in% GappedT2[ GappedT2$Keep == 1,2])])
  
  # Stich together T1 regions
  TADs <- data.frame()
  for (i in 1:(nrow(T1)-1)) {
    TADs[i,1] <- T1[i,2]
    TADs[i,2] <- T1[(i+1),2]
  }
  
  # Remove TADs that doesnt cross a T2 boundary or multiple
  TADs$Keep <- 0
  TADs$T2 <- 0
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,2]
    Tmp <- data.frame(T2 = T2[ T2[,2] >= Min & T2[,2] <= Max,])
    if (nrow(Tmp) == 1) { 
      TADs[i,"Keep"] <- 1
      TADs[i,"T2"] <- Tmp[1,1]
    }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have negative bias between first T1 and T2 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,4]
    Tmp <- data.frame(NegativeNonbias = NegativeNonbias[ NegativeNonbias > Min & NegativeNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have positive bias between T2 and second T1 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,4]
    Max <- TADs[i,2]
    Tmp <- data.frame(PositiveNonbias = PositiveNonbias[ PositiveNonbias > Min & PositiveNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate TAD boundaries
  TADs$chr <- as.character(Chr)
  for (i in 1:nrow(TADs)) {
    Tmp <- Subset[ Subset$Region == TADs[i,1],]
    TADs[i,"start"] <- Tmp[ which.max(Tmp[,5]),3]
    Tmp <- Subset[ Subset$Region == (TADs[i,2]-1),]
    TADs[i,"end"] <- Tmp[ which.min(Tmp[,5]),3]+50000
  }
  
  # Filter TADs in low coverage areas (+/- 20 bins)
  Subset$Forward <- Subset[2,3] - Subset[1,3]
  for (nr in 2:nrow(Subset)) { Subset[nr,"Forward"] <- Subset[nr,3] - Subset[(nr-1),3] }
  
  for (m in 1:nrow(TADs)) {
    Start <- as.numeric(TADs[m,6])-50000
    End <- as.numeric(TADs[m,7])+50000
    Tmp <- Subset
    Tmp <- Tmp[ Tmp[,3] >= Start & Tmp[,3] < End,]
    if (max(Tmp$Forward) > 5000) { TADs[m,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate T1 boundaries
  T1$chr <- as.character(Chr)
  for (i in 1:nrow(T1)) {
    Tmp <- Subset[ Subset$Region == T1[i,2],]
    T1[i,"boundary"] <- min(Tmp[,3])
  }
  
  # Finalize TADs
  TADs <- TADs[,c(5,6,7)]
  
  # Combine all TADs
  All <- rbind(All, TADs)
  
  # Save T1 Boundaries
  T1 <- T1[,c(4,5)]
  Bounds <- rbind(Bounds, T1)
}

# Rename the TAD calls
Zero_TADs <- All
Zero_Bounds <- Bounds

## Call TADs on 4 hour data
# Import the data
Zero <- read.table("Data/Interactions/TADs/H4.DI", quote="\"", comment.char="", skip=1)

## Call TADs on sDI
# Convert sDI to z-scores
Mean <- mean(Zero$V4)
SD <- sd(Zero$V4)
Zero$V4 <- (Zero$V4 - Mean)/SD
rm(Mean)
rm(SD)

Bounds <- data.frame()
All <- data.frame()
Chromosomes <- c(paste("chr",seq(1,19,by=1),sep=""),"chrX")
for (Chr in Chromosomes) {
  # Subset to a single chromosome
  Subset <- Zero[ Zero$V1 == Chr,]
  
  # Identify transition points (between positive and negative bias)
  Subset$Transition <- 0
  for (i in 1:(nrow(Subset)-1)) {
    SignCurr <- sign(Subset[i,4])
    SignNext <- sign(Subset[(i+1),4])
    if (SignCurr != SignNext) { Subset[c(i),"Transition"] <- 1}
  }
  
  # Identify each region of bias
  Region <- 1
  Subset$Region <- 0
  for (i in 1:nrow(Subset)) {
    Subset[i,"Region"] <- Region
    if (Subset[i,"Transition"] == 1) { Region <- Region + 1}
  }
  
  # Calculate the length of each region
  Subset$Count <- 1
  Counts <- aggregate(Subset$Count, by=list(Subset$Region), FUN="sum")
  colnames(Counts) <- c("Region","Length")
  Subset <- merge(Subset, Counts, by="Region")
  
  # Identify regions longer or equal to 10 bins
  Subset$Long <- 0
  Subset[ Subset$Length >= 10,"Long"] <- 1
  
  # Mark each region for sufficient bias and set the sign
  Subset$Bias <- 0
  Subset$Sign <- 0
  for (i in unique(Subset$Region)) {
    Tmp <- Subset[ Subset$Region == i,]
    if (max(abs(Tmp$V4)) >= 0.5) {
      Subset[ Subset$Region == i, "Bias"] <- 1
      Subset[ Subset$Region == i, "Sign"] <- sign(max(Tmp$V4))
    }
  }
  
  # Remove regions that are too short (< 10 bins = 50kb)
  Subset[ Subset$Length <= 10, "Bias"] <- 0
  
  # Identify biased regions
  PositiveRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == 1, "Region"])
  NegativeRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == -1, "Region"])
  
  ## Identify type 1 boundaries (left > right)
  # For each negative region, get the closest positive region
  T1 <- data.frame()
  for (i in 1:length(NegativeRegions)) {
    if (length(PositiveRegions[ PositiveRegions >= NegativeRegions[i]]) > 0) {
      T1[i,1] <- NegativeRegions[i]
      T1[i,2] <- min(PositiveRegions[ PositiveRegions >= NegativeRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T1 <- T1[ order(T1[,2], -T1[,1]),]
  T1 <- T1[ duplicated(T1[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T1$Keep <- 0
  T1[ T1[,2] == (T1[,1] + 1),"Keep"] <- 1
  Proceeding <- T1[ T1$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT1 <- T1[ T1$Keep == 0,]
  for (i in 1:nrow(GappedT1)) {
    Test <- Subset[ Subset$Region < GappedT1[i,2] & Subset$Region > GappedT1[i,1],]
    if (nrow(Test) <= 8) { GappedT1[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT1[i,"Keep"] <- 0 }
    if (GappedT1[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT1[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT1[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT1[ GappedT1$Keep == 1,]) > 0) {
    GappedT1Keep <- GappedT1[ GappedT1$Keep == 1, c(4,5,3)]
    colnames(GappedT1Keep) <- colnames(Proceeding)
    T1 <- rbind(GappedT1Keep, Proceeding)
  }
  T1 <- T1[ order(T1[,1]),]
  
  # Identify type 2 boundaries (right > left)
  T2 <- data.frame()
  for (i in 1:length(PositiveRegions)) {
    if (length(NegativeRegions[ NegativeRegions >= PositiveRegions[i]]) > 0) {
      T2[i,1] <- PositiveRegions[i]
      T2[i,2] <- min(NegativeRegions[ NegativeRegions >= PositiveRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T2 <- T2[ order(T2[,2], -T2[,1]),]
  T2 <- T2[ duplicated(T2[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T2$Keep <- 0
  T2[ T2[,2] == (T2[,1] + 1),"Keep"] <- 1
  Proceeding <- T2[ T2$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT2 <- T2[ T2$Keep == 0,]
  for (i in 1:nrow(GappedT2)) {
    Test <- Subset[ Subset$Region < GappedT2[i,2] & Subset$Region > GappedT2[i,1],]
    if (nrow(Test) <= 8) { GappedT2[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT2[i,"Keep"] <- 0 }
    if (GappedT2[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT2[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT2[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT2[ GappedT2$Keep == 1,]) > 0) {
    GappedT2Keep <- GappedT2[ GappedT2$Keep == 1, c(4,5,3)]
    colnames(GappedT2Keep) <- colnames(Proceeding)
    T2 <- rbind(GappedT2Keep, Proceeding)
  }
  T2 <- T2[ order(T2[,1]),]
  
  # Get non-boundary biased regions
  PositiveNonbias <- data.frame(Region = PositiveRegions[ !(PositiveRegions %in% T1[,2]) & !(PositiveRegions %in% GappedT1[ GappedT1$Keep == 1,2])])
  NegativeNonbias <- data.frame(Region = NegativeRegions[ !(NegativeRegions %in% T2[,2]) & !(NegativeRegions %in% GappedT2[ GappedT2$Keep == 1,2])])
  
  # Stich together T1 regions
  TADs <- data.frame()
  for (i in 1:(nrow(T1)-1)) {
    TADs[i,1] <- T1[i,2]
    TADs[i,2] <- T1[(i+1),2]
  }
  
  # Remove TADs that doesnt cross a T2 boundary or multiple
  TADs$Keep <- 0
  TADs$T2 <- 0
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,2]
    Tmp <- data.frame(T2 = T2[ T2[,2] >= Min & T2[,2] <= Max,])
    if (nrow(Tmp) == 1) { 
      TADs[i,"Keep"] <- 1
      TADs[i,"T2"] <- Tmp[1,1]
    }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have negative bias between first T1 and T2 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,4]
    Tmp <- data.frame(NegativeNonbias = NegativeNonbias[ NegativeNonbias > Min & NegativeNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have positive bias between T2 and second T1 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,4]
    Max <- TADs[i,2]
    Tmp <- data.frame(PositiveNonbias = PositiveNonbias[ PositiveNonbias > Min & PositiveNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate TAD boundaries
  TADs$chr <- as.character(Chr)
  for (i in 1:nrow(TADs)) {
    Tmp <- Subset[ Subset$Region == TADs[i,1],]
    TADs[i,"start"] <- Tmp[ which.max(Tmp[,5]),3]
    Tmp <- Subset[ Subset$Region == (TADs[i,2]-1),]
    TADs[i,"end"] <- Tmp[ which.min(Tmp[,5]),3]+50000
  }
  
  # Filter TADs in low coverage areas (+/- 20 bins)
  Subset$Forward <- Subset[2,3] - Subset[1,3]
  for (nr in 2:nrow(Subset)) { Subset[nr,"Forward"] <- Subset[nr,3] - Subset[(nr-1),3] }
  
  for (m in 1:nrow(TADs)) {
    Start <- as.numeric(TADs[m,6])-50000
    End <- as.numeric(TADs[m,7])+50000
    Tmp <- Subset
    Tmp <- Tmp[ Tmp[,3] >= Start & Tmp[,3] < End,]
    if (max(Tmp$Forward) > 5000) { TADs[m,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate T1 boundaries
  T1$chr <- as.character(Chr)
  for (i in 1:nrow(T1)) {
    Tmp <- Subset[ Subset$Region == T1[i,2],]
    T1[i,"boundary"] <- min(Tmp[,3])
  }
  
  # Finalize TADs
  TADs <- TADs[,c(5,6,7)]
  
  # Combine all TADs
  All <- rbind(All, TADs)
  
  # Save T1 Boundaries
  T1 <- T1[,c(4,5)]
  Bounds <- rbind(Bounds, T1)
}

# Rename the TAD calls
Four_TADs <- All
Four_Bounds <- Bounds

## Call TADs on D2 data
# Import the data
Zero <- read.table("Data/Interactions/TADs/D2.DI", quote="\"", comment.char="", skip=1)

## Call TADs on sDI
# Convert sDI to z-scores
Mean <- mean(Zero$V4)
SD <- sd(Zero$V4)
Zero$V4 <- (Zero$V4 - Mean)/SD
rm(Mean)
rm(SD)

Bounds <- data.frame()
All <- data.frame()
Chromosomes <- c(paste("chr",seq(1,19,by=1),sep=""),"chrX")
for (Chr in Chromosomes) {
  # Subset to a single chromosome
  Subset <- Zero[ Zero$V1 == Chr,]
  
  # Identify transition points (between positive and negative bias)
  Subset$Transition <- 0
  for (i in 1:(nrow(Subset)-1)) {
    SignCurr <- sign(Subset[i,4])
    SignNext <- sign(Subset[(i+1),4])
    if (SignCurr != SignNext) { Subset[c(i),"Transition"] <- 1}
  }
  
  # Identify each region of bias
  Region <- 1
  Subset$Region <- 0
  for (i in 1:nrow(Subset)) {
    Subset[i,"Region"] <- Region
    if (Subset[i,"Transition"] == 1) { Region <- Region + 1}
  }
  
  # Calculate the length of each region
  Subset$Count <- 1
  Counts <- aggregate(Subset$Count, by=list(Subset$Region), FUN="sum")
  colnames(Counts) <- c("Region","Length")
  Subset <- merge(Subset, Counts, by="Region")
  
  # Identify regions longer or equal to 10 bins
  Subset$Long <- 0
  Subset[ Subset$Length >= 10,"Long"] <- 1
  
  # Mark each region for sufficient bias and set the sign
  Subset$Bias <- 0
  Subset$Sign <- 0
  for (i in unique(Subset$Region)) {
    Tmp <- Subset[ Subset$Region == i,]
    if (max(abs(Tmp$V4)) >= 0.5) {
      Subset[ Subset$Region == i, "Bias"] <- 1
      Subset[ Subset$Region == i, "Sign"] <- sign(max(Tmp$V4))
    }
  }
  
  # Remove regions that are too short (< 10 bins = 50kb)
  Subset[ Subset$Length <= 10, "Bias"] <- 0
  
  # Identify biased regions
  PositiveRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == 1, "Region"])
  NegativeRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == -1, "Region"])
  
  ## Identify type 1 boundaries (left > right)
  # For each negative region, get the closest positive region
  T1 <- data.frame()
  for (i in 1:length(NegativeRegions)) {
    if (length(PositiveRegions[ PositiveRegions >= NegativeRegions[i]]) > 0) {
      T1[i,1] <- NegativeRegions[i]
      T1[i,2] <- min(PositiveRegions[ PositiveRegions >= NegativeRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T1 <- T1[ order(T1[,2], -T1[,1]),]
  T1 <- T1[ duplicated(T1[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T1$Keep <- 0
  T1[ T1[,2] == (T1[,1] + 1),"Keep"] <- 1
  Proceeding <- T1[ T1$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT1 <- T1[ T1$Keep == 0,]
  for (i in 1:nrow(GappedT1)) {
    Test <- Subset[ Subset$Region < GappedT1[i,2] & Subset$Region > GappedT1[i,1],]
    if (nrow(Test) <= 8) { GappedT1[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT1[i,"Keep"] <- 0 }
    if (GappedT1[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT1[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT1[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT1[ GappedT1$Keep == 1,]) > 0) {
    GappedT1Keep <- GappedT1[ GappedT1$Keep == 1, c(4,5,3)]
    colnames(GappedT1Keep) <- colnames(Proceeding)
    T1 <- rbind(GappedT1Keep, Proceeding)
  }
  T1 <- T1[ order(T1[,1]),]
  
  # Identify type 2 boundaries (right > left)
  T2 <- data.frame()
  for (i in 1:length(PositiveRegions)) {
    if (length(NegativeRegions[ NegativeRegions >= PositiveRegions[i]]) > 0) {
      T2[i,1] <- PositiveRegions[i]
      T2[i,2] <- min(NegativeRegions[ NegativeRegions >= PositiveRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T2 <- T2[ order(T2[,2], -T2[,1]),]
  T2 <- T2[ duplicated(T2[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T2$Keep <- 0
  T2[ T2[,2] == (T2[,1] + 1),"Keep"] <- 1
  Proceeding <- T2[ T2$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT2 <- T2[ T2$Keep == 0,]
  for (i in 1:nrow(GappedT2)) {
    Test <- Subset[ Subset$Region < GappedT2[i,2] & Subset$Region > GappedT2[i,1],]
    if (nrow(Test) <= 8) { GappedT2[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT2[i,"Keep"] <- 0 }
    if (GappedT2[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT2[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT2[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT2[ GappedT2$Keep == 1,]) > 0) {
    GappedT2Keep <- GappedT2[ GappedT2$Keep == 1, c(4,5,3)]
    colnames(GappedT2Keep) <- colnames(Proceeding)
    T2 <- rbind(GappedT2Keep, Proceeding)
  }
  T2 <- T2[ order(T2[,1]),]
  
  # Get non-boundary biased regions
  PositiveNonbias <- data.frame(Region = PositiveRegions[ !(PositiveRegions %in% T1[,2]) & !(PositiveRegions %in% GappedT1[ GappedT1$Keep == 1,2])])
  NegativeNonbias <- data.frame(Region = NegativeRegions[ !(NegativeRegions %in% T2[,2]) & !(NegativeRegions %in% GappedT2[ GappedT2$Keep == 1,2])])
  
  # Stich together T1 regions
  TADs <- data.frame()
  for (i in 1:(nrow(T1)-1)) {
    TADs[i,1] <- T1[i,2]
    TADs[i,2] <- T1[(i+1),2]
  }
  
  # Remove TADs that doesnt cross a T2 boundary or multiple
  TADs$Keep <- 0
  TADs$T2 <- 0
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,2]
    Tmp <- data.frame(T2 = T2[ T2[,2] >= Min & T2[,2] <= Max,])
    if (nrow(Tmp) == 1) { 
      TADs[i,"Keep"] <- 1
      TADs[i,"T2"] <- Tmp[1,1]
    }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have negative bias between first T1 and T2 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,4]
    Tmp <- data.frame(NegativeNonbias = NegativeNonbias[ NegativeNonbias > Min & NegativeNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have positive bias between T2 and second T1 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,4]
    Max <- TADs[i,2]
    Tmp <- data.frame(PositiveNonbias = PositiveNonbias[ PositiveNonbias > Min & PositiveNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate TAD boundaries
  TADs$chr <- as.character(Chr)
  for (i in 1:nrow(TADs)) {
    Tmp <- Subset[ Subset$Region == TADs[i,1],]
    TADs[i,"start"] <- Tmp[ which.max(Tmp[,5]),3]
    Tmp <- Subset[ Subset$Region == (TADs[i,2]-1),]
    TADs[i,"end"] <- Tmp[ which.min(Tmp[,5]),3]+50000
  }
  
  # Filter TADs in low coverage areas (+/- 20 bins)
  Subset$Forward <- Subset[2,3] - Subset[1,3]
  for (nr in 2:nrow(Subset)) { Subset[nr,"Forward"] <- Subset[nr,3] - Subset[(nr-1),3] }
  
  for (m in 1:nrow(TADs)) {
    Start <- as.numeric(TADs[m,6])-50000
    End <- as.numeric(TADs[m,7])+50000
    Tmp <- Subset
    Tmp <- Tmp[ Tmp[,3] >= Start & Tmp[,3] < End,]
    if (max(Tmp$Forward) > 5000) { TADs[m,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate T1 boundaries
  T1$chr <- as.character(Chr)
  for (i in 1:nrow(T1)) {
    Tmp <- Subset[ Subset$Region == T1[i,2],]
    T1[i,"boundary"] <- min(Tmp[,3])
  }
  
  # Finalize TADs
  TADs <- TADs[,c(5,6,7)]
  
  # Combine all TADs
  All <- rbind(All, TADs)
  
  # Save T1 Boundaries
  T1 <- T1[,c(4,5)]
  Bounds <- rbind(Bounds, T1)
}

# Rename the TAD calls
Two_TADs <- All
Two_Bounds <- Bounds

## Call TADs on D2 data
# Import the data
Zero <- read.table("Data/Interactions/TADs/ESC.DI", quote="\"", comment.char="", skip=1)

## Call TADs on sDI
# Convert sDI to z-scores
Mean <- mean(Zero$V4)
SD <- sd(Zero$V4)
Zero$V4 <- (Zero$V4 - Mean)/SD
rm(Mean)
rm(SD)

Bounds <- data.frame()
All <- data.frame()
Chromosomes <- c(paste("chr",seq(1,19,by=1),sep=""),"chrX")
for (Chr in Chromosomes) {
  # Subset to a single chromosome
  Subset <- Zero[ Zero$V1 == Chr,]
  
  # Identify transition points (between positive and negative bias)
  Subset$Transition <- 0
  for (i in 1:(nrow(Subset)-1)) {
    SignCurr <- sign(Subset[i,4])
    SignNext <- sign(Subset[(i+1),4])
    if (SignCurr != SignNext) { Subset[c(i),"Transition"] <- 1}
  }
  
  # Identify each region of bias
  Region <- 1
  Subset$Region <- 0
  for (i in 1:nrow(Subset)) {
    Subset[i,"Region"] <- Region
    if (Subset[i,"Transition"] == 1) { Region <- Region + 1}
  }
  
  # Calculate the length of each region
  Subset$Count <- 1
  Counts <- aggregate(Subset$Count, by=list(Subset$Region), FUN="sum")
  colnames(Counts) <- c("Region","Length")
  Subset <- merge(Subset, Counts, by="Region")
  
  # Identify regions longer or equal to 10 bins
  Subset$Long <- 0
  Subset[ Subset$Length >= 10,"Long"] <- 1
  
  # Mark each region for sufficient bias and set the sign
  Subset$Bias <- 0
  Subset$Sign <- 0
  for (i in unique(Subset$Region)) {
    Tmp <- Subset[ Subset$Region == i,]
    if (max(abs(Tmp$V4)) >= 0.5) {
      Subset[ Subset$Region == i, "Bias"] <- 1
      Subset[ Subset$Region == i, "Sign"] <- sign(max(Tmp$V4))
    }
  }
  
  # Remove regions that are too short (< 10 bins = 50kb)
  Subset[ Subset$Length <= 10, "Bias"] <- 0
  
  # Identify biased regions
  PositiveRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == 1, "Region"])
  NegativeRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == -1, "Region"])
  
  ## Identify type 1 boundaries (left > right)
  # For each negative region, get the closest positive region
  T1 <- data.frame()
  for (i in 1:length(NegativeRegions)) {
    if (length(PositiveRegions[ PositiveRegions >= NegativeRegions[i]]) > 0) {
      T1[i,1] <- NegativeRegions[i]
      T1[i,2] <- min(PositiveRegions[ PositiveRegions >= NegativeRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T1 <- T1[ order(T1[,2], -T1[,1]),]
  T1 <- T1[ duplicated(T1[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T1$Keep <- 0
  T1[ T1[,2] == (T1[,1] + 1),"Keep"] <- 1
  Proceeding <- T1[ T1$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT1 <- T1[ T1$Keep == 0,]
  for (i in 1:nrow(GappedT1)) {
    Test <- Subset[ Subset$Region < GappedT1[i,2] & Subset$Region > GappedT1[i,1],]
    if (nrow(Test) <= 8) { GappedT1[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT1[i,"Keep"] <- 0 }
    if (GappedT1[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT1[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT1[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT1[ GappedT1$Keep == 1,]) > 0) {
    GappedT1Keep <- GappedT1[ GappedT1$Keep == 1, c(4,5,3)]
    colnames(GappedT1Keep) <- colnames(Proceeding)
    T1 <- rbind(GappedT1Keep, Proceeding)
  }
  T1 <- T1[ order(T1[,1]),]
  
  # Identify type 2 boundaries (right > left)
  T2 <- data.frame()
  for (i in 1:length(PositiveRegions)) {
    if (length(NegativeRegions[ NegativeRegions >= PositiveRegions[i]]) > 0) {
      T2[i,1] <- PositiveRegions[i]
      T2[i,2] <- min(NegativeRegions[ NegativeRegions >= PositiveRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T2 <- T2[ order(T2[,2], -T2[,1]),]
  T2 <- T2[ duplicated(T2[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T2$Keep <- 0
  T2[ T2[,2] == (T2[,1] + 1),"Keep"] <- 1
  Proceeding <- T2[ T2$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT2 <- T2[ T2$Keep == 0,]
  for (i in 1:nrow(GappedT2)) {
    Test <- Subset[ Subset$Region < GappedT2[i,2] & Subset$Region > GappedT2[i,1],]
    if (nrow(Test) <= 8) { GappedT2[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT2[i,"Keep"] <- 0 }
    if (GappedT2[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT2[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT2[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT2[ GappedT2$Keep == 1,]) > 0) {
    GappedT2Keep <- GappedT2[ GappedT2$Keep == 1, c(4,5,3)]
    colnames(GappedT2Keep) <- colnames(Proceeding)
    T2 <- rbind(GappedT2Keep, Proceeding)
  }
  T2 <- T2[ order(T2[,1]),]
  
  # Get non-boundary biased regions
  PositiveNonbias <- data.frame(Region = PositiveRegions[ !(PositiveRegions %in% T1[,2]) & !(PositiveRegions %in% GappedT1[ GappedT1$Keep == 1,2])])
  NegativeNonbias <- data.frame(Region = NegativeRegions[ !(NegativeRegions %in% T2[,2]) & !(NegativeRegions %in% GappedT2[ GappedT2$Keep == 1,2])])
  
  # Stich together T1 regions
  TADs <- data.frame()
  for (i in 1:(nrow(T1)-1)) {
    TADs[i,1] <- T1[i,2]
    TADs[i,2] <- T1[(i+1),2]
  }
  
  # Remove TADs that doesnt cross a T2 boundary or multiple
  TADs$Keep <- 0
  TADs$T2 <- 0
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,2]
    Tmp <- data.frame(T2 = T2[ T2[,2] >= Min & T2[,2] <= Max,])
    if (nrow(Tmp) == 1) { 
      TADs[i,"Keep"] <- 1
      TADs[i,"T2"] <- Tmp[1,1]
    }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have negative bias between first T1 and T2 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,4]
    Tmp <- data.frame(NegativeNonbias = NegativeNonbias[ NegativeNonbias > Min & NegativeNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have positive bias between T2 and second T1 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,4]
    Max <- TADs[i,2]
    Tmp <- data.frame(PositiveNonbias = PositiveNonbias[ PositiveNonbias > Min & PositiveNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate TAD boundaries
  TADs$chr <- as.character(Chr)
  for (i in 1:nrow(TADs)) {
    Tmp <- Subset[ Subset$Region == TADs[i,1],]
    TADs[i,"start"] <- Tmp[ which.max(Tmp[,5]),3]
    Tmp <- Subset[ Subset$Region == (TADs[i,2]-1),]
    TADs[i,"end"] <- Tmp[ which.min(Tmp[,5]),3]+50000
  }
  
  # Filter TADs in low coverage areas (+/- 20 bins)
  Subset$Forward <- Subset[2,3] - Subset[1,3]
  for (nr in 2:nrow(Subset)) { Subset[nr,"Forward"] <- Subset[nr,3] - Subset[(nr-1),3] }
  
  for (m in 1:nrow(TADs)) {
    Start <- as.numeric(TADs[m,6])-50000
    End <- as.numeric(TADs[m,7])+50000
    Tmp <- Subset
    Tmp <- Tmp[ Tmp[,3] >= Start & Tmp[,3] < End,]
    if (max(Tmp$Forward) > 5000) { TADs[m,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate T1 boundaries
  T1$chr <- as.character(Chr)
  for (i in 1:nrow(T1)) {
    Tmp <- Subset[ Subset$Region == T1[i,2],]
    T1[i,"boundary"] <- min(Tmp[,3])
  }
  
  # Finalize TADs
  TADs <- TADs[,c(5,6,7)]
  
  # Combine all TADs
  All <- rbind(All, TADs)
  
  # Save T1 Boundaries
  T1 <- T1[,c(4,5)]
  Bounds <- rbind(Bounds, T1)
}

# Rename the TAD calls
ESC_TADs <- All
ESC_Bounds <- Bounds

## Call TADs on D2 data
# Import the data
Zero <- read.table("Data/Interactions/TADs/Cortex.DI", quote="\"", comment.char="", skip=1)

## Call TADs on sDI
# Convert sDI to z-scores
Mean <- mean(Zero$V4)
SD <- sd(Zero$V4)
Zero$V4 <- (Zero$V4 - Mean)/SD
rm(Mean)
rm(SD)

All <- data.frame()
Bounds <- data.frame()
Chromosomes <- c(paste("chr",seq(1,19,by=1),sep=""),"chrX")
for (Chr in Chromosomes) {
  # Subset to a single chromosome
  Subset <- Zero[ Zero$V1 == Chr,]
  
  # Identify transition points (between positive and negative bias)
  Subset$Transition <- 0
  for (i in 1:(nrow(Subset)-1)) {
    SignCurr <- sign(Subset[i,4])
    SignNext <- sign(Subset[(i+1),4])
    if (SignCurr != SignNext) { Subset[c(i),"Transition"] <- 1}
  }
  
  # Identify each region of bias
  Region <- 1
  Subset$Region <- 0
  for (i in 1:nrow(Subset)) {
    Subset[i,"Region"] <- Region
    if (Subset[i,"Transition"] == 1) { Region <- Region + 1}
  }
  
  # Calculate the length of each region
  Subset$Count <- 1
  Counts <- aggregate(Subset$Count, by=list(Subset$Region), FUN="sum")
  colnames(Counts) <- c("Region","Length")
  Subset <- merge(Subset, Counts, by="Region")
  
  # Identify regions longer or equal to 10 bins
  Subset$Long <- 0
  Subset[ Subset$Length >= 10,"Long"] <- 1
  
  # Mark each region for sufficient bias and set the sign
  Subset$Bias <- 0
  Subset$Sign <- 0
  for (i in unique(Subset$Region)) {
    Tmp <- Subset[ Subset$Region == i,]
    if (max(abs(Tmp$V4)) >= 0.5) {
      Subset[ Subset$Region == i, "Bias"] <- 1
      Subset[ Subset$Region == i, "Sign"] <- sign(max(Tmp$V4))
    }
  }
  
  # Remove regions that are too short (< 10 bins = 50kb)
  Subset[ Subset$Length <= 10, "Bias"] <- 0
  
  # Identify biased regions
  PositiveRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == 1, "Region"])
  NegativeRegions <- unique(Subset[ Subset$Bias == 1 & Subset$Sign == -1, "Region"])
  
  ## Identify type 1 boundaries (left > right)
  # For each negative region, get the closest positive region
  T1 <- data.frame()
  for (i in 1:length(NegativeRegions)) {
    if (length(PositiveRegions[ PositiveRegions >= NegativeRegions[i]]) > 0) {
      T1[i,1] <- NegativeRegions[i]
      T1[i,2] <- min(PositiveRegions[ PositiveRegions >= NegativeRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T1 <- T1[ order(T1[,2], -T1[,1]),]
  T1 <- T1[ duplicated(T1[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T1$Keep <- 0
  T1[ T1[,2] == (T1[,1] + 1),"Keep"] <- 1
  Proceeding <- T1[ T1$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT1 <- T1[ T1$Keep == 0,]
  for (i in 1:nrow(GappedT1)) {
    Test <- Subset[ Subset$Region < GappedT1[i,2] & Subset$Region > GappedT1[i,1],]
    if (nrow(Test) <= 8) { GappedT1[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT1[i,"Keep"] <- 0 }
    if (GappedT1[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT1[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT1[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT1[ GappedT1$Keep == 1,]) > 0) {
    GappedT1Keep <- GappedT1[ GappedT1$Keep == 1, c(4,5,3)]
    colnames(GappedT1Keep) <- colnames(Proceeding)
    T1 <- rbind(GappedT1Keep, Proceeding)
  }
  T1 <- T1[ order(T1[,1]),]
  
  # Identify type 2 boundaries (right > left)
  T2 <- data.frame()
  for (i in 1:length(PositiveRegions)) {
    if (length(NegativeRegions[ NegativeRegions >= PositiveRegions[i]]) > 0) {
      T2[i,1] <- PositiveRegions[i]
      T2[i,2] <- min(NegativeRegions[ NegativeRegions >= PositiveRegions[i]])  
    }
  }
  
  # Keep the negative that is closest to the positive border in case there are multiple
  T2 <- T2[ order(T2[,2], -T2[,1]),]
  T2 <- T2[ duplicated(T2[,2])==F,]
  
  # Keep the ones that are directly preceeding each other
  T2$Keep <- 0
  T2[ T2[,2] == (T2[,1] + 1),"Keep"] <- 1
  Proceeding <- T2[ T2$Keep == 1,]
  
  # Analyze the gap for gapped boundaries
  GappedT2 <- T2[ T2$Keep == 0,]
  for (i in 1:nrow(GappedT2)) {
    Test <- Subset[ Subset$Region < GappedT2[i,2] & Subset$Region > GappedT2[i,1],]
    if (nrow(Test) <= 8) { GappedT2[i,"Keep"] <- 1 }
    if (max(abs(Test[,5])) >= 0.25) { GappedT2[i,"Keep"] <- 0 }
    if (GappedT2[i,"Keep"] == 1) {
      Test <- Test[ order(-Test$V2),]
      GappedT2[i,4] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"] - 1
      GappedT2[i,5] <- Test[ which.min(abs(Test[,3] - mean(Test[,3]))),"Region"]
    }
  }
  if (nrow(GappedT2[ GappedT2$Keep == 1,]) > 0) {
    GappedT2Keep <- GappedT2[ GappedT2$Keep == 1, c(4,5,3)]
    colnames(GappedT2Keep) <- colnames(Proceeding)
    T2 <- rbind(GappedT2Keep, Proceeding)
  }
  T2 <- T2[ order(T2[,1]),]
  
  # Get non-boundary biased regions
  PositiveNonbias <- data.frame(Region = PositiveRegions[ !(PositiveRegions %in% T1[,2]) & !(PositiveRegions %in% GappedT1[ GappedT1$Keep == 1,2])])
  NegativeNonbias <- data.frame(Region = NegativeRegions[ !(NegativeRegions %in% T2[,2]) & !(NegativeRegions %in% GappedT2[ GappedT2$Keep == 1,2])])
  
  # Stich together T1 regions
  TADs <- data.frame()
  for (i in 1:(nrow(T1)-1)) {
    TADs[i,1] <- T1[i,2]
    TADs[i,2] <- T1[(i+1),2]
  }
  
  # Remove TADs that doesnt cross a T2 boundary or multiple
  TADs$Keep <- 0
  TADs$T2 <- 0
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,2]
    Tmp <- data.frame(T2 = T2[ T2[,2] >= Min & T2[,2] <= Max,])
    if (nrow(Tmp) == 1) { 
      TADs[i,"Keep"] <- 1
      TADs[i,"T2"] <- Tmp[1,1]
    }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have negative bias between first T1 and T2 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,1]
    Max <- TADs[i,4]
    Tmp <- data.frame(NegativeNonbias = NegativeNonbias[ NegativeNonbias > Min & NegativeNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Remove TADs that have positive bias between T2 and second T1 boundary
  for (i in 1:nrow(TADs)) {
    Min <- TADs[i,4]
    Max <- TADs[i,2]
    Tmp <- data.frame(PositiveNonbias = PositiveNonbias[ PositiveNonbias > Min & PositiveNonbias <= Max,])
    if (nrow(Tmp) > 0 ) { TADs[i,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate TAD boundaries
  TADs$chr <- as.character(Chr)
  for (i in 1:nrow(TADs)) {
    Tmp <- Subset[ Subset$Region == TADs[i,1],]
    TADs[i,"start"] <- Tmp[ which.max(Tmp[,5]),3]
    Tmp <- Subset[ Subset$Region == (TADs[i,2]-1),]
    TADs[i,"end"] <- Tmp[ which.min(Tmp[,5]),3]+50000
  }
  
  # Filter TADs in low coverage areas (+/- 20 bins)
  Subset$Forward <- Subset[2,3] - Subset[1,3]
  for (nr in 2:nrow(Subset)) { Subset[nr,"Forward"] <- Subset[nr,3] - Subset[(nr-1),3] }
  
  for (m in 1:nrow(TADs)) {
    Start <- as.numeric(TADs[m,6])-50000
    End <- as.numeric(TADs[m,7])+50000
    Tmp <- Subset
    Tmp <- Tmp[ Tmp[,3] >= Start & Tmp[,3] < End,]
    if (max(Tmp$Forward) > 5000) { TADs[m,"Keep"] <- 0 }
  }
  TADs <- TADs[ TADs$Keep == 1, ]
  
  # Annotate T1 boundaries
  T1$chr <- as.character(Chr)
  for (i in 1:nrow(T1)) {
    Tmp <- Subset[ Subset$Region == T1[i,2],]
    T1[i,"boundary"] <- min(Tmp[,3])
  }
  
  # Finalize TADs
  TADs <- TADs[,c(5,6,7)]
  
  # Combine all TADs
  All <- rbind(All, TADs)
  
  # Save T1 Boundaries
  T1 <- T1[,c(4,5)]
  Bounds <- rbind(Bounds, T1)
}

# Rename the TAD calls
Cortex_TADs <- All
Cortex_Bounds <- Bounds

## Clean up
rm(list=setdiff(ls(),c("Zero_TADs","Four_TADs","Two_TADs","ESC_TADs","Cortex_TADs","Zero_Bounds","Four_Bounds","Two_Bounds","ESC_Bounds","Cortex_Bounds")))

## Import data and convert to Z-scores
# D0
Zero <- read.table("C:/Users/jgsm/Desktop/D0.DI", quote="\"", comment.char="", skip=1)
Mean <- mean(Zero$V4)
SD <- sd(Zero$V4)
Zero$V4 <- (Zero$V4 - Mean)/SD
rm(Mean)
rm(SD)

# H4
Four <- read.table("C:/Users/jgsm/Desktop/H4.DI", quote="\"", comment.char="", skip=1)
Mean <- mean(Four$V4)
SD <- sd(Four$V4)
Four$V4 <- (Four$V4 - Mean)/SD
rm(Mean)
rm(SD)

# D2
Two <- read.table("C:/Users/jgsm/Desktop/D2.DI", quote="\"", comment.char="", skip=1)
Mean <- mean(Two$V4)
SD <- sd(Two$V4)
Two$V4 <- (Two$V4 - Mean)/SD
rm(Mean)
rm(SD)

# ESC
ESC <- read.table("C:/Users/jgsm/Desktop/ESC.DI", quote="\"", comment.char="", skip=1)
Mean <- mean(ESC$V4)
SD <- sd(ESC$V4)
ESC$V4 <- (ESC$V4 - Mean)/SD
rm(Mean)
rm(SD)

# Cortex
Cortex <- read.table("C:/Users/jgsm/Desktop/Cortex.DI", quote="\"", comment.char="", skip=1)
Mean <- mean(Cortex$V4)
SD <- sd(Cortex$V4)
Cortex$V4 <- (Cortex$V4 - Mean)/SD
rm(Mean)
rm(SD)

### Overlap boundaries and group by identification of T1 boundary
## Combine all boundaries
Bounds <- rbind(Zero_Bounds, Four_Bounds, Two_Bounds, ESC_Bounds, Cortex_Bounds)

## Filter boundaries within +/- 25 kb
# Setup
Bounds <- Bounds[ order(Bounds$chr, Bounds$boundary),]
Bounds$Region <- 0
Counter <- 1
Bounds$ID <- seq(1, nrow(Bounds), by = 1)

# Loop
for (i in 1:nrow(Bounds)) {
  Bounds[i,"Region"] <- Counter
  CurrChr <- Bounds[i,1]
  CurrBoundary <- Bounds[i,2]
  CurrID <- Bounds[i,4]
  Tmp <- Bounds[ Bounds$chr == CurrChr & abs(Bounds$boundary - CurrBoundary) <= 25000,]
  if (max(Tmp$ID) == CurrID) { Counter <- Counter + 1 }
}

# Keep center of each region
Bounds$Collaps <- paste(Bounds$chr, Bounds$boundary, Bounds$Region, sep="-")
Bounds <- Bounds[ duplicated(Bounds$Collaps)==F,]
Center <- aggregate(Bounds$boundary, by=list(Bounds$Region), FUN="mean")
colnames(Center) <- c("Region","boundary")
Bounds <- merge(Bounds, Center, by="Region")
Bounds <- Bounds[ duplicated(Bounds$Region)==F,]
Bounds <- Bounds[,c(2,6)]
colnames(Bounds) <- c("chr","boundary")

### Group by cell-type boundaries
## Combine 3T3-L1 T1s
T3_Bounds <- rbind(Zero_Bounds,Four_Bounds,Two_Bounds)

## Group IDs:
# Group 1 == 3T3-L1
# Group 2 == ESC 
# Group 3 == Cortex 
# Group 23 == ESC + Cortex
# Group 12 == 3T3-L1 + ESC
# Group 13 == 3T3-L1 + Cortex
# Group 123 == 3T3-L1 + ESC + Cortex
Bounds$Group <- 0

for (i in 1:nrow(Bounds)) {
  CurrChr <- Bounds[i,1]
  CurrBoundary <- Bounds[i,2]
  TestT3 <- T3_Bounds[ T3_Bounds$chr == CurrChr,]
  TestT3 <- TestT3[ abs(TestT3$boundary - CurrBoundary) <= 50000,]
  TestESC <- ESC_Bounds[ ESC_Bounds$chr == CurrChr,]
  TestESC <- TestESC[ abs(TestESC$boundary - CurrBoundary) <= 50000,]
  TestCortex <- Cortex_Bounds[ Cortex_Bounds$chr == CurrChr,]
  TestCortex <- TestCortex[ abs(TestCortex$boundary - CurrBoundary) <= 50000,]
  if (nrow(TestT3) > 0 & nrow(TestCortex) == 0 & nrow(TestESC) == 0) { Bounds[i,"Group"] <- 1 }
  if (nrow(TestT3) == 0 & nrow(TestCortex) == 0 & nrow(TestESC) > 0) { Bounds[i,"Group"] <- 2 }
  if (nrow(TestT3) == 0 & nrow(TestCortex) > 0 & nrow(TestESC) == 0) { Bounds[i,"Group"] <- 3 }
  if (nrow(TestT3) == 0 & nrow(TestCortex) > 0 & nrow(TestESC) > 0) { Bounds[i,"Group"] <- 23 }
  if (nrow(TestT3) > 0 & nrow(TestCortex) == 0 & nrow(TestESC) > 0) { Bounds[i,"Group"] <- 12 }
  if (nrow(TestT3) > 0 & nrow(TestCortex) > 0 & nrow(TestESC) == 0) { Bounds[i,"Group"] <- 13 }
  if (nrow(TestT3) > 0 & nrow(TestCortex) > 0 & nrow(TestESC) > 0) { Bounds[i,"Group"] <- 123 }
}

## Clean up
rm(list=setdiff(ls(),c("Zero_TADs","Four_TADs","Two_TADs","ESC_TADs","Cortex_TADs","Zero_Bounds","Four_Bounds","Two_Bounds","ESC_Bounds","Cortex_Bounds","Bounds","Zero","Four","Two","ESC","Cortex")))

## For each boundary get +/- 10 bins
HeatmapZero <- matrix(nrow=nrow(Bounds), ncol=31)
HeatmapFour <- matrix(nrow=nrow(Bounds), ncol=31)
HeatmapTwo <- matrix(nrow=nrow(Bounds), ncol=31)
HeatmapESC <- matrix(nrow=nrow(Bounds), ncol=31)
HeatmapCortex <- matrix(nrow=nrow(Bounds), ncol=31)

for (q in 1:nrow(Bounds)) {
  Start <- as.numeric(Bounds[q,2])
  Chr <- as.character(Bounds[q,1])
  Tmp <- Zero[ Zero$V1 == Chr,]
  Tmp$Dist <- abs(Tmp$V2 - Start)
  HeatmapZero[q,] <- Tmp[c((which.min(Tmp$Dist)-15):(which.min(Tmp$Dist)+15)),4]
  Tmp <- Four[ Four$V1 == Chr,]
  Tmp$Dist <- abs(Tmp$V2 - Start)
  HeatmapFour[q,] <- Tmp[c((which.min(Tmp$Dist)-15):(which.min(Tmp$Dist)+15)),4]
  Tmp <- Two[ Two$V1 == Chr,]
  Tmp$Dist <- abs(Tmp$V2 - Start)
  HeatmapTwo[q,] <- Tmp[c((which.min(Tmp$Dist)-15):(which.min(Tmp$Dist)+15)),4]
  Tmp <- ESC[ ESC$V1 == Chr,]
  Tmp$Dist <- abs(Tmp$V2 - Start)
  HeatmapESC[q,] <- Tmp[c((which.min(Tmp$Dist)-15):(which.min(Tmp$Dist)+15)),4]
  Tmp <- Cortex[ Cortex$V1 == Chr,]
  Tmp$Dist <- abs(Tmp$V2 - Start)
  HeatmapCortex[q,] <- Tmp[c((which.min(Tmp$Dist)-15):(which.min(Tmp$Dist)+15)),4]
}

# Combine the results
Result <- cbind(HeatmapZero, HeatmapFour, HeatmapTwo, HeatmapESC, HeatmapCortex)

# Set boundaries to max out at 1
Result[ Result < -1 ] <- -1
Result[ Result > 1 ] <- 1
Result <- as.data.frame(Result)
Result$ID <- seq(1, nrow(Result), by=1)

# Plot all 7 heatmaps
Bounds$ID <- seq(1, nrow(Bounds), by=1)
Group <- Result[ Result$ID %in% Bounds[ Bounds$Group == 123, "ID"],c(1:155)]
pheatmap(Group, cluster_cols = F, cluster_rows = T, color = colorRampPalette(c("red","white","green"))(40), breaks = seq(-1,1,length.out = 41), show_rownames = F, show_colnames = F)

Group <- Result[ Result$ID %in% Bounds[ Bounds$Group == 12, "ID"],c(1:155)]
pheatmap(Group, cluster_cols = F, cluster_rows = T, color = colorRampPalette(c("red","white","green"))(40), breaks = seq(-1,1,length.out = 41), show_rownames = F, show_colnames = F)

Group <- Result[ Result$ID %in% Bounds[ Bounds$Group == 13, "ID"],c(1:155)]
pheatmap(Group, cluster_cols = F, cluster_rows = T, color = colorRampPalette(c("red","white","green"))(40), breaks = seq(-1,1,length.out = 41), show_rownames = F, show_colnames = F)

Group <- Result[ Result$ID %in% Bounds[ Bounds$Group == 23, "ID"],c(1:155)]
pheatmap(Group, cluster_cols = F, cluster_rows = T, color = colorRampPalette(c("red","white","green"))(40), breaks = seq(-1,1,length.out = 41), show_rownames = F, show_colnames = F)

Group <- Result[ Result$ID %in% Bounds[ Bounds$Group == 1, "ID"],c(1:155)]
pheatmap(Group, cluster_cols = F, cluster_rows = T, color = colorRampPalette(c("red","white","green"))(40), breaks = seq(-1,1,length.out = 41), show_rownames = F, show_colnames = F)

Group <- Result[ Result$ID %in% Bounds[ Bounds$Group == 2, "ID"],c(1:155)]
pheatmap(Group, cluster_cols = F, cluster_rows = T, color = colorRampPalette(c("red","white","green"))(40), breaks = seq(-1,1,length.out = 41), show_rownames = F, show_colnames = F)

Group <- Result[ Result$ID %in% Bounds[ Bounds$Group == 3, "ID"],c(1:155)]
pheatmap(Group, cluster_cols = F, cluster_rows = T, color = colorRampPalette(c("red","white","green"))(40), breaks = seq(-1,1,length.out = 41), show_rownames = F, show_colnames = F)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 1](../Links/Figure1.md)
