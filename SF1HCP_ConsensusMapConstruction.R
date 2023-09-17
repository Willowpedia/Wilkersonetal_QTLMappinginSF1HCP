### DUSTIN G WILKERSON ###
# PRODUCE CONSENSUS LINKAGE MAPS FOR S. PURPUREA 94006 AND 94001
# USING BACKCROSS MARKERS IENTIFIED IN THE SALIX F1 HCP POPULATION 

###### PACKAGES AND FUNCTIONS ###### 

library(qtl)
library(ASMap)
library(stringr)
library(LPmerge)

NarrowBinCts <- function(key) { 
  # adds a column to the key file of freq of that marker in all for maps  
  output <- NULL
  key$MARKDAT <- NA
  for (chr in sort(unique(key$CHR))) {
    pushchr <- NULL
    pullchr <- droplevels(key[which(key$CHR == chr), ])
    for (bin in sort(unique(pullchr$BIN_NEW))) {
      pullbin <- droplevels(pullchr[which(pullchr$BIN_NEW == bin), ])
      if (nrow(pullbin) == 1) {
        pullbin$MARKDAT <- pullbin$MARK
        pushchr <- rbind(pushchr, pullbin)
      } else {
        
        maxfreq <- max(pullbin$FREQ)
        
        if (maxfreq == 1) {
          pullbin <- droplevels(pullbin[grepl("_cM", pullbin$TYPE), ])
          pullbin$MARKDAT <- pullbin$MARK
          pushchr <- rbind(pushchr, pullbin)
        }
        
        if (maxfreq == 2) {
          pullbin <- droplevels(pullbin[unique(c(which(grepl("_cM", 
                                                             pullbin$TYPE)), 
                                                 which(pullbin$FREQ == 2))), ])
          freqck <- pullbin$FREQ[which(grepl("_cM", pullbin$TYPE))]
          cmMark <- pullbin$MARK[which(grepl("_cM", pullbin$TYPE))]
          pullbin$MARKDAT <- cmMark
          if (freqck == maxfreq) { 
            pushchr <- rbind(pushchr, pullbin)
          } else { 
            pullbin <- droplevels(pullbin[which(pullbin$FREQ == 2), ])
            pushchr <- rbind(pushchr, pullbin)
          }
        }
        
        if (maxfreq == 3) {
          pullbin <- droplevels(pullbin[unique(c(which(grepl("_cM", 
                                                             pullbin$TYPE)), 
                                                 which(pullbin$FREQ == 3))), ])
          freqck <- pullbin$FREQ[which(grepl("_cM", pullbin$TYPE))]
          cmMark <- pullbin$MARK[which(grepl("_cM", pullbin$TYPE))]
          pullbin$MARKDAT <- cmMark
          if (freqck == maxfreq) { 
            pushchr <- rbind(pushchr, pullbin)
          } else { 
            pullbin <- droplevels(pullbin[which(pullbin$FREQ == 3), ])
            pushchr <- rbind(pushchr, pullbin)
          }
        }
        
        if (maxfreq == 4) {
          pullbin <- droplevels(pullbin[unique(c(which(grepl("_cM", 
                                                             pullbin$TYPE)), 
                                                 which(pullbin$FREQ == 4))), ])
          freqck <- pullbin$FREQ[which(grepl("_cM", pullbin$TYPE))]
          cmMark <- pullbin$MARK[which(grepl("_cM", pullbin$TYPE))]
          pullbin$MARKDAT <- cmMark
          if (freqck == maxfreq) { 
            pushchr <- rbind(pushchr, pullbin)
          } else { 
            pullbin <- droplevels(pullbin[which(pullbin$FREQ == 4), ])
            pushchr <- rbind(pushchr, pullbin)
          }
        }
      }
    }
    output <- rbind(output, pushchr)
  }
  return(output)
}

FindBinReps <- function(key) { 
  
  # Compares the list of markers that were retained in the linkage map with those dropped 
  # for being co-located with those markers and links these markers together within each map
  
  output <- NULL
  key$BINREP <- "TBD"
  for (chr in unique(key$CHR)) {
    pullchr <- droplevels(key[which(key$CHR == chr), ])
    pushchr <- NULL
    for (bin in unique(pullchr$BIN)) { 
      pullbin <- droplevels(pullchr[which(pullchr$BIN == bin), ])
      if (nrow(pullbin) == 1) { # only 1 mark in the bin
        pullbin$BINREP <- "YES"
        pushchr <- rbind(pushchr, pullbin)
      } else {# more than 1 mark in the bin
        if (sum(is.na(pullbin$OTHERSOLO)) == nrow(pullbin)) { # all marks in bin as only in this map
          pullbin$RAND <- rnorm(nrow(pullbin))
          pullbin$BINREP <- ifelse(pullbin$RAND == max(pullbin$RAND), "YES", "NO")
          pullbin$RAND <- NULL
          pushchr <- rbind(pushchr, pullbin)
        } else {# at least one mark in the bin is found in another map
          if (sum(is.na(pullbin$OTHERSOLO)) == 0) {# all marks in bin are found in other maps
            if (sum(pullbin$OTHERSOLO == "YES") == 1) {# only one mark in bin is solobin in other maps
              pullbin$BINREP <- ifelse(pullbin$OTHERSOLO == "YES", "YES", "NO")
              pushchr <- rbind(pushchr, pullbin)
            } else {# either more than one mark in bin is solobin or none are
              if (sum(pullbin$OTHERSOLO == "NO") == nrow(pullbin)) { # all marks not solobin in other maps
                pullbin$RAND <- rnorm(nrow(pullbin))
                pullbin$BINREP <- ifelse(pullbin$RAND == max(pullbin$RAND), "YES", "NO")
                pullbin$RAND <- NULL
                pushchr <- rbind(pushchr, pullbin)
              } else {# at least two marks are solobin in other maps
                pullbin$RAND <- rnorm(nrow(pullbin))
                pullbin$RAND <- ifelse(pullbin$OTHERSOLO == "NO", -100, pullbin$RAND)
                pullbin$BINREP <- ifelse(pullbin$RAND == max(pullbin$RAND), "YES", "NO")
                pullbin$RAND <- NULL
                pushchr <- rbind(pushchr, pullbin)
              }
            }
          } else {# bin mixed between mark only in this map and those in other maps
            if (sum(pullbin$OTHERSOLO == "YES", na.rm = T) == 1) {# only 1 mark is solobin in other maps 
              pullbin$BINREP <- ifelse(pullbin$OTHERSOLO == "YES", "YES", "NO")
              pullbin$BINREP <- ifelse(is.na(pullbin$BINREP), "NO", pullbin$BINREP)
              pushchr <- rbind(pushchr, pullbin)
            } else {# more than one mark is solobin in other maps
              pullbin$RAND <- rnorm(nrow(pullbin))
              pullbin$RAND <- ifelse((pullbin$OTHERSOLO == "NO" | 
                                        is.na(pullbin$OTHERSOLO)), -100, pullbin$RAND)
              pullbin$BINREP <- ifelse(pullbin$RAND == max(pullbin$RAND), "YES", "NO")
              pullbin$RAND <- NULL
              pushchr <- rbind(pushchr, pullbin)
            }
          }
        }  
      }
    }
    output <- rbind(output, pushchr)
  }
  return(output)
}

FindLGMaxInterval <- function(key1, key2, key3, key4, keyids, pop.size, purp) {
  # runs LP merge with a range of max intervals to find the most appropriate for
  # each consensis map. The output is captured and saves to a csv file 
  
  key1 <- droplevels(key1[which(key1$BINREP == "YES"), c(1, 3, 4)])
  key2 <- droplevels(key2[which(key2$BINREP == "YES"), c(1, 3, 4)])
  key3 <- droplevels(key3[which(key3$BINREP == "YES"), c(1, 3, 4)])
  key4 <- droplevels(key4[which(key4$BINREP == "YES"), c(1, 3, 4)])
  
  rownames(key1) <- NULL
  rownames(key2) <- NULL
  rownames(key3) <- NULL
  rownames(key4) <- NULL
  
  for (chr in unique(key1$CHR)) {
    
    lg1 <- droplevels(key1[which(key1$CHR == chr), ])
    lg1$CHR <- NULL
    
    lg2 <- droplevels(key2[which(key2$CHR == chr), ])
    lg2$CHR <- NULL
    
    lg3 <- droplevels(key3[which(key3$CHR == chr), ])
    lg3$CHR <- NULL
    
    lg4 <- droplevels(key4[which(key4$CHR == chr), ])
    lg4$CHR <- NULL
    
    lglist <- list(lg1, lg2, lg3, lg4)
    names(lglist) <- keyids
    
    if (purp == "94006") {
      capture.output(LPmerge(lglist, max.interval = 1:10, pop.size),
                     file = paste("../LPMerge_MaxInt/LGs_94006/CHR", chr,".csv", sep = ""))
    } 
    if (purp == "94001") {
      capture.output(LPmerge(lglist, max.interval = 1:10, pop.size),
                     file = paste("../LPMerge_MaxInt/LGs_94001/CHR", chr,".csv", sep = ""))
    }
  }
}

GetConsensusLGs <- function(key1, key2, key3, key4, keyids, pop.size, maxint, purp) {
  
  # uses the max interval determined in previous step to generate the consensus
  # postitions for each chromosome in each linkage map 
  
  key1 <- droplevels(key1[which(key1$BINREP == "YES"), c(1, 3, 4)])
  key2 <- droplevels(key2[which(key2$BINREP == "YES"), c(1, 3, 4)])
  key3 <- droplevels(key3[which(key3$BINREP == "YES"), c(1, 3, 4)])
  key4 <- droplevels(key4[which(key4$BINREP == "YES"), c(1, 3, 4)])
  
  rownames(key1) <- NULL
  rownames(key2) <- NULL
  rownames(key3) <- NULL
  rownames(key4) <- NULL
  
  for (chr in unique(key1$CHR)) {
    
    lg1 <- droplevels(key1[which(key1$CHR == chr), ])
    lg1$CHR <- NULL
    lg2 <- droplevels(key2[which(key2$CHR == chr), ])
    lg2$CHR <- NULL
    lg3 <- droplevels(key3[which(key3$CHR == chr), ])
    lg3$CHR <- NULL
    lg4 <- droplevels(key4[which(key4$CHR == chr), ])
    lg4$CHR <- NULL
    
    lglist <- list(lg1, lg2, lg3, lg4)
    names(lglist) <- keyids
    
    if (purp == "94006") {
      capture.output(LPmerge(lglist, max.interval = maxint[chr], pop.size),
                     file = paste("../LPMerge_Consensus/CLG_94006/CHR", 
                                  as.character(chr),"_CONSOLE.csv", sep = ""))
      consen <- LPmerge(lglist, max.interval = maxint[chr], pop.size)[[1]]
      write.csv(consen, paste("../LPMerge_Consensus/CLG_94006/CHR",
                              as.character(chr),"_CONSENSUS.csv", sep = ""), 
                quote = F, row.names = F)
    } 
    
    if (purp == "94001") {
      capture.output(LPmerge(lglist, max.interval = maxint[chr], pop.size),
                     file = paste("../LPMerge_Consensus/CLG_94001/CHR",
                                  as.character(chr),"_CONSOLE.csv", sep = ""))
      consen <- LPmerge(lglist, max.interval = maxint[chr], pop.size)[[1]]
      write.csv(consen, paste("../LPMerge_Consensus/CLG_94001/CHR", 
                              as.character(chr),"_CONSENSUS.csv", sep = ""), 
                quote = F, row.names = F)
    }
  }
}

RebuildConsensus <- function(dir) {
  
  # Each chromosome is treated seperatly. The function brings them all together

  con <- list.files(dir, pattern = "_CONSENSUS.csv")
  output <- NULL
  for (chr in con) {
    pullchr <- read.csv(paste(dir, "/", chr, sep = ""))
    colnames(pullchr) <- c("MARK", "CONcM", "AxP", "IxP", "SxP", "VxP")
    pullchr$CHR <- as.numeric(str_replace(word(chr, 1, sep = "_"), "CHR", ""))
    pullchr <- pullchr[ , c(1, 7, 2:6)]
    output <- rbind(output, pullchr)
  }
  output <- output[order(output$CHR, output$CONcM) , ]
  return(output)
}

RefineConsensusMap <- function(condf) { 
  
  output <- NULL
  
  for (chr in unique(condf$CHR)) {
    
    pullchr <- droplevels(condf[which(condf$CHR == chr), ])
    
    for (bin in unique(pullchr$CONcM)) { 
      
      pullbin <- droplevels(pullchr[which(pullchr$CONcM == bin), ])
      
      if (nrow(pullbin) == 1) { 
        
        output <- rbind(output, pullbin)
        
      } else { 
        
        if (nrow(pullbin) > 1) { 
          
          pullbin <- pullbin[which(pullbin$FREQ == max(pullbin$FREQ, na.rm = T)), ]
          
          output <- rbind(output, pullbin)
          
        }
      }
    }
  }
  
  return(output)
}

###### MAXIMIZE COSEGREGATING MARKERS AND CHOOSING BIN REPS ###### 

setwd("CompMaps_94006")

PxU <- read.csv("PxU_94006_FINAL.csv")
PxS <- read.csv("PxS_94006_FINAL.csv")
PxV <- read.csv("PxV_94006_FINAL.csv")
PxK <- read.csv("PxK_94006_FINAL.csv")

PxU_key <- read.csv("PxU_94006_CoLoBins.csv")
PxS_key <- read.csv("PxS_94006_CoLoBins.csv")
PxV_key <- read.csv("PxV_94006_CoLoBins.csv")
PxK_key <- read.csv("PxK_94006_CoLoBins.csv")

table(table(c(colnames(PxU)[-1], colnames(PxS)[-1], 
              colnames(PxV)[-1], colnames(PxK)[-1])))
#    1    2    3    4 
# 6789 2245  479   82 

table(table(c(PxU_key$MARK, PxS_key$MARK, 
              PxV_key$MARK, PxK_key$MARK)))
#    1    2    3    4 
# 5350 4710 4363 5344 

CoSegCts <- as.data.frame(table(c(PxU_key$MARK, PxS_key$MARK, 
                                  PxV_key$MARK, PxK_key$MARK)))
colnames(CoSegCts) <- c("MARK", "FREQ")
PxU_key <- merge(PxU_key, CoSegCts, by = "MARK", all.x = T, all.y = F)
PxS_key <- merge(PxS_key, CoSegCts, by = "MARK", all.x = T, all.y = F)
PxV_key <- merge(PxV_key, CoSegCts, by = "MARK", all.x = T, all.y = F)
PxK_key <- merge(PxK_key, CoSegCts, by = "MARK", all.x = T, all.y = F)

PxU_key <- NarrowBinCts(PxU_key)
PxS_key <- NarrowBinCts(PxS_key)
PxV_key <- NarrowBinCts(PxV_key)
PxK_key <- NarrowBinCts(PxK_key)

table(table(c(PxU_key$MARK, PxS_key$MARK, 
              PxV_key$MARK, PxK_key$MARK)))
#    1    2    3    4 
# 3575 1139  313 5344 

PxU_key <- PxU_key[order(PxU_key$CHR, PxU_key$BIN_NEW), ]
PxS_key <- PxS_key[order(PxS_key$CHR, PxS_key$BIN_NEW), ]
PxV_key <- PxV_key[order(PxV_key$CHR, PxV_key$BIN_NEW), ]
PxK_key <- PxK_key[order(PxK_key$CHR, PxK_key$BIN_NEW), ]

rownames(PxU_key) <- seq.int(1, nrow(PxU_key))
rownames(PxS_key) <- seq.int(1, nrow(PxS_key))
rownames(PxV_key) <- seq.int(1, nrow(PxV_key))
rownames(PxK_key) <- seq.int(1, nrow(PxK_key))

PxU_key <- PxU_key[ , c(1, 9, 2, 5, 7, 8)]
PxS_key <- PxS_key[ , c(1, 9, 2, 5, 7, 8)]
PxV_key <- PxV_key[ , c(1, 9, 2, 5, 7, 8)]
PxK_key <- PxK_key[ , c(1, 9, 2, 5, 7, 8)]

colnames(PxU_key)[5] <- "BIN"
colnames(PxS_key)[5] <- "BIN"
colnames(PxV_key)[5] <- "BIN"
colnames(PxK_key)[5] <- "BIN"

CoSegCts <- as.data.frame(table(c(PxU_key$MARK, PxS_key$MARK, 
                                  PxV_key$MARK, PxK_key$MARK)))
colnames(CoSegCts) <- c("MARK", "FREQ")

PxU_key$FREQ <- unlist(lapply(seq(1, nrow(PxU_key)), function(x) { 
  CoSegCts$FREQ[which(CoSegCts$MARK == PxU_key$MARK[x])]}))
PxS_key$FREQ <- unlist(lapply(seq(1, nrow(PxS_key)), function(x) { 
  CoSegCts$FREQ[which(CoSegCts$MARK == PxS_key$MARK[x])]}))
PxV_key$FREQ <- unlist(lapply(seq(1, nrow(PxV_key)), function(x) { 
  CoSegCts$FREQ[which(CoSegCts$MARK == PxV_key$MARK[x])]}))
PxK_key$FREQ <- unlist(lapply(seq(1, nrow(PxK_key)), function(x) { 
  CoSegCts$FREQ[which(CoSegCts$MARK == PxK_key$MARK[x])]}))

PxU_key$SOLOBIN <- unlist(lapply(seq(1, nrow(PxU_key)), function(x) {
  binsize <- nrow(PxU_key[unique(which((PxU_key$CHR == PxU_key$CHR[x]) & 
                                         (PxU_key$BIN == PxU_key$BIN[x]))), ])
  if (binsize == 1) {"YES"} else {"NO"}
}))
PxS_key$SOLOBIN <- unlist(lapply(seq(1, nrow(PxS_key)), function(x) {
  binsize <- nrow(PxS_key[unique(which((PxS_key$CHR == PxS_key$CHR[x]) & 
                                         (PxS_key$BIN == PxS_key$BIN[x]))), ])
  if (binsize == 1) {"YES"} else {"NO"}
}))
PxV_key$SOLOBIN <- unlist(lapply(seq(1, nrow(PxV_key)), function(x) {
  binsize <- nrow(PxV_key[unique(which((PxV_key$CHR == PxV_key$CHR[x]) & 
                                         (PxV_key$BIN == PxV_key$BIN[x]))), ])
  if (binsize == 1) {"YES"} else {"NO"}
}))
PxK_key$SOLOBIN <- unlist(lapply(seq(1, nrow(PxK_key)), function(x) {
  binsize <- nrow(PxK_key[unique(which((PxK_key$CHR == PxK_key$CHR[x]) & 
                                         (PxK_key$BIN == PxK_key$BIN[x]))), ])
  if (binsize == 1) {"YES"} else {"NO"}
}))

pxu_solo <- PxU_key[ , c(1, 7)]
pxs_solo <- PxS_key[ , c(1, 7)]
pxv_solo <- PxV_key[ , c(1, 7)]
pxk_solo <- PxK_key[ , c(1, 7)]

colnames(pxu_solo)[2] <- "SOLO_PXU"
colnames(pxs_solo)[2] <- "SOLO_PXS"
colnames(pxv_solo)[2] <- "SOLO_PXV"
colnames(pxk_solo)[2] <- "SOLO_PXK"

SOLO <- merge(pxu_solo,
              merge(pxs_solo, 
                    merge(pxv_solo, pxk_solo, 
                          by = "MARK", all = T), 
                    by = "MARK", all = T),
              by = "MARK", all = T)

SOLO$NOTPxU <- ifelse((SOLO$SOLO_PXS == "YES" | 
                         SOLO$SOLO_PXV == "YES" |
                         SOLO$SOLO_PXK == "YES"), "YES", "NO")
SOLO$NOTPxS <- ifelse((SOLO$SOLO_PXU == "YES" | 
                         SOLO$SOLO_PXV == "YES" |
                         SOLO$SOLO_PXK == "YES"), "YES", "NO")
SOLO$NOTPxV <- ifelse((SOLO$SOLO_PXS == "YES" | 
                         SOLO$SOLO_PXU == "YES" |
                         SOLO$SOLO_PXK == "YES"), "YES", "NO")
SOLO$NOTPxK <- ifelse((SOLO$SOLO_PXS == "YES" | 
                         SOLO$SOLO_PXV == "YES" |
                         SOLO$SOLO_PXU == "YES"), "YES", "NO")

PxU_key$OTHERSOLO <- unlist(lapply(seq(1, nrow(PxU_key)), function(x) {
  SOLO$NOTPxU[which(SOLO$MARK == PxU_key$MARK[x])]}))
PxS_key$OTHERSOLO <- unlist(lapply(seq(1, nrow(PxS_key)), function(x) {
  SOLO$NOTPxS[which(SOLO$MARK == PxS_key$MARK[x])]}))
PxV_key$OTHERSOLO <- unlist(lapply(seq(1, nrow(PxV_key)), function(x) {
  SOLO$NOTPxV[which(SOLO$MARK == PxV_key$MARK[x])]}))
PxK_key$OTHERSOLO <- unlist(lapply(seq(1, nrow(PxK_key)), function(x) {
  SOLO$NOTPxK[which(SOLO$MARK == PxK_key$MARK[x])]}))

rm(pxk_solo, pxu_solo, pxv_solo, pxs_solo, SOLO, CoSegCts)

PxU_key <- FindBinReps(PxU_key)
PxS_key <- FindBinReps(PxS_key)
PxV_key <- FindBinReps(PxV_key)
PxK_key <- FindBinReps(PxK_key)

table(table(c(PxU_key$MARK[which(PxU_key$BINREP == "YES")], 
              PxS_key$MARK[which(PxS_key$BINREP == "YES")], 
              PxV_key$MARK[which(PxV_key$BINREP == "YES")], 
              PxK_key$MARK[which(PxK_key$BINREP == "YES")])))
#    1    2    3    4 
# 4325 1592 1105  555 

write.csv(PxU_key, "PxU_94006_BinReps.csv", 
          quote = F, row.names = F)
write.csv(PxS_key, "PxS_94006_BinReps.csv", 
          quote = F, row.names = F)
write.csv(PxV_key, "PxV_94006_BinReps.csv", 
          quote = F, row.names = F)
write.csv(PxK_key, "PxK_94006_BinReps.csv", 
          quote = F, row.names = F)

###### LPMERGE - FINDING MAX.INT ###### 

setwd("CompMaps_94006")

PxU <- read.csv("PxU_94006_FINAL.csv")
PxS <- read.csv("PxS_94006_FINAL.csv")
PxV <- read.csv("PxV_94006_FINAL.csv")
PxK <- read.csv("PxK_94006_FINAL.csv")

PxU_key <- read.csv("PxU_94006_BinReps.csv")
PxS_key <- read.csv("PxS_94006_BinReps.csv")
PxV_key <- read.csv("PxV_94006_BinReps.csv")
PxK_key <- read.csv("PxK_94006_BinReps.csv")

keyids <- c("PxU", "PxS", "PxV", "PxK")
pop.size <- c((nrow(PxU) - 2), (nrow(PxS) - 2), (nrow(PxV) - 2), (nrow(PxK) - 2))

FindLGMaxInterval(PxU_key, PxS_key, PxV_key, PxK_key, keyids, pop.size, "94006")

###### LPMERGE - CONSENSUS LINKAGE GROUPS ###### 

setwd("CompMaps_94006")

PxU <- read.csv("PxU_94006_FINAL.csv")
PxS <- read.csv("PxS_94006_FINAL.csv")
PxV <- read.csv("PxV_94006_FINAL.csv")
PxK <- read.csv("PxK_94006_FINAL.csv")

PxU_key <- read.csv("PxU_94006_BinReps.csv")
PxS_key <- read.csv("PxS_94006_BinReps.csv")
PxV_key <- read.csv("PxV_94006_BinReps.csv")
PxK_key <- read.csv("PxK_94006_BinReps.csv")

keyids <- c("PxU", "PxS", "PxV", "PxK")
pop.size <- c((nrow(PxU) - 2), (nrow(PxS) - 2), (nrow(PxV) - 2), (nrow(PxK) - 2))
maxint <- c(3, 1, 2, 2, 1, 2, 1, 2, 5, 7, 6, 4, 3, 4, 3, 2, 3, 3, 5)

GetConsensusLGs(PxU_key, PxS_key, PxV_key, PxK_key, keyids, pop.size, maxint, "94006")

###### RECONSTRUCTING THE CONSENSUS MAP ######

setwd("CompMaps_94006")

PxU <- read.csv("PxU_94006_FINAL.csv")
PxS <- read.csv("PxS_94006_FINAL.csv")
PxV <- read.csv("PxV_94006_FINAL.csv")
PxK <- read.csv("PxK_94006_FINAL.csv")

PxU_key <- read.csv("PxU_94006_BinReps.csv")
PxS_key <- read.csv("PxS_94006_BinReps.csv")
PxV_key <- read.csv("PxV_94006_BinReps.csv")
PxK_key <- read.csv("PxK_94006_BinReps.csv")

dir <- "../LPMerge_Consensus/CLG_94006"
consensus <- RebuildConsensus(dir)

PxU_key$CONBINREP <- ifelse(PxU_key$MARK %in% 
                              unique(c(PxU_key$MARK[which(PxU_key$BINREP == "YES")],
                                       PxS_key$MARK[which(PxS_key$BINREP == "YES")],
                                       PxV_key$MARK[which(PxV_key$BINREP == "YES")],
                                       PxK_key$MARK[which(PxK_key$BINREP == "YES")])), 
                            "YES", "NO")
PxS_key$CONBINREP <- ifelse(PxS_key$MARK %in% 
                              unique(c(PxU_key$MARK[which(PxU_key$BINREP == "YES")],
                                       PxS_key$MARK[which(PxS_key$BINREP == "YES")],
                                       PxV_key$MARK[which(PxV_key$BINREP == "YES")],
                                       PxK_key$MARK[which(PxK_key$BINREP == "YES")])), 
                            "YES", "NO")
PxV_key$CONBINREP <- ifelse(PxV_key$MARK %in% 
                              unique(c(PxU_key$MARK[which(PxU_key$BINREP == "YES")],
                                       PxS_key$MARK[which(PxS_key$BINREP == "YES")],
                                       PxV_key$MARK[which(PxV_key$BINREP == "YES")],
                                       PxK_key$MARK[which(PxK_key$BINREP == "YES")])), 
                            "YES", "NO")
PxK_key$CONBINREP <- ifelse(PxK_key$MARK %in% 
                              unique(c(PxU_key$MARK[which(PxU_key$BINREP == "YES")],
                                       PxS_key$MARK[which(PxS_key$BINREP == "YES")],
                                       PxV_key$MARK[which(PxV_key$BINREP == "YES")],
                                       PxK_key$MARK[which(PxK_key$BINREP == "YES")])), 
                            "YES", "NO")

PxU_key <- droplevels(PxU_key[which(PxU_key$CONBINREP == "YES"), ]) # 5157
PxS_key <- droplevels(PxS_key[which(PxS_key$CONBINREP == "YES"), ]) # 4350
PxV_key <- droplevels(PxV_key[which(PxV_key$CONBINREP == "YES"), ]) # 4635
PxK_key <- droplevels(PxK_key[which(PxK_key$CONBINREP == "YES"), ]) # 5920

PxU_key <- PxU_key[ , c(1:3, 5)]
PxS_key <- PxS_key[ , c(1:3, 5)]
PxV_key <- PxV_key[ , c(1:3, 5)]
PxK_key <- PxK_key[ , c(1:3, 5)]

CoSegCts <- as.data.frame(table(c(PxU_key$MARK, PxS_key$MARK, 
                                  PxV_key$MARK, PxK_key$MARK)))
colnames(CoSegCts) <- c("MARK", "FREQ")
PxU_key$FREQ <- unlist(lapply(seq(1, nrow(PxU_key)), function(x) { 
  CoSegCts$FREQ[which(CoSegCts$MARK == PxU_key$MARK[x])]}))
PxS_key$FREQ <- unlist(lapply(seq(1, nrow(PxS_key)), function(x) { 
  CoSegCts$FREQ[which(CoSegCts$MARK == PxS_key$MARK[x])]}))
PxV_key$FREQ <- unlist(lapply(seq(1, nrow(PxV_key)), function(x) { 
  CoSegCts$FREQ[which(CoSegCts$MARK == PxV_key$MARK[x])]}))
PxK_key$FREQ <- unlist(lapply(seq(1, nrow(PxK_key)), function(x) { 
  CoSegCts$FREQ[which(CoSegCts$MARK == PxK_key$MARK[x])]}))

rm(CoSegCts)

PxU <- PxU[ -c(1:2), ]
rownames(PxU) <- PxU$GENOTYPE
PxU$GENOTYPE <- NULL
PxU <- t(PxU)
PxU <- cbind(rownames(PxU), PxU)
rownames(PxU) <- NULL
colnames(PxU)[1] <- "MARK"

PxS <- PxS[ -c(1:2), ]
rownames(PxS) <- PxS$GENOTYPE
PxS$GENOTYPE <- NULL
PxS <- t(PxS)
PxS <- cbind(rownames(PxS), PxS)
rownames(PxS) <- NULL
colnames(PxS)[1] <- "MARK"

PxV <- PxV[ -c(1:2), ]
rownames(PxV) <- PxV$GENOTYPE
PxV$GENOTYPE <- NULL
PxV <- t(PxV)
PxV <- cbind(rownames(PxV), PxV)
rownames(PxV) <- NULL
colnames(PxV)[1] <- "MARK"

PxK <- PxK[ -c(1:2), ]
rownames(PxK) <- PxK$GENOTYPE
PxK$GENOTYPE <- NULL
PxK <- t(PxK)
PxK <- cbind(rownames(PxK), PxK)
rownames(PxK) <- NULL
colnames(PxK)[1] <- "MARK"

PxU_key <- merge(PxU_key, PxU, by.x = "MARKDAT", by.y = "MARK", all = T)
PxU_key <- PxU_key[ order(PxU_key$CHR, PxU_key$BIN), ]
PxS_key <- merge(PxS_key, PxS, by.x = "MARKDAT", by.y = "MARK", all = T)
PxS_key <- PxS_key[ order(PxS_key$CHR, PxS_key$BIN), ]
PxV_key <- merge(PxV_key, PxV, by.x = "MARKDAT", by.y = "MARK", all = T)
PxV_key <- PxV_key[ order(PxV_key$CHR, PxV_key$BIN), ]
PxK_key <- merge(PxK_key, PxK, by.x = "MARKDAT", by.y = "MARK", all = T)
PxK_key <- PxK_key[ order(PxK_key$CHR, PxK_key$BIN), ]

PxU_key <- merge(consensus[ , c(1:3)], PxU_key, by = c("MARK", "CHR"), all.x = F, all.y = T)
PxU_key <- PxU_key[ order(PxU_key$CHR, PxU_key$CONcM), ]
PxS_key <- merge(consensus[ , c(1:3)], PxS_key, by = c("MARK", "CHR"), all.x = F, all.y = T)
PxS_key <- PxS_key[ order(PxS_key$CHR, PxS_key$CONcM), ]
PxV_key <- merge(consensus[ , c(1:3)], PxV_key, by = c("MARK", "CHR"), all.x = F, all.y = T)
PxV_key <- PxV_key[ order(PxV_key$CHR, PxV_key$CONcM), ]
PxK_key <- merge(consensus[ , c(1:3)], PxK_key, by = c("MARK", "CHR"), all.x = F, all.y = T)
PxK_key <- PxK_key[ order(PxK_key$CHR, PxK_key$CONcM), ]

conmap <- merge(PxU_key[ , -c(4:5)],
                merge(PxS_key[ , -c(4:5)], 
                      merge(PxV_key[ , -c(4:5)], PxK_key[ , -c(4:5)],
                            by = c("MARK", "CHR", "CONcM", "FREQ"), all = T), 
                      by = c("MARK", "CHR", "CONcM", "FREQ"), all = T), 
                by = c("MARK", "CHR", "CONcM", "FREQ"), all = T)

conmap <- conmap[order(conmap$CHR, conmap$CONcM, -conmap$FREQ) , ]

conmap <- RefineConsensusMap(conmap)

table(conmap[which(conmap$CHR == 1), ]$FREQ)
#   1   2   3   4 
# 117  61   9 232 

table(conmap[which(conmap$CHR == 2), ]$FREQ)
#   1   2   3   4 
# 111  35   9 288

table(conmap[which(conmap$CHR == 3), ]$FREQ)
#  1   2   3   4 
# 72  23   4 174 

table(conmap[which(conmap$CHR == 4), ]$FREQ)
#   1   2   3   4 
# 100  29   9 212 

table(conmap[which(conmap$CHR == 5), ]$FREQ)
#   1   2   3   4 
# 123  38  16 277 

table(conmap[which(conmap$CHR == 6), ]$FREQ)
#   1   2   3   4 
# 125  32  16 269 

table(conmap[which(conmap$CHR == 7), ]$FREQ)
#  1   2   3   4 
# 92  39  11 161 

table(conmap[which(conmap$CHR == 8), ]$FREQ)
#   1   2   3   4 
# 108  47   8 190 

table(conmap[which(conmap$CHR == 9), ]$FREQ)
#  1   2   3   4 
# 80  30  14 151 

table(conmap[which(conmap$CHR == 10), ]$FREQ)
#   1   2   3   4 
# 118  60  14 241 

table(conmap[which(conmap$CHR == 11), ]$FREQ)
#  1   2   3   4 
# 61  33   7 145 

table(conmap[which(conmap$CHR == 12), ]$FREQ)
# 1   2   3   4 
# 93  30  13 109

table(conmap[which(conmap$CHR == 13), ]$FREQ)
#  1   2   3   4 
# 88  21  17 169 

table(conmap[which(conmap$CHR == 14), ]$FREQ)
#  1   2   3   4 
# 60  24   3 153 

table(conmap[which(conmap$CHR == 15), ]$FREQ)
#  1   2   3   4 
# 85  28  13 152

table(conmap[which(conmap$CHR == 16), ]$FREQ)
#   1   2   3   4 
# 151  58  15 331 

table(conmap[which(conmap$CHR == 17), ]$FREQ)
#  1   2   3   4 
# 82  32  11 133 

table(conmap[which(conmap$CHR == 18), ]$FREQ)
#  1   2   3   4 
# 73  44   4 188 

table(conmap[which(conmap$CHR == 19), ]$FREQ)
#  1  2  3  4 
# 67 35 12 97 

consensus <- droplevels(consensus[which(consensus$MARK %in% conmap$MARK), ]) 

write.csv(conmap, "../94006_ConsensusMap.csv", 
          quote = F, row.names = F)
write.csv(consensus, "../94006_ComponentPos.csv", 
          quote = F, row.names = F)

