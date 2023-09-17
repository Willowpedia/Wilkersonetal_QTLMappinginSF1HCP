### DUSTIN G WILKERSON
### Phenotypic Analysis and QTL Mapping of the Salix F1HCP ### 

### Please reach out to Dustin if you have any questions: dgw65@cornell.edu

###### PACKAGES AND FUNCTIONS ###### 

library(stringr)
library(agricolae)
library(ggplot2)
library(MASS)
library(conover.test)
library(corrplot)
library(gridExtra)
library(VCA)
library(ASMap)
library(qtl)
library(circlize)
library(plyr)
library(ggpattern)
library(viridis)

DoBCT <- function(popdat, hsdat, fsdat) {
  # This function will take the untransformaed phenotype data and perform box-cox
  # transformation at the population, half-sib family, and f1 family level for each
  # phenotype. I ran this twice to see which phenotypes could be normalized and set aside 
  # those that couldn't be normalized in (nptrts)
  
  nptrts <- c("WLB_17", "WLB_18", "WLB_19", "PLHC_17",
              "PLHC_19", "PLHN_17", "PLHN_19", "PLHS_17",
              "PLHS_19", "RST_17", "RST_19", "HT_17", "HT_18")
  
  bctrts <- colnames(popdat)[which(!(colnames(popdat) %in% nptrts))][-c(1:9)]
  
  colnames(popdat) <- unlist(lapply(colnames(popdat), function(x) {
    if (x %in% nptrts) {
      colnames(popdat)[which(colnames(popdat) == x)] <- paste("NP_", x, sep = "")
    } else {
      x 
    }
  }))
  
  # population #
  for (trt in bctrts) { 
    bc <- boxcox((popdat[ , trt] + 1) ~ popdat$FAMILY + popdat$REP, lambda = seq(-4, 4, 0.1), plotit = F)
    lambda <- bc$x[which.max(bc$y)]
    lambda <- ifelse(lambda <= -3.5, -4,
                     ifelse(lambda > -3.5 & lambda <= -2.5, -3,
                            ifelse(lambda > -2.5 & lambda <= -1.5, -2,
                                   ifelse(lambda > -1.5 & lambda <= -0.75, -1,
                                          ifelse(lambda > -0.75 & lambda <= -0.25, -0.5,
                                                 ifelse(lambda > -0.25 & lambda < 0.25, 0,
                                                        ifelse(lambda >= 0.25 & lambda < 0.75, 0.5,
                                                               ifelse(lambda >= 0.75 & lambda < 1.5, 1,
                                                                      ifelse(lambda >= 1.5 & lambda < 2.5, 2,
                                                                             ifelse(lambda >= 2.5 & lambda < 3.5, 3,
                                                                                    ifelse(lambda >= 3.5, 4, 99)))))))))))
    if (lambda != 0 & lambda != 99) {
      popdat[ , trt] <- ((popdat[ , trt] + 1) ^ lambda - 1) / lambda
      colnames(popdat)[which(colnames(popdat) == trt)] <- paste("BCT_", trt, sep = "")
    } else {
      if (lambda == 0) {
        popdat[ , trt] <- log2(popdat[ , trt] + 1)
        colnames(popdat)[which(colnames(popdat) == trt)] <- paste("BCT_", trt, sep = "")
      } else {
        print(paste("We have a problem at", trt))
      }
    }
  }
  
  # half-sib # 
  hsout <- NULL
  colnames(hsdat) <- unlist(lapply(colnames(hsdat), function(x) {
    if (x %in% nptrts) {
      colnames(hsdat)[which(colnames(hsdat) == x)] <- paste("NP_", x, sep = "")
    } else {
      x 
    }
  }))
  hsdat$REP <- as.factor(hsdat$REP)
  hsdat$COMPAR <- ifelse(is.na(hsdat$COMPAR), "OFFPAR", hsdat$COMPAR)
  hsdat$COMPAR <- as.factor(hsdat$COMPAR)
  for (compar in unique(hsdat$COMPAR)) {
    pullhs <- droplevels(hsdat[which(hsdat$COMPAR == compar), ])
    for (trt in bctrts) { 
      bc <- boxcox((pullhs[ , trt] + 1) ~ pullhs$FAMILY + pullhs$REP, lambda = seq(-4, 4, 0.1), plotit = F)
      lambda <- bc$x[which.max(bc$y)]
      lambda <- ifelse(lambda <= -3.5, -4,
                       ifelse(lambda > -3.5 & lambda <= -2.5, -3,
                              ifelse(lambda > -2.5 & lambda <= -1.5, -2,
                                     ifelse(lambda > -1.5 & lambda <= -0.75, -1,
                                            ifelse(lambda > -0.75 & lambda <= -0.25, -0.5,
                                                   ifelse(lambda > -0.25 & lambda < 0.25, 0,
                                                          ifelse(lambda >= 0.25 & lambda < 0.75, 0.5,
                                                                 ifelse(lambda >= 0.75 & lambda < 1.5, 1,
                                                                        ifelse(lambda >= 1.5 & lambda < 2.5, 2,
                                                                               ifelse(lambda >= 2.5 & lambda < 3.5, 3,
                                                                                      ifelse(lambda >= 3.5, 4, 99)))))))))))
      if (lambda != 0 & lambda != 99) {
        pullhs[ , trt] <- ((pullhs[ , trt] + 1) ^ lambda - 1) / lambda
        colnames(pullhs)[which(colnames(pullhs) == trt)] <- paste("BCT_", trt, sep = "")
      } else {
        if (lambda == 0) {
          pullhs[ , trt] <- log2(pullhs[ , trt] + 1)
          colnames(pullhs)[which(colnames(pullhs) == trt)] <- paste("BCT_", trt, sep = "")
        } else {
          print(paste("We have a problem at", trt))
        }
      }
    }
    hsout <- rbind(hsout, pullhs)
  }
  
  # full-sib # 
  fsout <- NULL
  colnames(fsdat) <- unlist(lapply(colnames(fsdat), function(x) {
    if (x %in% nptrts) {
      colnames(fsdat)[which(colnames(fsdat) == x)] <- paste("NP_", x, sep = "")
    } else {
      x 
    }
  }))
  fsdat$REP <- as.factor(fsdat$REP)
  
  for (fam in unique(fsdat$FAMILY)) {
    pullfs <- droplevels(fsdat[which(fsdat$FAMILY == fam), ])
    for (trt in bctrts) { 
      bc <- boxcox((pullfs[ , trt] + 1) ~ pullfs$REP, lambda = seq(-4, 4, 0.1), plotit = F)
      lambda <- bc$x[which.max(bc$y)]
      lambda <- ifelse(lambda <= -3.5, -4,
                       ifelse(lambda > -3.5 & lambda <= -2.5, -3,
                              ifelse(lambda > -2.5 & lambda <= -1.5, -2,
                                     ifelse(lambda > -1.5 & lambda <= -0.75, -1,
                                            ifelse(lambda > -0.75 & lambda <= -0.25, -0.5,
                                                   ifelse(lambda > -0.25 & lambda < 0.25, 0,
                                                          ifelse(lambda >= 0.25 & lambda < 0.75, 0.5,
                                                                 ifelse(lambda >= 0.75 & lambda < 1.5, 1,
                                                                        ifelse(lambda >= 1.5 & lambda < 2.5, 2,
                                                                               ifelse(lambda >= 2.5 & lambda < 3.5, 3,
                                                                                      ifelse(lambda >= 3.5, 4, 99)))))))))))
      if (lambda != 0 & lambda != 99) {
        pullfs[ , trt] <- ((pullfs[ , trt] + 1) ^ lambda - 1) / lambda
        colnames(pullfs)[which(colnames(pullfs) == trt)] <- paste("BCT_", trt, sep = "")
      } else {
        if (lambda == 0) {
          pullfs[ , trt] <- log2(pullfs[ , trt] + 1)
          colnames(pullfs)[which(colnames(pullfs) == trt)] <- paste("BCT_", trt, sep = "")
        } else {
          print(paste("We have a problem at", trt))
        }
      }
    }
    fsout <- rbind(fsout, pullfs)
  }
  
  return(list(popdat, hsout, fsout))
  
}

GetMeanSeps <- function(popdat, type) {
  # Performs mean seperation analysis splitting the phenotypes based on 
  # transformation flag into non-parametric or parametric methods 
  
  if (type == "NP") { # non-parametric
    
      nonpar <- popdat[ , c(2, 4, 7, which(grepl("NP_", colnames(popdat))))]
      trtnames <- word(colnames(nonpar)[-c(1:3)], start = 2, end = 3, sep = "_")
      
      # full sib # 
      fams <- lapply(4:ncol(nonpar), function(x){
        aggregate(nonpar[ , x] ~ FAMILY, data = nonpar, mean)[ , 2] })
      fams <- as.data.frame(do.call(cbind, fams))
      fams <- round(fams, digits = 3)
      fams$FAMILY <- levels(nonpar$FAMILY)
      colnames(fams) <- c(trtnames, "FAMILY")
      fams <- fams[ , c(ncol(fams), 1:(ncol(fams) - 1))]
      
      fpvs <- lapply(4:ncol(nonpar), function(x){
        kruskal.test(nonpar[ , x] ~ FAMILY, data = nonpar)$p.value })
      
      stats <- rbind(c("FS K-W PVAL", unlist(fpvs)), fams)
      
      # half sib # 
      hsbs <- lapply(4:ncol(nonpar), function(x){
        aggregate(nonpar[ , x] ~ COMPAR, data = nonpar, mean)[ , 2] })
      hsbs <- as.data.frame(do.call(cbind, hsbs))
      hsbs <- round(hsbs, digits = 3)
      hsbs$FAMILY <- levels(nonpar$COMPAR)
      colnames(hsbs) <- c(trtnames, "FAMILY")
      hsbs <- hsbs[ , c(ncol(hsbs), 1:(ncol(hsbs) - 1))]
      
      hpvs <- lapply(4:ncol(nonpar), function(x){
        wilcox.test(nonpar[ , x] ~ COMPAR, data = nonpar)$p.value })
      
      stats <- rbind(stats, c("HS WILCOX PVAL", unlist(hpvs)), hsbs)
      
      # mean and cv # 
      trtmeans <- lapply(4:ncol(nonpar), function(x){
        mean(nonpar[ , x], na.rm = T) })
      stats <- rbind(stats, c("MEAN", round(unlist(trtmeans), 2)))
      
      trtcv <- lapply(4:ncol(nonpar), function(x){
        (sd(nonpar[ , x], na.rm = T) / mean(nonpar[ , x], na.rm = T)) * 100 })
      stats <- rbind(stats, c("CV", round(unlist(trtcv), 2)))
      
      # comparison test # 
      
      fscomp <- conover.test(x = nonpar[ , 4], g = nonpar$FAMILY, method = "hochberg",
                                  kw = F, table = F, list = F)$comparisons
      cipvals <- lapply(4:ncol(nonpar), function(x){
        conover.test(x = nonpar[ , x], g = nonpar$FAMILY, method = "hochberg", 
                     kw = F, table = F, list = F)$P.adjusted })
      
      cipvals <- as.data.frame(do.call(cbind, cipvals))
      cipvals$COMPARISON <- fscomp
      cipvals <- cipvals[ , c(ncol(cipvals), 1:(ncol(cipvals) - 1))]
      colnames(cipvals) <- c("COMPARISONS", trtnames)
      
    return(list(stats, cipvals))
  }
  
  if (type == "BCT") { # parametric
    
    par <- popdat[ , c(2, 4, 7, which(grepl("BCT_", colnames(popdat))))]
    trtnames <- word(colnames(par)[-c(1:3)], start = 2, end = 3, sep = "_")
    
    # full sib # 
    fams <- lapply(4:ncol(par), function(x){
      aggregate(par[ , x] ~ FAMILY, data = par, mean)[ , 2] })
    fams <- as.data.frame(do.call(cbind, fams))
    fams <- round(fams, digits = 3)
    fams$FAMILY <- levels(par$FAMILY)
    colnames(fams) <- c(trtnames, "FAMILY")
    fams <- fams[ , c(ncol(fams), 1:(ncol(fams) - 1))]
    pvls <- lapply(4:ncol(par), function(x) {
      summary(aov(par[ , x] ~ REP + FAMILY, data = par))[[1]]$`Pr(>F)`[2]})
    
    stats <- rbind(c("FS AVOVA PVAL", unlist(pvls)), fams)
    
    # half sib # 
    hsbs <- lapply(4:ncol(par), function(x){
      aggregate(par[ , x] ~ COMPAR, data = par, mean)[ , 2] })
    hsbs <- as.data.frame(do.call(cbind, hsbs))
    hsbs <- round(hsbs, digits = 3)
    hsbs$FAMILY <- levels(par$COMPAR)
    colnames(hsbs) <- c(trtnames, "FAMILY")
    hsbs <- hsbs[ , c(ncol(hsbs), 1:(ncol(hsbs) - 1))]
    hpvs <- lapply(4:ncol(par), function(x){
      t.test(par[ , x] ~ COMPAR, data = par)$p.value })
    
    stats <- rbind(stats, c("HS T-Test PVAL", unlist(hpvs)), hsbs)
    
    # mean and cv # 
    
    trtmeans <- lapply(4:ncol(par), function(x) {
      mean(par[ , x], na.rm = T) })
    stats <- rbind(stats, c("MEAN", round(unlist(trtmeans), 3)))
    
    trtcv <- lapply(4:ncol(par), function(x) {
      (sd(par[ , x], na.rm = T) / mean(par[ , x], na.rm = T)) * 100 })
    stats <- rbind(stats, c("CV", round(unlist(trtcv), 3)))
    
    ### 
    
    comparisons <- rownames(TukeyHSD(aov(par[ , 4] ~ REP + FAMILY, data = par), which = "FAMILY")$FAMILY)
    pvals <- lapply(4:ncol(par), function(x) {
      as.numeric(TukeyHSD(aov(par[ , x] ~ REP + FAMILY, data = par), which = "FAMILY")$FAMILY[ , 4])
      })
    pvals <- as.data.frame(do.call(cbind, pvals))
    pvals$COMPARISON <- comparisons
    pvals <- pvals[ , c(ncol(pvals), 1:(ncol(pvals) - 1))]
    colnames(pvals) <- c("COMPARISONS", trtnames)
    
    return(list(stats, pvals))
  }
}

GetVarCompsAndH2 <- function(hcp, half, full) {
  # extracts variation components and calculates broad-sense heritability 
  nomissingpheno <- function(varcomp) { 
    nona <- NULL
    for (fam in unique(varcomp$FAMILY)) {
      pullfam <- droplevels(varcomp[which(varcomp$FAMILY == fam), ])
      nomiss <- lapply(seq(6, ncol(pullfam)), function(trt) {
        unlist(lapply(seq(1, nrow(pullfam)), function(geno) {
          ifelse(!is.na(pullfam[geno, trt]),
                 pullfam[geno, trt], 
                 round(mean(pullfam[which(pullfam$GENOTYPE == pullfam$GENOTYPE[geno]), trt], na.rm = T), 2))
        }))
      })
      names(nomiss) <- colnames(pullfam)[6:ncol(pullfam)]
      nomiss <- as.data.frame(do.call("cbind", nomiss))
      nomiss <- cbind(pullfam[ , c(1:2)], nomiss)
      nomiss <- lapply(seq(3, ncol(nomiss)), function(trt) {
        unlist(lapply(seq(1, nrow(nomiss)), function(geno) {
          ifelse(!is.nan(nomiss[geno, trt]),
                 nomiss[geno, trt], 
                 round(mean(nomiss[ , trt], na.rm = T), 2))
        }))
      })
      names(nomiss) <- colnames(pullfam)[6:ncol(pullfam)]
      nomiss <- as.data.frame(do.call("cbind", nomiss))
      nomiss <- cbind(pullfam[ , c(1:5)], nomiss)
      nona <- rbind(nona, nomiss)
    }
    return(nona)
  }
  
  hcp$FAMSPEC <- as.factor(ifelse(hcp$FPAR == "94006", as.character(hcp$MPAR), as.character(hcp$FPAR)))
  hcp <- hcp[ , c(2, 4, 7, 51, 8, 10:50)]
  hcp <- nomissingpheno(hcp)
  
  half$FAMSPEC <- as.factor(ifelse(half$FPAR == "94006", as.character(half$MPAR), as.character(half$FPAR)))
  half <- half[ , c(2, 4, 7, 51, 8, 10:50)]
  half <- nomissingpheno(half)
  
  full$FAMSPEC <- as.factor(ifelse(full$FPAR == "94006", as.character(full$MPAR), as.character(full$FPAR)))
  full <- full[ , c(2, 4, 7, 51, 8, 10:50)]
  full <- nomissingpheno(full)
  
  popout <- NULL
  halfout <- NULL
  fullout <- NULL
  
  for (trt in seq(6, ncol(hcp))) {
    
    trtID <- colnames(hcp)[trt]
    pulltrt <- hcp[ , c(1:5, which(colnames(hcp) == trtID))]
    colnames(pulltrt)[6] <- "Y"
    
    vca <- fitVCA(Y ~ REP + COMPAR/FAMSPEC/GENOTYPE, Data = pulltrt)
    vca <- VCAinference(vca, VarVC = T, excludeNeg = F)
    
    out <- as.data.frame(vca$VCAobj$aov.tab)
    out$TERMS <- rownames(out)
    
    Vr <- round(out$VC[which(out$TERMS == "REP")], 4)
    Vcp <- round(out$VC[which(out$TERMS == "COMPAR")], 4)
    Vfs <- round(out$VC[which(out$TERMS == "COMPAR:FAMSPEC")], 4)
    Vf1 <- round(out$VC[which(out$TERMS == "COMPAR:FAMSPEC:GENOTYPE")], 4)
    Vg <- sum(Vcp, Vfs, Vf1)
    Ve <- round(out$VC[which(out$TERMS == "error")], 4)
    Vp <- round(Vg + (Ve/4), 4)
    
    H2 = round((Vg/Vp), 4)
    
    vcs <- list(c(trtID, Vr, Vcp, Vfs, Vf1, Vg, Ve, Vp, H2))
    
    popout <- append(popout, vcs)
    
  }
  
  for (trt in seq(6, ncol(half))) {
    
    trtID <- colnames(half)[trt]
    pulltrt <- half[ , c(1:5, which(colnames(half) == trtID))]
    colnames(pulltrt)[6] <- "Y"
    
    for (halfsib in c("94006", "94001")) { 
      
      pullhalfsib <- droplevels(pulltrt[which(pulltrt$COMPAR == halfsib), ])[ , c(1, 2, 5, 6)]
    
      vca <- fitVCA(Y ~ REP + FAMILY/GENOTYPE, Data = pullhalfsib)
      vca <- VCAinference(vca, VarVC = T, excludeNeg = F)
      
      out <- as.data.frame(vca$VCAobj$aov.tab)
      out$TERMS <- rownames(out)
      
      Vr <- round(out$VC[which(out$TERMS == "REP")], 4)
      Vfam <- round(out$VC[which(out$TERMS == "FAMILY")], 4)
      Vf1 <- round(out$VC[which(out$TERMS == "FAMILY:GENOTYPE")], 4)
      Vg <- sum(Vfam, Vf1)
      Ve <- round(out$VC[which(out$TERMS == "error")], 4)
      Vp <- round(Vg + (Ve/4), 4)
      
      H2 = round((Vg/Vp), 4)
      
      vcs <- list(c(trtID, halfsib, Vr, Vfam, Vf1, Vg, Ve, Vp, H2))
      
      halfout <- append(halfout, vcs)
    }
  }
  
  for (trt in seq(6, ncol(full))) {
    
    trtID <- colnames(full)[trt]
    pulltrt <- full[ , c(1:5, which(colnames(full) == trtID))]
    colnames(pulltrt)[6] <- "Y"
    
    for (fam in unique(as.character(pulltrt$FAMILY))) { 
      
      pullfam <- droplevels(pulltrt[which(pulltrt$FAMILY == fam), ])[ , c(1, 5, 6)]
      
      vca <- fitVCA(Y ~ REP + GENOTYPE, Data = pullfam)
      vca <- VCAinference(vca, VarVC = T, excludeNeg = F)
      
      out <- as.data.frame(vca$VCAobj$aov.tab)
      out$TERMS <- rownames(out)
      
      Vr <- round(out$VC[which(out$TERMS == "REP")], 4)
      Vg <- round(out$VC[which(out$TERMS == "GENOTYPE")], 4)
      Ve <- round(out$VC[which(out$TERMS == "error")], 4)
      Vp <- round(Vg + (Ve/4), 4)
      
      H2 = round(Vg/Vp, 4)
      
      vcs <- list(c(trtID, fam, Vr, Vg, Ve, Vp, H2))
      
      fullout <- append(fullout, vcs)
      
    }
  }
  
  popout <- as.data.frame(do.call("rbind", popout))
  colnames(popout) <- c("PHENOTYPE", "Vr", "Vcp", 
                        "Vfs", "Vf1", "Vg", 
                        "Ve", "Vp", "H2")
  
  halfout <- as.data.frame(do.call("rbind", halfout))
  colnames(halfout) <- c("PHENOTYPE", "HALFSIB", "Vr",
                        "Vfam", "Vf1", "Vg", 
                        "Ve", "Vp", "H2")
  
  fullout <- as.data.frame(do.call("rbind", fullout))
  colnames(fullout) <- c("PHENOTYPE", "FAMILY", "Vr", 
                        "Vg", "Ve", "Vp", "H2")
  
  return(list(popout, halfout, fullout))
}

GetMidPH <- function(df) {
  # Calculates mid-parent heterosis # 
  
  nomissingpheno <- function(df) { 
      nona <- NULL
      for (fam in unique(df$FAMILY)) {
        pullfam <- droplevels(df[which(df$FAMILY == fam), ])
        nomiss <- lapply(seq(8, ncol(pullfam)), function(trt) {
          unlist(lapply(seq(1, nrow(pullfam)), function(geno) {
            ifelse(!is.na(pullfam[geno, trt]),
                   pullfam[geno, trt], 
                   round(mean(pullfam[which(pullfam$GENOTYPE == pullfam$GENOTYPE[geno]), trt], na.rm = T), 2))
          }))
        })
        names(nomiss) <- colnames(pullfam)[8:ncol(pullfam)]
        nomiss <- as.data.frame(do.call("cbind", nomiss))
        nomiss <- cbind(pullfam[ , c(1:2)], nomiss)
        nomiss <- lapply(seq(3, ncol(nomiss)), function(trt) {
          unlist(lapply(seq(1, nrow(nomiss)), function(geno) {
            ifelse(!is.nan(nomiss[geno, trt]),
                   nomiss[geno, trt], 
                   round(mean(nomiss[ , trt], na.rm = T), 2))
          }))
        })
        names(nomiss) <- colnames(pullfam)[8:ncol(pullfam)]
        nomiss <- as.data.frame(do.call("cbind", nomiss))
        nomiss <- cbind(pullfam[ , c(1:7)], nomiss)
        nona <- rbind(nona, nomiss)
      }
      return(nona)
    }

  midout <- NULL
  
  pars <- droplevels(df[which(df$PED == "Parent"), ])
  f1 <- droplevels(df[which(df$PED == "F1"), ])
  f1 <- nomissingpheno(f1)
  
  for (fam in as.character(unique(f1$FAMILY))) {
    
    if (fam == "10X-400") {
      ps <- droplevels(pars[which(pars$GENOTYPE == "94006" | pars$GENOTYPE == "P63"), ])
    }
    if (fam == "11X-407") {
      ps <- droplevels(pars[which(pars$GENOTYPE == "94006" | pars$GENOTYPE == "Jorr"), ])
    }
    if (fam == "12X-421") {
      ps <- droplevels(pars[which(pars$GENOTYPE == "94001" | pars$GENOTYPE == "07-MBG-5027"), ])
    }
    if (fam == "13X-358") {
      ps <- droplevels(pars[which(pars$GENOTYPE == "94006" | pars$GENOTYPE == "04-BN-051"), ])
    }
    if (fam == "13X-426") {
      ps <- droplevels(pars[which(pars$GENOTYPE == "94001" | pars$GENOTYPE == "P336"), ])
    }
    if (fam == "13X-438") {
      ps <- droplevels(pars[which(pars$GENOTYPE == "94006" | pars$GENOTYPE == "04-FF-016"), ])
    }
    if (fam == "13X-440") {
      ps <- droplevels(pars[which(pars$GENOTYPE == "94001" | pars$GENOTYPE == "P295"), ])
    }
    if (fam == "13X-443") {
      ps <- droplevels(pars[which(pars$GENOTYPE == "94001" | pars$GENOTYPE == "P294"), ])
    }
    
    gs <- droplevels(f1[which(f1$FAMILY == fam), ])
    
    
    mph <- lapply(8:48, function(trt) {
      mid <- mean(ps[ , trt], na.rm = T)
      if (mid == 0) { 
        round(unlist(lapply(1:nrow(gs), function(geno) {
          gs[ geno, trt] * 100
        })), 1)
      } else {
        round(unlist(lapply(1:nrow(gs), function(geno) {
          ((gs[ geno, trt] - mid)/mid) * 100
        })), 1)
      }
    })
    mph <- as.data.frame(do.call("cbind", mph))
    mph <- cbind(gs$FAMILY, gs$GENOTYPE, mph)
    colnames(mph) <- colnames(gs)[c(2, 3, 8:48)]
    midout <- rbind(midout, mph)
    
  }
  
  return(midout)
}

DoMidPHSigTest <- function(df) { 
  # Performs a significance test for mid-parent heterosis 
  
  output <- data.frame(list("TRAIT" = colnames(df)[3:ncol(df)]))
  
  for (fam in as.character(unique(df$FAMILY))) { 
    
    family <- droplevels(df[which(df$FAMILY == fam), ])
    
    mean <- unlist(lapply(3:ncol(family), function(x) {
      round(mean(family[ , x], na.rm = T), 1)
    }))
    pvals <- unlist(lapply(3:ncol(family), function(x) {
      t.test(family[ , x])$p.value
    }))
    mn <- unlist(lapply(3:ncol(family), function(x) {
      round(min(family[ , x], na.rm = T), 1)
    }))
    med <- unlist(lapply(3:ncol(family), function(x) {
      round(median(family[ , x], na.rm = T), 1)
    }))
    mx <- unlist(lapply(3:ncol(family), function(x) {
      round(max(family[ , x], na.rm = T), 1)
    }))
    
    res <- data.frame(cbind(colnames(family)[3:ncol(family)], 
                            mn, mx, med, mean, pvals))
    res$pvals <- as.numeric(res$pvals)
    res$siglvl <- ifelse(res$pvals < 0.0001, 0.0001,
                         ifelse(res$pvals >= 0.0001 & res$pvals < 0.001, 0.001, 
                                ifelse(res$pvals >= 0.001 & res$pvals < 0.01, 0.01,
                                       ifelse(res$pvals >= 0.01 & res$pvals < 0.05, 0.05, 0.99))))
    
    colnames(res) <- c("TRAIT",
                       paste(fam, "MIN", sep = "_"),
                       paste(fam, "MAX", sep = "_"),
                       paste(fam, "MEDIAN", sep = "_"),
                       paste(fam, "MEAN", sep = "_"),
                       paste(fam, "PVAL", sep = "_"),
                       paste(fam, "SIGLVL", sep = "_"))
    
    output <- merge(output, res, by = "TRAIT", all = T)
    
  }
  
  return(output)
}

SaveQTLResults <- function(maplist, pop, hs) {
  # for each linkage map, performs QTL mapping and saves a results file with LOD scores 
  
  nomissingpheno <- function(df) { 
    nona <- NULL
    for (fam in unique(df$FAMILY)) {
      pullfam <- droplevels(df[which(df$FAMILY == fam), ])
      nomiss <- lapply(seq(8, ncol(pullfam)), function(trt) {
        unlist(lapply(seq(1, nrow(pullfam)), function(geno) {
          ifelse(!is.na(pullfam[geno, trt]),
                 pullfam[geno, trt], 
                 round(mean(pullfam[which(pullfam$GENOTYPE == pullfam$GENOTYPE[geno]), trt], na.rm = T), 2))
        }))
      })
      names(nomiss) <- colnames(pullfam)[8:ncol(pullfam)]
      nomiss <- as.data.frame(do.call("cbind", nomiss))
      nomiss <- cbind(pullfam[ , c(1:2)], nomiss)
      nomiss <- lapply(seq(3, ncol(nomiss)), function(trt) {
        unlist(lapply(seq(1, nrow(nomiss)), function(geno) {
          ifelse(!is.nan(nomiss[geno, trt]),
                 nomiss[geno, trt], 
                 round(mean(nomiss[ , trt], na.rm = T), 2))
        }))
      })
      names(nomiss) <- colnames(pullfam)[8:ncol(pullfam)]
      nomiss <- as.data.frame(do.call("cbind", nomiss))
      nomiss <- cbind(pullfam[ , c(1:7)], nomiss)
      nona <- rbind(nona, nomiss)
    }
    return(nona)
  }
  
  pop <- droplevels(pop[which(pop$PED == "F1"), ])
  pop <- nomissingpheno(pop)
  
  hs <- droplevels(hs[which(hs$PED == "F1"), ])
  hs <- nomissingpheno(hs)
  
  for (mapID in maplist) {
    
    cat("Starting Map: ", mapID, "\n", sep = "")
    start <- Sys.time()
    Sys.time()
    map <- read.cross("csv", "LinkageMaps", mapID, estimate.map = F, na.strings = c("NA", "-"), 
                      genotypes = c("AA", "AB"), alleles = c("A", "B"))
    fam <- word(mapID, start = 1, sep = "_")
    
    if (fam == "CON") {
      par <- word(word(mapID, start = 2, sep = "_"), start = 1, sep = "[.]")
      mean <- droplevels(hs[which(hs$COMPAR == par), ])
    } else {
      mean <- droplevels(pop[grepl(fam, pop$FAMILY), ])
    }
    
    mean <- mean[ , c(3, 8:48)]
    mean <- droplevels(mean[which(mean$GENOTYPE %in% map$pheno$GENOTYPE), ])
    mean <- mean[match(map$pheno$GENOTYPE, mean$GENOTYPE), ]
    map$pheno <- mean
    
    cgp <- calc.genoprob(map, step = 5, off.end = 0, 
                         error.prob = 0.01, map.function = "kosambi", 
                         stepwidth = "fixed")
    
    qtl <- scanone(cgp, pheno.col = 2, model = "np")
    perm <- scanone(cgp, pheno.col = 2, model = "np", n.perm = 1000)
    
    results <- as.data.frame(list("SUBJECT" = c("PERM05", "PERM10", rownames(qtl)),
                                  "MARKORDER" = c(1:(length(rownames(qtl)) + 2)), 
                                  "CHR" = c(NA, NA, qtl$chr), 
                                  "POS" = c(NA, NA, qtl$pos), 
                                  "NP_WLB_17" = c(summary(perm, alpha = 0.05)[1], 
                                                  summary(perm, alpha = 0.1)[1],
                                                  qtl$lod)))
    
    for (tr in c(3:ncol(map$pheno))) {
      
      if (tr >= 3 & tr <= 10 | tr == 27 | tr == 28 | tr == 33 | tr == 34) { # nonparametric 
        
        qtl <- scanone(cgp, pheno.col = tr, model = "np")
        perm <- scanone(cgp, pheno.col = tr, model = "np", n.perm = 1000)
        
        out <- as.data.frame(list("SUBJECT" = c("PERM05", "PERM10", rownames(qtl)),
                                  "LOD" = c(summary(perm, alpha = 0.05)[1], 
                                            summary(perm, alpha = 0.1)[1],
                                            qtl$lod)))
        colnames(out)[2] <- colnames(map$pheno)[tr]
        
      } else { # parametric
        
        qtl <- cim(cgp, pheno.col = tr, method = "ehk")
        perm <- cim(cgp, pheno.col = tr, method = "ehk", n.perm = 1000)
        
        out <- as.data.frame(list("SUBJECT" = c("PERM05", "PERM10", rownames(qtl)),
                                  "LOD" = c(summary(perm, alpha = 0.05)[1], 
                                            summary(perm, alpha = 0.1)[1],
                                            qtl$lod)))
        colnames(out)[2] <- colnames(map$pheno)[tr]
        
      }
      
      results <- merge(results, out, by = "SUBJECT", all = T)
      
      if (tr %in% c(6, 10, 14, 18, 22, 26, 30, 34, 38)) {
        cat("Run Update for ", mapID, ": ", tr, " of ", ncol(map$pheno), " (", 
            round(((tr - 2)/(ncol(map$pheno) - 2))*100, digits = 2), "%)", "\n", sep = "")
        print(Sys.time() - start)
      }
    }
    
    if (fam == "CON") {
      fileID <- word(mapID, start = 1, sep = "[.]")
    } else {
      fileID <- word(mapID, start = 1, end = 2, sep = "_")
    }
    
    cat("Saving Map: ", mapID, "\n", sep = "")
    print(Sys.time() - start)
    
    write.table(results, 
                file = paste("QTLMapping/Results_", fileID, ".csv", sep = ""), 
                quote = F, sep = ",", na = "NA", row.names = F, col.names = T, eol = "\r")
  }
  
}

GetQTLSummary <- function(reslist) {
  # Reads in all the results files for each map and brings them all together 
  # into a summary file, calculating summary statistics for each significant QTL 
  
  CalcPhys4Pseudo <- function(pseudo) {
    psechr <- res[which(res$CHR == res$CHR[which(res$SUBJECT == pseudo)]), 1:4]
    psechr$PROX <- abs(psechr$POS[which(psechr$SUBJECT == pseudo)] - psechr$POS)
    psecm <- psechr$POS[which(psechr$SUBJECT == pseudo)]
    psechr <- psechr[which(!grepl("loc", psechr$SUBJECT)), ]
    nstcm <- psechr$POS[which.min(psechr$PROX)]
    nstmb <- as.numeric(word(psechr$SUBJECT[which.min(psechr$PROX)], 2, sep = "_"))/1e6
    return((nstmb/nstcm)*psecm)
  }
  sumout <- NULL
  herb <- c("WLB", "PLH")
  arch <- c("DLW", "LA", "SLA", "LL", 
            "LW", "LP", "LF", "LR")
  rust <- c("RST")
  comp <- c("MEANDIAM", "STEMCT", "HT", "PLTAREA", 
            "PLTVOL", "SPAD", "CROWN")
  for (lods in reslist) {
    fam <- word(lods, start = 2, sep = "_")
    par <- str_replace(word(lods, start = 3, sep = "_"), 
                       pattern = ".csv", "")
    res <- read.csv(paste("QTLMapping/", lods, sep = ""))
    res <- res[order(res$MARKORDER), ]
    for (trt in seq(5, 45)) {
      sig10 <- res[which(res$SUBJECT == "PERM10"), trt]
      if (sum(res[ -c(1, 2), trt] >= sig10) > 0) {
        subqtl <- res[ -c(1, 2), ]
        subqtl <- subqtl[which(subqtl[ , trt] >= sig10), ]
        if (length(unique(subqtl$CHR)) == 1) { # qtl on one chr
          deets <- as.data.frame(list("FAM" = fam, 
                                      "PAR" = par, 
                                      "GRP" = ifelse(grepl(paste(herb, collapse = "|"),
                                                           colnames(res)[trt]), "HERB",
                                                     ifelse(grepl(paste(arch, collapse = "|"),
                                                                  colnames(res)[trt]), "LEAFARCH",
                                                            ifelse(grepl(paste(rust, collapse = "|"), 
                                                                         colnames(res)[trt]), "LEAFRUST",
                                                                   ifelse(grepl(paste(comp, collapse = "|"), 
                                                                                colnames(res)[trt]), "YLDCOMP",
                                                                          "PROB")))),
                                      "TRT" = colnames(res)[trt], 
                                      "CHR" = unique(subqtl$CHR),
                                      "CHRID" = paste(fam, par, "_", unique(subqtl$CHR), sep = ""),
                                      "QTLID" = paste(fam, par, "_", unique(subqtl$CHR), "_", 
                                                      str_replace_all(colnames(res)[trt], "_", ""), sep = ""),
                                      "PEAK" = subqtl$SUBJECT[which.max(subqtl[ , trt])],
                                      "PEAKLOD" = subqtl[ , trt][which.max(subqtl[ , trt])], 
                                      "SIG10" = sig10, 
                                      "MIN_cM" = min(subqtl$POS),
                                      "PEAK_cM" = subqtl$POS[which.max(subqtl[ , trt])],
                                      "MAX_cM" = max(subqtl$POS), 
                                      "MIN_Mb" = ifelse(grepl("loc", subqtl$SUBJECT[which.min(subqtl$POS)]),
                                                        CalcPhys4Pseudo(subqtl$SUBJECT[which.min(subqtl$POS)]),
                                                        as.numeric(word(subqtl$SUBJECT[which.min(subqtl$POS)],
                                                                        2, sep = "_"))/1e6), 
                                      "PEAK_Mb" = ifelse(grepl("loc", subqtl$SUBJECT[which.max(subqtl[ , trt])]),
                                                         CalcPhys4Pseudo(subqtl$SUBJECT[which.max(subqtl[ , trt])]),
                                                         as.numeric(word(subqtl$SUBJECT[which.max(subqtl[ , trt])],
                                                                         2, sep = "_"))/1e6), 
                                      "MAX_Mb" = ifelse(grepl("loc", subqtl$SUBJECT[which.max(subqtl$POS)]),
                                                        CalcPhys4Pseudo(subqtl$SUBJECT[which.max(subqtl$POS)]),
                                                        as.numeric(word(subqtl$SUBJECT[which.max(subqtl$POS)],
                                                                        2, sep = "_"))/1e6)))
          sumout <- rbind(sumout, deets)
        } else { # qtl on more than one chr
          for (chr in unique(subqtl$CHR)) {
            subqtlchr <- subqtl[which(subqtl$CHR == chr), ]
            deets <- as.data.frame(list("FAM" = fam, 
                                        "PAR" = par, 
                                        "GRP" = ifelse(grepl(paste(herb, collapse = "|"),
                                                             colnames(res)[trt]), "HERB",
                                                       ifelse(grepl(paste(arch, collapse = "|"),
                                                                    colnames(res)[trt]), "LEAFARCH",
                                                              ifelse(grepl(paste(rust, collapse = "|"), 
                                                                           colnames(res)[trt]), "LEAFRUST",
                                                                     ifelse(grepl(paste(comp, collapse = "|"), 
                                                                                  colnames(res)[trt]), "YLDCOMP",
                                                                            "PROB")))),
                                        "TRT" = colnames(res)[trt], 
                                        "CHR" = unique(subqtlchr$CHR),
                                        "CHRID" = paste(fam, par, "_", unique(subqtlchr$CHR), sep = ""),
                                        "QTLID" = paste(fam, par, "_", unique(subqtlchr$CHR), "_", 
                                                        str_replace_all(colnames(res)[trt], "_", ""), sep = ""),
                                        "PEAK" = subqtlchr$SUBJECT[which.max(subqtlchr[ , trt])],
                                        "PEAKLOD" = subqtl[ , trt][which.max(subqtl[ , trt])], 
                                        "SIG10" = sig10, 
                                        "MIN_cM" = min(subqtlchr$POS),
                                        "PEAK_cM" = subqtlchr$POS[which.max(subqtlchr[ , trt])],
                                        "MAX_cM" = max(subqtlchr$POS), 
                                        "MIN_Mb" = ifelse(grepl("loc", subqtlchr$SUBJECT[which.min(subqtlchr$POS)]),
                                                          CalcPhys4Pseudo(subqtlchr$SUBJECT[which.min(subqtlchr$POS)]),
                                                          as.numeric(word(subqtlchr$SUBJECT[which.min(subqtlchr$POS)],
                                                                          2, sep = "_"))/1e6), 
                                        "PEAK_Mb" = ifelse(grepl("loc", subqtlchr$SUBJECT[which.max(subqtlchr[ , trt])]),
                                                           CalcPhys4Pseudo(subqtlchr$SUBJECT[which.max(subqtlchr[ , trt])]),
                                                           as.numeric(word(subqtlchr$SUBJECT[which.max(subqtlchr[ , trt])],
                                                                           2, sep = "_"))/1e6), 
                                        "MAX_Mb" = ifelse(grepl("loc", subqtlchr$SUBJECT[which.max(subqtlchr$POS)]),
                                                          CalcPhys4Pseudo(subqtlchr$SUBJECT[which.max(subqtlchr$POS)]),
                                                          as.numeric(word(subqtlchr$SUBJECT[which.max(subqtlchr$POS)],
                                                                          2, sep = "_"))/1e6)))
            sumout <- rbind(sumout, deets)
          }
        }
      }
    }
  }
  
  return(sumout)
}

GetQTLLSI <- function(qtlsummary, pop, hs, maplist, reslist) { 
  # calculates the 1.5 lod support interval for significant QTL 
  
  CalcPhys4Pseudo <- function(pseudo) {
    psechr <- res[which(res$CHR == res$CHR[which(res$SUBJECT == pseudo)]), 1:4]
    psechr$PROX <- abs(psechr$POS[which(psechr$SUBJECT == pseudo)] - psechr$POS)
    psecm <- psechr$POS[which(psechr$SUBJECT == pseudo)]
    psechr <- psechr[which(!grepl("loc", psechr$SUBJECT)), ]
    nstcm <- psechr$POS[which.min(psechr$PROX)]
    nstmb <- as.numeric(word(psechr$SUBJECT[which.min(psechr$PROX)], 2, sep = "_"))/1e6
    return((nstmb/nstcm)*psecm)
  }
  NoMissingPheno <- function(df) { 
    nona <- NULL
    for (fam in unique(df$FAMILY)) {
      pullfam <- droplevels(df[which(df$FAMILY == fam), ])
      nomiss <- lapply(seq(8, ncol(pullfam)), function(trt) {
        unlist(lapply(seq(1, nrow(pullfam)), function(geno) {
          ifelse(!is.na(pullfam[geno, trt]),
                 pullfam[geno, trt], 
                 round(mean(pullfam[which(pullfam$GENOTYPE == pullfam$GENOTYPE[geno]), trt], na.rm = T), 2))
        }))
      })
      names(nomiss) <- colnames(pullfam)[8:ncol(pullfam)]
      nomiss <- as.data.frame(do.call("cbind", nomiss))
      nomiss <- cbind(pullfam[ , c(1:2)], nomiss)
      nomiss <- lapply(seq(3, ncol(nomiss)), function(trt) {
        unlist(lapply(seq(1, nrow(nomiss)), function(geno) {
          ifelse(!is.nan(nomiss[geno, trt]),
                 nomiss[geno, trt], 
                 round(mean(nomiss[ , trt], na.rm = T), 2))
        }))
      })
      names(nomiss) <- colnames(pullfam)[8:ncol(pullfam)]
      nomiss <- as.data.frame(do.call("cbind", nomiss))
      nomiss <- cbind(pullfam[ , c(1:7)], nomiss)
      nona <- rbind(nona, nomiss)
    }
    return(nona)
  }
  lsiout <- NULL
  herb <- c("WLB", "PLH")
  arch <- c("DLW", "LA", "SLA", "LL", 
            "LW", "LP", "LF", "LR")
  rust <- c("RST")
  comp <- c("MEANDIAM", "STEMCT", "HT", "PLTAREA", 
            "PLTVOL", "SPAD", "CROWN")
  
  pop <- droplevels(pop[which(pop$PED == "F1"), ])
  pop <- NoMissingPheno(pop)
  
  hs <- droplevels(hs[which(hs$PED == "F1"), ])
  hs <- NoMissingPheno(hs)
  
  for (mapID in maplist) { 
    map <- read.cross("csv", "LinkageMaps", mapID, estimate.map = F, na.strings = c("NA", "-"), 
                      genotypes = c("AA", "AB"), alleles = c("A", "B"))
    fam <- word(mapID, start = 1, sep = "_")

    if (fam == "CON") {
      par <- word(word(mapID, start = 2, sep = "_"), start = 1, sep = "[.]")
      mean <- droplevels(hs[which(hs$COMPAR == par), ])
    } else {
      par <- word(word(mapID, start = 2, sep = "_"), start = 1, sep = "[.]")
      mean <- droplevels(pop[grepl(fam, pop$FAMILY), ])
    }
    
    res <- read.csv(paste("QTLMapping/Results_", fam, "_", par, ".csv", sep = ""))
    
    mean <- mean[ , c(3, 8:48)]
    mean <- droplevels(mean[which(mean$GENOTYPE %in% map$pheno$GENOTYPE), ])
    mean <- mean[match(map$pheno$GENOTYPE, mean$GENOTYPE), ]
    map$pheno <- mean
    
    map <- jittermap(map)
    cgp <- calc.genoprob(map, step = 5, off.end = 0, 
                         error.prob = 0.01, map.function = "kosambi", 
                         stepwidth = "fixed")
    mapqtl <- qtlsummary[which(qtlsum$FAM == fam & qtlsum$PAR == par), ]
    
    for (qtl in seq(1, nrow(mapqtl))) {
      getqtl <- makeqtl(cgp, what = "prob", 
                        chr = mapqtl$CHR[qtl], 
                        pos = mapqtl$PEAK_cM[qtl])
      refqtl <- refineqtl(cgp, qtl = getqtl,
                          pheno.col = which(colnames(cgp$pheno) == mapqtl$TRT[qtl]),
                          method = "hk", model = "normal", verbose = F)
      fit <- fitqtl(cgp, qtl = refqtl,
                    pheno.col = which(colnames(cgp$pheno) == mapqtl$TRT[qtl]),
                    method = "hk", model = "normal")
      lsi <- lodint(refqtl, chr = mapqtl$CHR[qtl], expandtomarkers = T)
      
      deets <- as.data.frame(list("QTLID" = mapqtl$QTLID[qtl], 
                                  "PVE" = round((1 - 10^(-(2/nind(map))*mapqtl$PEAKLOD[qtl]))*100, 2), 
                                  "PEAK1.5" = rownames(lsi)[2],
                                  "PEAK1.5_LOD" = lsi$lod[2],
                                  "PVE1.5" = fit$result.full["Model", "%var"],
                                  "MIN1.5_cM" = lsi$pos[1],
                                  "PEAK1.5_cM" = lsi$pos[2],
                                  "MAX1.5_cM" = lsi$pos[3], 
                                  "MIN1.5_Mb" = ifelse(grepl("loc", rownames(lsi)[1]),
                                                       CalcPhys4Pseudo(rownames(lsi)[1]),
                                                       as.numeric(word(rownames(lsi)[1],
                                                                       2, sep = "_"))/1e6),
                                  "PEAK1.5_Mb" = ifelse(grepl("loc", rownames(lsi)[2]),
                                                        CalcPhys4Pseudo(rownames(lsi)[2]),
                                                        as.numeric(word(rownames(lsi)[2],
                                                                        2, sep = "_"))/1e6),
                                  "MAX1.5_Mb" = ifelse(grepl("loc", rownames(lsi)[3]),
                                                       CalcPhys4Pseudo(rownames(lsi)[3]),
                                                       as.numeric(word(rownames(lsi)[3],
                                                                       2, sep = "_"))/1e6)))
      lsiout <- rbind(lsiout, deets)
    }
  }
  lsisum <- merge(qtlsummary, lsiout, by = "QTLID", all = T)
  lsisum <- lsisum[ , c(2:7, 1, 8:ncol(lsisum))]
  return(lsisum)
}

GetLinkageMaps <- function(qtlmeta, trt, mappath) {
  # Used in generating QTL figures; produces a dataframe of map positions
  # for linkage groups with QTL 
  
  qtlmeta <- droplevels(qtlmeta[which(grepl(trt, qtlmeta$TRT)), ])
  qtlmeta <- cbind(paste(qtlmeta$FAM, qtlmeta$PAR, sep = "_"), qtlmeta)
  colnames(qtlmeta)[1] <- "MAP"
  
  fampar <- unique(qtlmeta$MAP)
  
  maplist <- list.files(mappath)
  maplist <- maplist[which(grepl(paste(fampar, collapse = "|"), maplist))]
  
  chrcatch <- NULL
  
  for (mapID in maplist) { 
    
    map <- read.cross("csv", mappath, mapID, estimate.map = F, na.strings = c("NA", "-"), 
                      genotypes = c("AA", "AB"), alleles = c("A", "B"))
    
    if (grepl("CON", mapID)) { 
      mapqtl <- qtlmeta[which(qtlmeta$MAP == word(mapID, 1, sep = "[.]")), ]
    } else { 
      mapqtl <- qtlmeta[which(qtlmeta$MAP == word(mapID, 1, 2, sep = "_")), ]
    }
    
    for (mapchr in nrow(mapqtl)) { 
      
      getchr <- map$geno[as.numeric(word(mapqtl$CHRID[mapchr], 2, sep = "_"))][[1]]$map
      
      markid <- names(getchr)
      getchr <- as.data.frame(list("MAP" = mapqtl$MAP[mapchr], 
                                   "FAM" = mapqtl$FAM[mapchr], 
                                   "PARENT" = mapqtl$PAR[mapchr], 
                                   "CHR" = mapqtl$CHR[mapchr],
                                   "PLOTCHRID" = mapqtl$CHR[mapchr],
                                   "MARKER" = markid, 
                                   "cM" = as.numeric(getchr)))
      
      chrcatch <- rbind(chrcatch, getchr)
    }
  }
  return(chrcatch)
}

###### REMOVING OUTLIERS AND TRANSFORM ###### 

### NORMALIZING WITH BOX-COX ### 
# large sample size gives normality tests too much power; SW becomes 
# hypersensitive to deviations from normality # 
# better to determine normality visually through hist and qq # 

POP_ROT <- read.csv("CPNoOut_POP.csv", stringsAsFactors = T)
HS_ROT <- read.csv("CPNoOut_HS.csv", stringsAsFactors = T)
FS_ROT <- read.csv("CPNoOut_FS.csv", stringsAsFactors = T)

catch <- DoBCT(POP_ROT, HS_ROT, FS_ROT)

POP_ROT <- catch[[1]]
HS_ROT <- catch[[2]]
FS_ROT <- catch[[3]]

# population check
par(mfrow = c(3, 3), mar = c(7, 4, 2, 1))
for (i in seq.int(10, ncol(POP_ROT))) {
  hist(POP_ROT[, i], 
       main = colnames(POP_ROT[i]), 
       col = "darkorchid3") 
} 

par(mfrow = c(3, 3), mar = c(7, 4, 2, 1))
for (compar in unique(HS_ROT$COMPAR)) {
  pullcompar <- droplevels(HS_ROT[which(HS_ROT$COMPAR == compar), ])
  for (i in seq.int(10, ncol(HS_ROT))) {
    hist(pullcompar[, i], 
         main = paste(compar, colnames(pullcompar)[i], sep = ":"),
         col = "darkorchid3") 
  } 
}

par(mfrow = c(3, 4), mar = c(7, 4, 2, 1))
for (fam in c("10X-400", "11X-407", "12X-421", "13X-443", 
              "13X-440", "13X-426", "13X-358", "13X-438")) {
  pullfam <- droplevels(FS_ROT[which(FS_ROT$FAMILY == fam), ])
  for (i in seq.int(10, ncol(FS_ROT))) {
    hist(pullfam[, i], 
         main = paste(fam, colnames(pullfam)[i], sep = ":"),
         col = "darkorchid3") 
  } 
}

write.table(POP_ROT, "CPPOP_ROT.txt", quote = F, row.names = F, sep = "\t")
write.table(HS_ROT, "CPHS_ROT.txt", quote = F, row.names = F, sep = "\t")
write.table(FS_ROT, "CPFS_ROT.txt", quote = F, row.names = F, sep = "\t")

###### GENOTYPE-LEVEL MEANS ######  

pop <- read.table("CPPOP_ROT.txt", header = T, sep = "\t")
hs <- read.table("CPHS_ROT.txt", header = T, sep = "\t")
fs <- read.table("CPFS_ROT.txt", header = T, sep = "\t")

for (c in 1:9) { 
  pop[ , c] <- as.factor(pop[ , c])
  hs[ , c] <- as.factor(hs[ , c])
  fs[ , c] <- as.factor(fs[ , c])
}

pop_gm <- aggregate.data.frame(x = pop[ , c(10:50)], by = pop[ , c(3, 4, 8)], mean, na.rm = T)
pop_gm <- merge(pop[ , c(3:9)], pop_gm, by = c("PED", "FAMILY", "GENOTYPE"), all = T)
pop_gm <- pop_gm[!duplicated(pop_gm), ] # 1048

hs_gm <- aggregate.data.frame(x = hs[ , c(10:50)], by = hs[ , c(3, 4, 8)], mean, na.rm = T)
hs_gm <- merge(hs[ , c(3:9)], hs_gm, by = c("PED", "FAMILY", "GENOTYPE"), all = T)
hs_gm <- hs_gm[!duplicated(hs_gm), ] # 1048

fs_gm <- aggregate.data.frame(x = fs[ , c(10:50)], by = fs[ , c(3, 4, 8)], mean, na.rm = T)
fs_gm <- merge(fs[ , c(3:9)], fs_gm, by = c("PED", "FAMILY", "GENOTYPE"), all = T)
fs_gm <- fs_gm[!duplicated(fs_gm), ] # 1048

for (c in 8:48) { 
  pop_gm[ , c] <- ifelse(is.nan(pop_gm[ , c]), NA, pop_gm[, c])
  hs_gm[ , c] <- ifelse(is.nan(hs_gm[ , c]), NA, hs_gm[, c])
  fs_gm[ , c] <- ifelse(is.nan(fs_gm[ , c]), NA, fs_gm[, c])
}

write.table(pop_gm, "CPPOP_GLM.txt", quote = F, row.names = F, sep = "\t")
write.table(hs_gm, "CPHS_GLM.txt", quote = F, row.names = F, sep = "\t")
write.table(fs_gm, "CPFS_GLM.txt", quote = F, row.names = F, sep = "\t")

###### CORRELATION ANALYSIS ###### 

pop_gm <- read.table("CPPOP_GLM.txt", header = T, sep = "\t")
hs_gm <- read.table("CPHS_GLM.txt", header = T, sep = "\t")
fs_gm <- read.table("CPFS_GLM.txt", header = T, sep = "\t")

# no parents 
pop_gm <- droplevels(pop_gm[which(pop_gm$PED == "F1"), ]) # 1038
hs_gm <- droplevels(hs_gm[which(hs_gm$PED == "F1"), ]) # 1038
fs_gm <- droplevels(fs_gm[which(fs_gm$PED == "F1"), ]) # 1038

# population # 

pop_gm <- pop_gm[ , -c(1:7)] #only trait data
colnames(pop_gm) <- word(colnames(pop_gm), start = 2, end = 3, sep = "_")

popcorrs <- cor(pop_gm, method = "spearman", use = "complete.obs")
poppvs <- cor.mtest(pop_gm, use = "complete.obs", method = "spearman", exact = FALSE)$p

pdf("POP_corrplot.pdf", width = 8.25, height = 9)

corrplot(popcorrs, p.mat = poppvs, sig.level = (0.05/(41^2)),
         method = "square", type = "upper", 
         insig = "blank", diag = FALSE, tl.col = "black", 
         tl.cex = 0.8, cl.pos = "b", cl.ratio = 0.05, 
         mar = c(0, 0, 2, 0),
         order = "hclust", title = "Salix F1 HCP")

dev.off()

write.csv(popcorrs, "POP_CORRVALUES.csv", quote = F, row.names = T)
write.csv(poppvs, "POP_PVALUES.csv", quote = F, row.names = T)

# half-sib # 

hs_gm <- hs_gm[ , -c(1:5, 7)] # only compar and traits
colnames(hs_gm)[2:42] <- word(colnames(hs_gm)[2:42], start = 2, end = 3, sep = "_")

for (compar in unique(hs_gm$COMPAR)) { 
  
  pullhs <- droplevels(hs_gm[which(hs_gm$COMPAR == compar), ])
  pullhs$COMPAR <- NULL
  
  corrs <- cor(pullhs, method = "spearman", use = "complete.obs")
  pvls <- cor.mtest(pullhs, use = "complete.obs", method = "spearman", exact = FALSE)$p
  
  pdf(paste("HS_", compar, "_corrplot.pdf", sep = ""), width = 8.25, height = 9)
  
  corrplot(corrs, p.mat = pvls, sig.level = (0.05/(41^2)),
           method = "square", type = "upper", 
           insig = "blank", diag = FALSE, tl.col = "black", 
           tl.cex = 0.8, cl.pos = "b", cl.ratio = 0.05, 
           mar = c(0, 0, 2, 0),
           order = "hclust", title = paste(compar, "Half-Sib Family", sep = " "))
  
  dev.off()
  
  write.csv(corrs, paste("HS_", compar, "_CORRVALUES.csv", sep = ""), quote = F, row.names = T)
  write.csv(pvls, paste("HS_", compar, "_PVALUES.csv", sep = ""), quote = F, row.names = T)
  
}

# full-sib # 

fs_gm <- fs_gm[ , -c(1, 3:7)] # only compar and traits
colnames(fs_gm)[2:42] <- word(colnames(fs_gm)[2:42], start = 2, end = 3, sep = "_")

for (family in unique(fs_gm$FAMILY)) { 
  
  pullfs <- droplevels(fs_gm[which(fs_gm$FAMILY == family), ])
  pullfs$FAMILY <- NULL
  
  corrs <- cor(pullfs, method = "spearman", use = "complete.obs")
  pvls <- cor.mtest(pullfs, use = "complete.obs", method = "spearman", exact = FALSE)$p
  
  pdf(paste("CorrelationAnalysis/FS_", family, "_corrplot.pdf", sep = ""), width = 8.25, height = 9)
  
  corrplot(corrs, p.mat = pvls, sig.level = (0.05/(41^2)),
           method = "square", type = "upper", 
           insig = "blank", diag = FALSE, tl.col = "black", 
           tl.cex = 0.8, cl.pos = "b", cl.ratio = 0.05, 
           mar = c(0, 0, 2, 0),
           order = "hclust", title = paste(family, "Full-Sib Family", sep = " "))
  
  dev.off()
  
  write.csv(corrs, paste("CorrelationAnalysis/FS_", family, "_CORRVALUES.csv", sep = ""), quote = F, row.names = T)
  write.csv(pvls, paste("CorrelationAnalysis/FS_",family, "_PVALUES.csv", sep = ""), quote = F, row.names = T)
  
}

###### MEAN SEPARATIONS & VIOLINS ###### 

pop <- read.table("CPPOP_ROT.txt", header = T, sep = "\t")

pop <- droplevels(pop[which(pop$PED == "F1"), ])
pop$REP <- as.factor(pop$REP)
pop$FAMILY <- as.factor(pop$FAMILY)
pop$COMPAR <- as.factor(pop$COMPAR)

nonpar <- GetMeanSeps(pop, "NP")
para <- GetMeanSeps(pop, "BCT")

popdat <- pop

write.csv(nonpar[[1]], "NONPAR_MEANS.csv", quote = F, row.names = F)
write.csv(nonpar[[2]], "NONPAR_COMPARISONS.csv", quote = F, row.names = F)

write.csv(para[[1]], "PARA_MEANS.csv", quote = F, row.names = F)
write.csv(para[[2]], "PARA_COMPARISONS.csv", quote = F, row.names = F)

# VIOLINS # 

pop <- read.table("CPPOP_GLM.txt", header = T, sep = "\t")

for (c in 1:7) { 
  pop[ , c] <- as.factor(pop[ , c])
}

pop$SET <- ifelse(pop$FAMILY == "13X-358" | pop$FAMILY == "04-BN-051", "S.p x S.u",
                  ifelse(pop$FAMILY == "10X-400" | pop$FAMILY == "P63", "S.p x S.s",
                         ifelse(pop$FAMILY == "11X-407" | pop$FAMILY == "Jorr", "S.p x S.v",
                                ifelse(pop$FAMILY == "12X-421" | pop$FAMILY == "07-MBG-5027", "S.v x S.p",
                                       ifelse(pop$FAMILY == "13X-426" | pop$FAMILY == "P336", "S.i x S.p",
                                              ifelse(pop$FAMILY == "13X-438" | pop$FAMILY == "04-FF-016", "S.p x S.k",
                                                     ifelse(pop$FAMILY == "13X-440" | pop$FAMILY == "P295","S.s(P295) x S.p",
                                                            ifelse(pop$FAMILY == "13X-443" | pop$FAMILY == "P294", "S.s(P294) x S.p", "COMPAR"))))))))

purp <- droplevels(pop[which(pop$SET == "COMPAR"), ])
pop <- droplevels(pop[which(pop$SET != "COMPAR"), ])

purp <- rbind(purp, purp, purp, purp)
purp$SET <- c("S.p x S.s", "S.p x S.v", "S.v x S.p", "S.i x S.p", 
              "S.p x S.u", "S.p x S.k", "S.s(P295) x S.p", "S.s(P294) x S.p")

pop <- rbind(purp, pop)
pop$POPCOL <- ifelse(pop$SET == "S.p x S.u", "#006cf3",
                     ifelse(pop$SET == "S.p x S.s", "#00b036",
                            ifelse(pop$SET == "S.s(P295) x S.p", "#00b036",
                                   ifelse(pop$SET == "S.s(P294) x S.p", "#00b036",
                                          ifelse(pop$SET == "S.v x S.p", "#ff2b36",
                                                 ifelse(pop$SET == "S.p x S.v", "#ff2b36",
                                                        ifelse(pop$SET == "S.i x S.p", "#29d3b9",
                                                               ifelse(pop$SET == "S.p x S.k", "#a39e91",
                                                                      ifelse(pop$SET == "COMPAR", "#c764d9", "PROB")))))))))

pop$COMPAR <- ifelse(pop$PED == "F1", as.character(pop$COMPAR), 
                     ifelse(pop$PED == "Parent" & pop$GENOTYPE %in% c("94001", "P294", "P295", "P336", "07-MBG-5027"), "94001", 
                            ifelse(pop$PED == "Parent" & pop$GENOTYPE %in% c("94006", "P63", "04-BN-051", "04-FF-016", "Jorr"), "94006", "PROB")))


fs <- lapply(seq.int(8, 48), function(trt) {
  ggplot(pop[which(pop$PED == "F1"), ],
         aes(x = reorder(SET, -pop[which(pop$PED == "F1"), ][ , trt], mean, na.rm = T), 
             y = pop[which(pop$PED == "F1"), ][ , trt], 
             fill = POPCOL)) + 
         geom_violin(alpha = 0.4) +
         geom_jitter(alpha = 0.25, width = 0.15) +
         stat_summary(fun = mean, na.rm = T, geom = "point", shape = 23, size = 5, fill = "black") +
         stat_summary(fun.data = mean_se, na.rm = T, geom = "errorbar", size = 1, width = 0.8) +
         geom_point(data = pop[which(pop$PED != "F1"), ],
                    aes(y = pop[which(pop$PED != "F1"), ][ , trt], 
                        x = reorder(SET, -pop[which(pop$PED != "F1"), ][ , trt], mean, na.rm = T)),
                    fill = ifelse(pop[which(pop$PED != "F1"), ]$GENOTYPE %in% c("94001", "94006"), "#c764d9", pop[which(pop$PED != "F1"), ]$POPCOL),
                    size = 5, shape = 23, color = "black") + 
         scale_fill_manual(values = c("#006cf3", "#00b036", "#29d3b9", "#a39e91", "#ff2b36")) + 
         labs(x = "FULL-SIB FAMILY", y = colnames(pop)[trt]) + 
         theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                          vjust = 0.5, size = 10, colour = "black"),
               axis.text.y = element_text(size = 10, colour = "black"),
               axis.title = element_text(size = 10, face = "bold"), 
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(fill = NA),
               plot.title = element_text(hjust = 0.5), 
               legend.position = "none")
  })

hs <- lapply(seq.int(8, 48), function(trt) {
  ggplot(pop[which(pop$PED == "F1"), ], 
         aes(x = reorder(COMPAR, -pop[which(pop$PED == "F1"), ][ , trt], mean, na.rm = T), 
             y = pop[which(pop$PED == "F1"), ][ , trt])) + 
         geom_violin(alpha = 0.4, fill = "#c764d9") +
         geom_jitter(alpha = 0.25, width = 0.15) +
         stat_summary(fun = mean, na.rm = T, geom = "point", shape = 23, size = 5, fill = "black") +
         stat_summary(fun.data = mean_se, na.rm = T, geom = "errorbar", size = 1, width = 0.8) + 
         geom_point(data = pop[which(pop$PED != "F1"), ],
                    aes(y = pop[which(pop$PED != "F1"), ][ , trt], 
                        x = reorder(COMPAR, -pop[which(pop$PED != "F1"), ][ , trt], mean, na.rm = T)),
                    fill = ifelse(pop[which(pop$PED != "F1"), ]$GENOTYPE %in% c("94001", "94006"), "#c764d9", pop[which(pop$PED != "F1"), ]$POPCOL),
                    size = 5, shape = 23, color = "black") + 
         scale_fill_manual(values = c("#006cf3", "#00b036", "#29d3b9", "#a39e91", "#ff2b36")) + 
         labs(x = "HALF-SIB FAMILY", y = colnames(pop)[trt]) + 
         theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                          vjust = 0.5, size = 10, colour = "black"),
               axis.text.y = element_text(size = 10, colour = "black"),
               axis.title = element_text(size = 10, face = "bold"), 
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(fill = NA),
               plot.title = element_text(hjust = 0.5), 
               legend.position = "none")
  })

test <- list(fs, hs)
test <- t(do.call(cbind, test))

catch <- NULL
for (c in seq.int(1, ncol(test))) {
  catch <- append(catch, test[ , c])
}

mg <- marrangeGrob(catch, nrow = 2, ncol = 2, widths = c(3, 1), align = "h", layout_matrix = matrix(1:4, 2, 2, TRUE))

ggsave("CP_MeanSeparationViolins.pdf", mg, width = 11, height = 8.5, units = "in")

###### VARIANCE COMPONENTS & HERITABILITY ###### 

pop <- read.table("CPPOP_ROT.txt", header = T, sep = "\t")
hs <- read.table("CPHS_ROT.txt", header = T, sep = "\t")
fs <- read.table("CPFS_ROT.txt", header = T, sep = "\t")

for (c in 1:9) { 
  pop[ , c] <- as.factor(pop[ , c])
  hs[ , c] <- as.factor(hs[ , c])
  fs[ , c] <- as.factor(fs[ , c])
}
rm(c)

pop <- droplevels(pop[which(pop$PED == "F1"), ])
hs <- droplevels(hs[which(hs$PED == "F1"), ])
fs <- droplevels(fs[which(fs$PED == "F1"), ])

catch <- GetVarCompsAndH2(pop, hs, fs)

pop <- catch[[1]]
half <- catch[[2]]
fam <- catch[[3]]

write.csv(pop, "VCH2_POP.csv", 
            quote = F, row.names = F)
write.csv(half, "VCH2_HSIB.csv", 
          quote = F, row.names = F)
write.csv(fam, "VCH2_FSIB.csv", 
          quote = F, row.names = F)

###### MID PARENT HETEROSIS ###### 

pop_gm <- read.table("CPPOP_GLM.txt", header = T, stringsAsFactors = T)
pop_gm$COMPAR <- as.factor(pop_gm$COMPAR)

mph <- GetMidPH(pop_gm)

write.csv(mph, "MidParHeterosis.csv", quote = F, row.names = F)

### Testing Significance ### 

mph <- read.csv("MidParHeterosis.csv", stringsAsFactors = T) 

mphsig <- DoMidPHSigTest(mph)

write.csv(mphsig, "MPHSignificance.csv", quote = F, row.names = F)

###### QTL MAPPING ######

pop_gm <- read.table("CPPOP_GLM.txt", header = T, stringsAsFactors = T)
pop_gm$COMPAR <- as.factor(pop_gm$COMPAR)

hs_gm <- read.table("CPHS_GLM.txt", header = T, stringsAsFactors = T)
hs_gm$COMPAR <- as.factor(hs_gm$COMPAR)

maplist <- list.files("LinkageMaps")

#SaveQTLResults(maplist, pop_gm, hs_gm)

reslist <- list.files("QTLMapping", pattern = ".csv")[-1]

qtlsum <- GetQTLSummary(reslist) # 267! 

qtlsum <- GetQTLLSI(qtlsum, pop_gm, hs_gm, maplist, reslist)

write.table(qtlsum, file = "QTLMapping/QTLSummary.csv",
            quote = F, sep = ",", na = "NA",
            row.names = F, col.names = T, eol = "\r")

# Refining the results # 

qtlsum <- read.csv("QTLMapping/QTLSummary.csv", stringsAsFactors = T)
qtlsum$CHR <- factor(ifelse(qtlsum$CHR <= 9, paste("CHR0", qtlsum$CHR, sep = ""), paste("CHR", qtlsum$CHR, sep = "")))

qtlsum$OGSIZE <- qtlsum$MAX_Mb - qtlsum$MIN_Mb
qtlsum$REFSIZE <- qtlsum$MAX1.5_Mb - qtlsum$MIN1.5_Mb

qtlsum$QTLMIN_Mb <- ifelse(qtlsum$OGSIZE >= qtlsum$REFSIZE, qtlsum$MIN1.5_Mb, qtlsum$MIN_Mb)
qtlsum$QTLPEAK_Mb <- ifelse(qtlsum$OGSIZE >= qtlsum$REFSIZE, qtlsum$PEAK1.5_Mb, qtlsum$PEAK_Mb)
qtlsum$QTLMAX_Mb <- ifelse(qtlsum$OGSIZE >= qtlsum$REFSIZE, qtlsum$MAX1.5_Mb, qtlsum$MAX_Mb)

qtlsum$QTLMIN_cM <- ifelse(qtlsum$OGSIZE >= qtlsum$REFSIZE, qtlsum$MIN1.5_cM, qtlsum$MIN_cM)
qtlsum$QTLPEAK_cM <- ifelse(qtlsum$OGSIZE >= qtlsum$REFSIZE, qtlsum$PEAK1.5_cM, qtlsum$PEAK_cM)
qtlsum$QTLMAX_cM <- ifelse(qtlsum$OGSIZE >= qtlsum$REFSIZE, qtlsum$MAX1.5_cM, qtlsum$MAX_cM)

qtlsum$QTLSIZE_Mb <- qtlsum$QTLMAX_Mb - qtlsum$QTLMIN_Mb 

curated <- qtlsum[ , c(1:9, 17, 35, 29:34)]
curated$QTLSIZE_Mb <- ifelse(curated$QTLSIZE_Mb < 0, 0, curated$QTLSIZE_Mb)

for (c in 11:14) { 
  curated[ , c] <- round(curated[ , c], 6)
}

for (c in 15:17) { 
  curated[ , c] <- round(curated[ , c], 1)
}

curated$PVE <- round(curated$PVE, 1)
curated$PEAKLOD <- round(curated$PEAKLOD, 2)

write.table(curated, file = "QTLMapping/QTLClean.csv",
            quote = F, sep = ",", na = "NA",
            row.names = F, col.names = T, eol = "\r")

curated <- read.csv("QTLMapping/QTLClean.csv", stringsAsFactors = T)

sum(curated$QTLSIZE_Mb == 0) # 51

curated <- droplevels(curated[which(curated$QTLSIZE_Mb != 0), ]) # 216
curated <- droplevels(curated[which(curated$PEAKLOD >= 3), ]) # 215

write.table(curated, file = "QTLMapping/QTLCurated.csv",
            quote = F, sep = ",", na = "NA",
            row.names = F, col.names = T, eol = "\r")


###### QTL CIRCOS PLOT ###### 

qtl <- read.csv("QTLCurated.csv", stringsAsFactors = T)
qtl <- qtl[ , c(1:5, 12:14)]
qtl[ , 6:8] <- (qtl[ , 6:8] * 1000000)


ref <- read.table("SpRefChrBounds.txt", colClasses = c("character", "numeric", "numeric"), sep = "\t")

png("QTLCircos_HeatMap.png", width = 5, height = 5, res = 1500, units = "in", pointsize = 5)

circos.clear()
circos.par("start.degree" = 90, 
           cell.padding = c(0.0008, 0.01, 0.01, 0.01), 
           track.margin = c(0.005, 0.005))

circos.initializeWithIdeogram(ref, plotType = c("axis", "labels"), 
                              axis.labels.cex = 0.7, labels.cex = 1, major.by = 10000000)

viri <- viridis(20)
circos.track(ylim = c(0, 1), track.height = 0.075, bg.col = alpha("grey", 0.5))

SI <- unique(qtl$CHR)[5]
for (SI in unique(qtl$CHR)) { 
  
  pullchr <- droplevels(qtl[which(qtl$CHR == SI), ])
  chrlen <- ref$V3[which(ref$V1 == SI)]
  binmins <- unlist(lapply(c(1:floor(chrlen/5e5), ceiling(chrlen/5e5)), function(x) {
    if (x == ceiling(chrlen/5e5)) {
      ((ceiling(chrlen/5e5) - 1)*5e5) + 1
    } else {
      ((x * 5e5) - 500000) + 1
    }
  }))
  
  binmids <- ceiling(unlist(lapply(c(1:floor(chrlen/5e5), ceiling(chrlen/5e5)), function(x) {
    if (x == ceiling(chrlen/5e5)) {
      ((((ceiling(chrlen/5e5) - 1)*5e5) + 1) + chrlen)/2
    } else {
      ((((x * 5e5) - 500000) + 1) + (x * 5e5))/2
    }
  })))
  
  binmaxs <- ceiling(unlist(lapply(c(1:floor(chrlen/5e5), ceiling(chrlen/5e5)), function(x) {
    if (x == ceiling(chrlen/5e5)) {
      chrlen
    } else {
      x * 5e5
    }
  })))
  
  bins <- as.data.frame(list("BIN" = c(1:floor(chrlen/5e5), ceiling(chrlen/5e5)), 
                             "MIN" = binmins, 
                             "MAX" = binmaxs, 
                             "MID" = binmids))
  
  bins$SCORE <- unlist(lapply(1:nrow(bins), function(x) {
    sum(pullchr$QTLMIN_Mb >= bins$MIN[x] & pullchr$QTLMAX_Mb <= bins$MAX[x] |
        pullchr$QTLMIN_Mb >= bins$MIN[x] & pullchr$QTLMIN_Mb <= bins$MAX[x] & 
          pullchr$QTLMAX_Mb >= bins$MAX[x] | 
        pullchr$QTLMIN_Mb <= bins$MIN[x] & pullchr$QTLMAX_Mb >= bins$MAX[x] |
        pullchr$QTLMIN_Mb <= bins$MIN[x] & 
          pullchr$QTLMAX_Mb >= bins$MIN[x] & pullchr$QTLMAX_Mb <= bins$MAX[x])
    }))
  
  bins$COLOR <- unlist(lapply(1:nrow(bins), function(x) {
    if (bins$SCORE[x] == 0) { 
      alpha("grey", 0.5)
    } else { 
      viri[bins$SCORE[x]]
    }
  }))
  
  circos.rect(sector.index = SI, 
              xleft = bins$MIN,
              xright = bins$MAX,
              ybottom = 0.01,
              ytop = 1.1,
              col = bins$COLOR, 
              border = NA)
  
}

# YIELD COMPONENTS # 
circos.track(ylim = c(0, 20), track.height = 0.2, bg.col = alpha("deepskyblue2", 0.2))
for (SI in unique(qtl$CHR)) { 
  pullchr <- droplevels(qtl[which(qtl$CHR == SI & qtl$GRP == "YLDCOMP"), ])
  qtlct <- nrow(pullchr)
  if (qtlct > 0) {
    x0 <- pullchr$QTLMIN_Mb
    x1 <- pullchr$QTLMAX_Mb
    y0 <-  unlist(lapply(seq(1, qtlct), function(x) {
      (19.5/qtlct*x)
      }))
    y1 <- y0
    circos.segments(sector.index = SI,
                    x0 = x0, y0 = y0, 
                    x1 = x1, y1 = y1, lwd = 1.5, 
                    col = ifelse(pullchr$FAM == "CON", "darkorchid4", "grey70"))
    circos.points(sector.index = SI, 
                  x = pullchr$QTLPEAK_Mb, 
                  y = y1, cex = 0.5, pch = 19, 
                  col = ifelse(pullchr$FAM == "358", "#006cf3", 
                               ifelse(pullchr$FAM %in% c("400", "440", "443"), "#00b036",  
                                      ifelse(pullchr$FAM == "426", "#29d3b9",
                                             ifelse(pullchr$FAM == "438", "#a39e91",
                                                    ifelse(pullchr$FAM %in% c("407", "421"), "#ff2b36",
                                                           ifelse(pullchr$FAM == "CON", "#c764d9", "white")))))))
  }
}

# LEAF ARCHITECTURE # 
circos.track(ylim = c(0, 14), track.height = 0.15, bg.col = alpha("forestgreen", 0.2))
for (SI in unique(qtl$CHR)) {
  pullchr <- droplevels(qtl[which(qtl$CHR == SI & qtl$GRP == "LEAFARCH"), ])
  qtlct <- nrow(pullchr)
  if (qtlct > 0) {
    x0 <- pullchr$QTLMIN_Mb
    x1 <- pullchr$QTLMAX_Mb
    y0 <-  unlist(lapply(seq(1, qtlct), function(x) {
      (13.5/qtlct*x)
    }))
    y1 <- y0
    circos.segments(sector.index = SI,
                    x0 = x0, y0 = y0, 
                    x1 = x1, y1 = y1, lwd = 1.5, 
                    col = ifelse(pullchr$FAM == "CON", "darkorchid4", "grey70"))
    circos.points(sector.index = SI, 
                  x = pullchr$QTLPEAK_Mb, 
                  y = y1, cex = 0.5, pch = 19,
                  col = ifelse(pullchr$FAM == "358", "#006cf3", 
                               ifelse(pullchr$FAM %in% c("400", "440", "443"), "#00b036",  
                                      ifelse(pullchr$FAM == "426", "#29d3b9",
                                             ifelse(pullchr$FAM == "438", "#a39e91",
                                                    ifelse(pullchr$FAM %in% c("407", "421"), "#ff2b36",
                                                           ifelse(pullchr$FAM == "CON", "#c764d9", "white")))))))
  }
}

# HERBIVORY # 
circos.track(ylim = c(0, 8), track.height = 0.1, bg.col = alpha("darkgoldenrod3", 0.2))
for (SI in unique(qtl$CHR)) {
  pullchr <- droplevels(qtl[which(qtl$CHR == SI & qtl$GRP == "HERB"), ])
  qtlct <- nrow(pullchr)
  if (qtlct > 0) {
    x0 <- pullchr$QTLMIN_Mb
    x1 <- pullchr$QTLMAX_Mb
    y0 <-  unlist(lapply(seq(1, qtlct), function(x) {
      (7.85/qtlct*x)
    }))
    y1 <- y0
    circos.segments(sector.index = SI,
                    x0 = x0, y0 = y0, 
                    x1 = x1, y1 = y1, lwd = 1.5, 
                    col = ifelse(pullchr$FAM == "CON", "darkorchid4", "grey70"))
    circos.points(sector.index = SI, 
                  x = pullchr$QTLPEAK_Mb, 
                  y = y1, cex = 0.5, pch = 19, 
                  col = ifelse(pullchr$FAM == "358", "#006cf3", 
                               ifelse(pullchr$FAM %in% c("400", "440", "443"), "#00b036",  
                                      ifelse(pullchr$FAM == "426", "#29d3b9",
                                             ifelse(pullchr$FAM == "438", "#a39e91",
                                                    ifelse(pullchr$FAM %in% c("407", "421"), "#ff2b36",
                                                           ifelse(pullchr$FAM == "CON", "#c764d9", "white")))))))
  }
}

# LEAF RUST # 
circos.track(ylim = c(0, 5), track.height = 0.075, bg.col = alpha("firebrick3", 0.2))
for (SI in unique(qtl$CHR)) { 
  pullchr <- droplevels(qtl[which(qtl$CHR == SI & qtl$GRP == "LEAFRUST"), ])
  qtlct <- nrow(pullchr)
  if (qtlct > 0) {
    x0 <- pullchr$QTLMIN_Mb
    x1 <- pullchr$QTLMAX_Mb
    y0 <-  unlist(lapply(seq(1, qtlct), function(x) {
      (4.75/qtlct*x)
    }))
    y1 <- y0
    circos.segments(sector.index = SI,
                    x0 = x0, y0 = y0, 
                    x1 = x1, y1 = y1, lwd = 1.5, 
                    col = ifelse(pullchr$FAM == "CON", "darkorchid4", "grey70"))
    circos.points(sector.index = SI, 
                  x = pullchr$QTLPEAK_Mb, 
                  y = y1, cex = 0.5, pch = 19, 
                  col = ifelse(pullchr$FAM == "358", "#006cf3", 
                               ifelse(pullchr$FAM %in% c("400", "440", "443"), "#00b036",  
                                      ifelse(pullchr$FAM == "426", "#29d3b9",
                                             ifelse(pullchr$FAM == "438", "#a39e91",
                                                    ifelse(pullchr$FAM %in% c("407", "421"), "#ff2b36",
                                                           ifelse(pullchr$FAM == "CON", "#c764d9", "white")))))))
  }
}

dev.off()

###### QTL FIGURES: PREP & STATS ######

### QTL figures were generated in R and finished with adjustments in photoshop 

qtl <- read.csv("QTLCurated.csv", stringsAsFactors = T)

ref <- read.table("SpRefChrBounds.txt", colClasses = c("character", "numeric", "numeric"), sep = "\t")
colnames(ref) <- c("CHR", "START", "END")
ref$END_Mb <- ref$END/1e6

THEME <- theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, 
                                          size = 11, colour = "black"),
               axis.text.y = element_text(size = 11, colour = "black"),
               axis.title = element_text(size = 12, face = "bold"), 
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(fill = NA),
               plot.title = element_text(hjust = 0.5), 
               legend.position = "none")

###### LEAF RUST ######

rstqtl <- droplevels(qtl[which(grepl("RST", qtl$TRT)), ])

rstref <- ref
rstref <- droplevels(rstref[which(rstref$CHR %in% rstqtl$CHR), ])
rstref$REFX <- c(15, 43, 71, 99, 127, 155, 183, 222.5, 250.5)

mbticks <- NULL
for (rc in 1:nrow(rstref)) {
  ticks <- as.data.frame(list("REFX" = rstref$REFX[rc],
                              "MBTS" = seq.int(0, round_any(rstref$END_Mb[rc], 0.5, f = floor), by = 0.5)))
  ticks$XEND <- ifelse(grepl("[.]5", as.character(ticks$MBTS)), 2.5, 5)
  mbticks <- rbind(mbticks, ticks)
}

rstqtl <- rstqtl[order(rstqtl$CHR, rstqtl$QTLMIN_Mb), ]
rstqtl$QTLX <- c(27.5, 55.5, 83.5, 111.5, 139.5, 167.5,
                 195.5, 195.5, 207, 
                 235, 
                 263, 274.5, 286, 263, 274.5)
rstqtl$TRT_YR <- word(rstqtl$TRT, 2, 3, sep = "_")

physgrad <- NULL
for (chr in unique(rstqtl$CHR)) {
  pullchr <- droplevels(rstqtl[which(rstqtl$CHR == chr), ])
  chrend <- round(rstref$END_Mb[which(rstref$CHR == chr)], 2)
  atphys <- 0
  splits <- NULL
  while (atphys != chrend) {
    if (atphys == 0) {
      atphys <- round(pullchr$QTLMIN_Mb[which.min(pullchr$QTLMIN_Mb)], 2)
      splits <- append(splits, atphys)
    } else {
      pos <- round(c(pullchr$QTLMIN_Mb, pullchr$QTLMAX_Mb, chrend), 2)
      pos <- sort(pos)
      atphys <- pos[which(pos > atphys)][1]
      splits <- append(splits, atphys)
    }
  }
  chrbks <- as.data.frame(list("CHR" = pullchr$CHR[1], 
                               "XMIN" = rep(rstref$REFX[which(rstref$CHR == chr)], length(splits)), 
                               "XMAX" = rep(rstref$REFX[which(rstref$CHR == chr)] + 5, length(splits)), 
                               "YMIN" = c(0, (splits[-length(splits)] + 0.01)), 
                               "YMAX" = c(splits)))
  physgrad <- rbind(physgrad, chrbks)
}

physgrad$DENS <- c(0, 1, 0, 
                   0, 1, 0,
                   0, 1, 0,
                   0, 1, 0,
                   0, 1, 0,
                   0, 1, 0,
                   0, 1, 0, 1, 2, 1, 0, 
                   0, 1, 0,
                   0, 1, 3, 2, 1, 0, 1, 2, 1, 0)

denscale <- as.data.frame(list("CHR" = rep("SCALE", 4), 
                               "XMIN" = c(155, (155 + 33.75), (155 + (33.75 * 2)), (155 + (33.75 * 3))), 
                               "XMAX" = c((155 + 33.75), (155 + (33.75 * 2)), (155 + (33.75 * 3)), 290), 
                               "YMIN" = rep(19.2, 4), "YMAX" = rep(19.8, 4), 
                               "DENS" = c(0, 1, 2, 3))) 

physgrad <- rbind(physgrad, denscale)
physgrad$COLTEST <- as.factor(ifelse(physgrad$DENS == 0, NA, physgrad$DENS))

ggplot( ) + 
  # set ploting dimensions # 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 295), 
                     position = "top", breaks = (rstref$REFX + 2.25), 
                     labels = rstref$CHR) +
  scale_y_reverse(expand = c(0.01, 0), limits = c(20, -0.01),
                  breaks = seq.int(0, 20, by = 1),
                  labels = seq.int(0, 20, by = 1)) +
  # qtl density heatmap # 
  geom_rect(aes(xmin = physgrad$XMIN, xmax = physgrad$XMAX, 
                ymin = physgrad$YMIN, ymax = physgrad$YMAX, 
                alpha = 0.55, fill = physgrad$COLTEST)) + 
  scale_fill_viridis(discrete = T, begin = .4) +
  # physical Mb scale for chromosomes # 
  geom_segment(aes(x = c(rstref$REFX, rstref$REFX + 5), xend =  c(rstref$REFX, rstref$REFX + 5),
                   y = c(rstref$START, rstref$START), yend = c(rstref$END_Mb, rstref$END_Mb)),
               linewidth = 1) + 
  geom_segment(aes(x = mbticks$REFX, xend = (mbticks$REFX + mbticks$XEND),
                   y = mbticks$MBTS, yend = mbticks$MBTS),
               linewidth = 1, alpha = 0.5) + 
  # plot qtl: sex inner border with species fill # 
  geom_rect(aes(xmin = (rstqtl$QTLX - 3.5), xmax = (rstqtl$QTLX + 3.5),
              ymin = (rstqtl$QTLMIN_Mb), ymax = (rstqtl$QTLMAX_Mb)),
          linewidth = 1.75, alpha = 0.33, 
          fill = ifelse(rstqtl$PAR == "04BN051", "#006cf3",
                        ifelse(rstqtl$PAR %in% c("P63", "P295", "P294"), "#00b036",
                               ifelse(rstqtl$PAR == "P336", "#29d3b9",
                                      ifelse(rstqtl$PAR == "04FF016", "#a39e91",
                                             ifelse(rstqtl$PAR %in% c("Jorr", "07MBG5027"), "#ff2b36",
                                                    ifelse(rstqtl$PAR %in% c("94006", "94001"), "#c764d9", "white")))))),
          color = ifelse(rstqtl$PAR %in% c("94001", "P63", "Jorr", "04FF016", "04BN051"), "cyan4",
                         ifelse(rstqtl$PAR %in% c("94006", "P294", "P295", "07MBG5027", "P336"), "deeppink3", "black"))) +
  # plot qtl: family outer border # 
  geom_rect(aes(xmin = (rstqtl$QTLX - 5), xmax = (rstqtl$QTLX + 5),
                ymin = rstqtl$QTLMIN_Mb, ymax = rstqtl$QTLMAX_Mb),
            linewidth = 1.5, fill = NA,
            color = ifelse(rstqtl$FAM == "358", "#006cf3", 
                           ifelse(rstqtl$FAM %in% c("400", "440", "443"), "#00b036",
                                  ifelse(rstqtl$FAM == "426", "#29d3b9",
                                         ifelse(rstqtl$FAM == "438", "#a39e91",
                                                ifelse(rstqtl$FAM %in% c("407", "421"), "#ff2b36",
                                                       ifelse(rstqtl$FAM == "CON", "#c764d9", "white"))))))) +
  # plot peak marker # 
  geom_segment(aes(x = (rstqtl$QTLX - 2.5), xend = (rstqtl$QTLX + 2.5),
                   y = (rstqtl$QTLPEAK_Mb), yend = (rstqtl$QTLPEAK_Mb)),
               linewidth = 1.5, col = "black") +
  # qtl text labels # 
  geom_text(aes(x = rstqtl$QTLX[-c(7)] - 0.5, 
                y = rstqtl$QTLMAX_Mb[-c(7)] + 0.2,
                label = rstqtl$TRT_YR[-c(7)]), 
            size = 3.75, fontface = "bold", angle = 90, hjust = "right") + 
  # creating legend: qtl denstiy scale # 
  geom_rect(aes(xmin = 155, xmax = 290, ymin = 19.2, ymax = 19.8), linewidth = 0.5, color = "black", fill = NA) + 
  # creating legend: parent species boxes # 
  geom_rect(aes(xmin = c(155, 155, (155 + 45), (155 + 45), (155 + (45 * 2)), (155 + (45 * 2))),  
                xmax = c(170, 170, (170 + 45), (170 + 45), (170 + (45 * 2)), (170 + (45 * 2))),
                ymin = rep(c(18, 18.6), 3),
                ymax = rep(c(18.4, 19), 3)), 
            linewidth = 0.5, color = "black", alpha = 0.33,
            fill = c("#29d3b9", "#00b036", "#a39e91", "#ff2b36", "#c764d9", "#006cf3")) + 
  # creating legend: family boxes # 
  geom_rect(aes(xmin = c(155, 155, 155, (155 + 45), (155 + 45), (155 + 45), (155 + (45 * 2)), (155 + (45 * 2)), (155 + (45 * 2))),  
                xmax = c(170, 170, 170, (170 + 45), (170 + 45), (170 + 45), (170 + (45 * 2)), (170 + (45 * 2)), (170 + (45 * 2))),
                ymin = rep(c(16.2, 16.8, 17.4), 3), 
                ymax = rep(c(16.6, 17.2, 17.8), 3)), 
            linewidth = 1.75, fill = NA, 
            color = c("#a39e91", "#00b036", "#ff2b36", "#006cf3", 
                      "#c764d9", 
                      "#29d3b9", "#00b036", "#00b036", "#ff2b36")) + 
  # creating legend: sex inner border # 
  geom_rect(aes(xmin = c((155 + 1.5), (155 + 45 + 1.5), 155, (155 + 45)), 
                xmax = c((170 - 1.5), (170 + 45 - 1.5), 170, (170 + 45)),
                ymin = c(15.6, 15.6, 15.6, 15.6), 
                ymax = c(16, 16, 16, 16)), 
            linewidth = 1.75, fill = NA, 
            color = c("deeppink3", "cyan4", "black", "black")) + 
  # creating legend: text labels # 
  geom_text(aes(x = c(120, (170 + 3), (170 + 45 + 3), 
                      (170 + 3), (170 + 45 + 3), (170 + (45 * 2) + 3),
                      120, (170 + 3), (170 + 45 + 3), (170 + (45 * 2) + 3),
                      (170 + 3), (170 + 45 + 3), (170 + (45 * 2) + 3), 
                      (170 + 3), (170 + 45 + 3), (170 + (45 * 2) + 3),
                      120,
                      (170 + 3), (170 + 45 + 3), (170 + (45 * 2) + 3), 
                      120, 171.85, 205.65, 239.45, 273.25), 
                y = c(15.75, 15.75, 15.75, 
                      16.35, 16.35, 16.35, 
                      16.95, 16.95, 16.95, 16.95, 
                      17.55, 17.55, 17.55, 
                      18.15, 18.15, 18.15,
                      18.45,
                      18.75, 18.75, 18.75, 
                      19.5, 19.5, 19.5, 19.5, 19.5), 
                label = c("Parent Sex", "Female", "Male", 
                          "S.p x S.k", "S.p x S.u", "S.s(4) x S.p", 
                          "Family", "S.p x S.s", "S.p Cons", "S.s(5) x S.p", 
                          "S.p x S.v", "S.i x S.p", "S.v x S.p", 
                          "S. inte.", "S. kori.", "S. purp.",
                          "Par. Species", "S. such.", "S. vimen", "S. uden", 
                          "QTL Density", "0", "1", "2", "3")),
            size = 4, fontface = "bold", hjust = "left") + 
  # axis titles and theme # 
  labs(x = "Chromosomes with Leaf Rust Severity QTL" , y = "Physical Distance in 94006 v5.1 Reference Genome (Mb)") + 
  THEME

rm(chrbks, mbticks, physgrad, pullchr, rstqtl, rstref, ticks, atphys, chr, chrend, pos, rc, splits)

###### HERBIVORY ######

herbqtl <- droplevels(qtl[which(qtl$GRP == "HERB"), ])

herbref <- ref
herbref <- droplevels(herbref[which(herbref$CHR %in% herbqtl$CHR), ])
herbref <- droplevels(herbref[which(herbref$CHR %in% c("CHR02", "CHR04", "CHR05", "CHR09", 
                                                       "CHR12", "CHR13", "CHR19")), ])

herbref$REFX <- c(15, 77.5, 128.5, 191, 288, 385, 424.5)

mbticks <- NULL
for (rc in 1:nrow(herbref)) {
  ticks <- as.data.frame(list("REFX" = herbref$REFX[rc],
                              "MBTS" = seq.int(0, round_any(herbref$END_Mb[rc], 0.5, f = floor), by = 0.5)))
  ticks$XEND <- ifelse(grepl("[.]5", as.character(ticks$MBTS)), 2.5, 5)
  mbticks <- rbind(mbticks, ticks)
}

herbqtl <- herbqtl[order(herbqtl$CHR, herbqtl$QTLMIN_Mb), ]
herbqtl <- droplevels(herbqtl[which(herbqtl$CHR %in% c("CHR02", "CHR04", "CHR05", "CHR09", 
                                                       "CHR12", "CHR13", "CHR19")), ])

herbqtl$QTLX <- c(27.5, 39, 50.5, 62, # 2 
                  90, 101.5, 113, 90, 90, 101.5, 113, # 4
                  141, 152.5, 141, 152.5, 164, 175.5, # 5
                  203.5, 203.5, 215, 226.5, 238, 249.5, 261, 272.5, # 9 
                  300.5, 312, 323.5, 335, 346.5, 358, 369.5, # 12
                  397.5, 409, 397.5, # 13
                  437, 437, 448.5) # 19 
herbqtl$TRT_YR <- word(herbqtl$TRT, 2, 3, sep = "_")

physgrad <- NULL
for (chr in unique(herbqtl$CHR)) {
  pullchr <- droplevels(herbqtl[which(herbqtl$CHR == chr), ])
  chrend <- round(herbref$END_Mb[which(herbref$CHR == chr)], 2)
  atphys <- 0
  splits <- NULL
  while (atphys != chrend) {
    if (atphys == 0) {
      atphys <- round(pullchr$QTLMIN_Mb[which.min(pullchr$QTLMIN_Mb)], 2)
      splits <- append(splits, atphys)
    } else {
      pos <- round(c(pullchr$QTLMIN_Mb, pullchr$QTLMAX_Mb, chrend), 2)
      pos <- sort(pos)
      atphys <- pos[which(pos > atphys)][1]
      splits <- append(splits, atphys)
    }
  }
  chrbks <- as.data.frame(list("CHR" = pullchr$CHR[1], 
                               "XMIN" = rep(herbref$REFX[which(herbref$CHR == chr)], length(splits)), 
                               "XMAX" = rep(herbref$REFX[which(herbref$CHR == chr)] + 5, length(splits)), 
                               "YMIN" = c(0, (splits[-length(splits)] + 0.01)), 
                               "YMAX" = c(splits)))
  physgrad <- rbind(physgrad, chrbks)
}

physgrad$DENS <- c(0, 1, 2, 3, 4, 1, 0,
                   0, 3, 1, 0, 1, 0, 1, 2, 3, 2, 1, 0,
                   0, 1, 2, 1, 0, 1, 2, 1, 2, 3, 2, 1, 0,
                   0, 1, 0, 1, 2, 4, 3, 4, 5, 4, 5, 4, 2, 1, 0,
                   0, 1, 2, 4, 5, 6, 5, 6, 4, 3, 2, 0,
                   0, 2, 1, 0, 1, 0,
                   0, 1, 0, 1, 2, 1, 0)

denscale <- as.data.frame(list("CHR" = rep("SCALE", 7), 
                               "XMIN" = c(250, (250 + 27.7), (250 + (27.7 * 2)), (250 + (27.7 * 3)), 
                                          (250 + (27.7 * 4)), (250 + (27.7 * 5)), (250 + (27.7 * 6))), 
                               "XMAX" = c((250 + 27.7), (250 + (27.7 * 2)), (250 + (27.7 * 3)),
                                          (250 + (27.7 * 4)), (250 + (27.7 * 5)), (250 + (27.7 * 6)), 444), 
                               "YMIN" = rep(19.2, 7), "YMAX" = rep(19.8, 7), 
                               "DENS" = c(0, 1, 2, 3, 4, 5, 6))) 

physgrad <- rbind(physgrad, denscale)

physgrad$COLTEST <- as.factor(ifelse(physgrad$DENS == 0, NA, physgrad$DENS))

ggplot( ) + 
  # setting plot dimensions # 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 455), 
                     position = "top", breaks = (herbref$REFX + 2.25), 
                     labels = herbref$CHR) +
  scale_y_reverse(expand = c(0.01, 0), limits = c(20, -0.01),
                  breaks = seq.int(0, 20, by = 1),
                  labels = seq.int(0, 20, by = 1)) +
  # qtl density heatmap # 
  geom_rect(aes(xmin = physgrad$XMIN, xmax = physgrad$XMAX,
                ymin = physgrad$YMIN, ymax = physgrad$YMAX,
                linewidth = 1, alpha = 0.55, fill = physgrad$COLTEST)) +
  scale_fill_viridis(discrete = T, begin = .4) +
  # physical Mb scale for chromosomes # 
  geom_segment(aes(x = c(herbref$REFX, herbref$REFX + 5), xend =  c(herbref$REFX, herbref$REFX + 5),
                   y = c(herbref$START, herbref$START), yend = c(herbref$END_Mb, herbref$END_Mb)),
               linewidth = 1) + 
  geom_segment(aes(x = mbticks$REFX, xend = (mbticks$REFX + mbticks$XEND),
                   y = mbticks$MBTS, yend = mbticks$MBTS),
               linewidth = 1, alpha = 0.5) + 
  # plot qtl: sex inner border with species fill # 
  geom_rect(aes(xmin = (herbqtl$QTLX - 3.75), xmax = (herbqtl$QTLX + 3.75),
                ymin = (herbqtl$QTLMIN_Mb), ymax = (herbqtl$QTLMAX_Mb)),
            linewidth = 1, alpha = 0.33, 
            fill = ifelse(herbqtl$PAR == "04BN051", "#006cf3",
                          ifelse(herbqtl$PAR %in% c("P63", "P295", "P294"), "#00b036",
                                 ifelse(herbqtl$PAR == "P336", "#29d3b9",
                                        ifelse(herbqtl$PAR == "04FF016", "#a39e91",
                                               ifelse(herbqtl$PAR %in% c("Jorr", "07MBG5027"), "#ff2b36",
                                                      ifelse(herbqtl$PAR %in% c("94006", "94001"), "#c764d9", "white")))))),
            color = ifelse(herbqtl$PAR %in% c("94001", "P63", "Jorr", "04FF016", "04BN051"), "cyan4",
                           ifelse(herbqtl$PAR %in% c("94006", "P294", "P295", "07MBG5027", "P336"), "deeppink3", "black"))) +
  # plot qtl: family outer border # 
  geom_rect(aes(xmin = (herbqtl$QTLX - 5), xmax = (herbqtl$QTLX + 5),
                        ymin = herbqtl$QTLMIN_Mb, ymax = herbqtl$QTLMAX_Mb),
                    linewidth = 1,  fill = NA,
                    color = ifelse(herbqtl$FAM == "358", "#006cf3", 
                                   ifelse(herbqtl$FAM %in% c("400", "440", "443"), "#00b036",
                                          ifelse(herbqtl$FAM == "426", "#29d3b9",
                                                 ifelse(herbqtl$FAM == "438", "#a39e91",
                                                        ifelse(herbqtl$FAM %in% c("407", "421"), "#ff2b36",
                                                               ifelse(herbqtl$FAM == "CON", "#c764d9", "white"))))))) +
  # # plot peak marker #  
  geom_segment(aes(x = (herbqtl$QTLX - 2.5), xend = (herbqtl$QTLX + 2.5),
                   y = (herbqtl$QTLPEAK_Mb), yend = (herbqtl$QTLPEAK_Mb)),
               linewidth = 1.5, col = "black") +
  # creating legend: text labels #  
  geom_text(aes(x = herbqtl$QTLX[-c(5, 18, 36)] - 0.5, 
                y = herbqtl$QTLMAX_Mb[-c(5, 18, 36)] + 0.2,
                label = herbqtl$TRT_YR[-c(5, 18, 36)]), 
            size = 3.75, fontface = "bold", angle = 90, hjust = "right") + 
  # creating legend: qtl denstiy scale # 
  geom_rect(aes(xmin = 250, xmax = 444, ymin = 19.2, ymax = 19.8), linewidth = 0.5, color = "black", fill = NA) + 
  # creating legend: parent species boxes #
  geom_rect(aes(xmin = c(250, 250, (250 + 64.7), (250 + 64.7), (250 + (64.7 * 2)), (250 + (64.7 * 2))),
                xmax = c(271.6, 271.6, (271.6 + 64.7), (271.6 + 64.7), (271.6 + (64.7 * 2)), (271.6 + (64.7 * 2))),
                ymin = rep(c(18, 18.6), 3),
                ymax = rep(c(18.4, 19), 3)),
            linewidth = 0.5, color = "black", alpha = 0.33,
            fill = c("#29d3b9", "#00b036", "#a39e91", "#ff2b36", "#c764d9", "#006cf3")) +
  # creating legend: family boxes #
  geom_rect(aes(xmin = c(250, 250, 250, (250 + 64.7), (250 + 64.7), (250 + 64.7), (250 + (64.7 * 2)), (250 + (64.7 * 2)), (250 + (64.7 * 2))),
                xmax = c(271.6, 271.6, 271.6, (271.6 + 64.7), (271.6 + 64.7), (271.6 + 64.7), (271.6 + (64.7 * 2)), (271.6 + (64.7 * 2)), (271.6 + (64.7 * 2))),
                ymin = rep(c(16.2, 16.8, 17.4), 3),
                ymax = rep(c(16.6, 17.2, 17.8), 3)),
            linewidth = 1.75, fill = NA,
            color = c("#a39e91", "#00b036", "#ff2b36", "#006cf3",
                      "#c764d9",
                      "#29d3b9", "#00b036", "#00b036", "#ff2b36")) +
  # creating legend: sex inner border #
  geom_rect(aes(xmin = c((250 + 1.5), (250 + 64.7 + 1.5), 250, (250 + 64.7)),
                xmax = c((271.6 - 1.5), (271.6 + 64.7 - 1.5), 271.6, (271.6 + 64.7)),
                ymin = c(15.6, 15.6, 15.6, 15.6),
                ymax = c(16, 16, 16, 16)),
            linewidth = 1.75, fill = NA,
            color = c("deeppink3", "cyan4", "black", "black")) +
  # creating legend: text labels #
  geom_text(aes(x = c(191, (271.6 + 3), (271.6 + 64.7 + 3),
                      (271.6 + 3), (271.6 + 64.7 + 3), (271.6 + (64.7 * 2) + 3),
                      191, (271.6 + 3), (271.6 + 64.7 + 3), (271.6 + (64.7 * 2) + 3),
                      (271.6 + 3), (271.6 + 64.7 + 3), (271.6 + (64.7 * 2) + 3),
                      (271.6 + 3), (271.6 + 64.7 + 3), (271.6 + (64.7 * 2) + 3),
                      191,
                      (271.6 + 3), (271.6 + 64.7 + 3), (271.6 + (64.7 * 2) + 3),
                      191, 263.9, 291.6, 319.3, 347, 374.7, 402.4, 430.1),
                y = c(15.75, 15.75, 15.75,
                      16.35, 16.35, 16.35,
                      16.95, 16.95, 16.95, 16.95,
                      17.55, 17.55, 17.55,
                      18.15, 18.15, 18.15,
                      18.45,
                      18.75, 18.75, 18.75,
                      19.5, 19.5, 19.5, 19.5, 19.5, 19.5, 19.5, 19.5),
                label = c("Parent Sex", "Female", "Male",
                          "S.p x S.k", "S.p x S.u", "S.s(4) x S.p",
                          "Family", "S.p x S.s", "S.p Cons", "S.s(5) x S.p",
                          "S.p x S.v", "S.i x S.p", "S.v x S.p",
                          "S. inte.", "S. kori.", "S. purp.",
                          "Par. Species", "S. such.", "S. vimen", "S. uden",
                          "QTL Density", "0", "1", "2", "3", "4", "5", "6")),
            size = 4, fontface = "bold", hjust = "left") +
  labs(x = "Chromosomes with Overlapping Herbivory QTL" , y = "Physical Distance in 94006 v5.1 Reference Genome (Mb)") + 
  THEME

rm(chrbks, mbticks, physgrad, pullchr, herbqtl, herbref, ticks, atphys, chr, chrend, pos, rc, splits)

###### LEAF ARCHITECTURE ######

archqtl <- droplevels(qtl[which(qtl$GRP == "LEAFARCH"), ])
archqtl <- droplevels(archqtl[which(archqtl$CHR %in% c("CHR03", "CHR05", "CHR11", 
                                                       "CHR12", "CHR13", 
                                                       "CHR17", "CHR19")), ])

archref <- ref
archref <- droplevels(archref[which(archref$CHR %in% archqtl$CHR), ])
archref <- droplevels(archref[which(archref$CHR %in% c("CHR03", "CHR05", "CHR11", 
                                                       "CHR12", "CHR13", 
                                                       "CHR17", "CHR19")), ])

archref$REFX <- c(15, 54.5, 105.5, 214, 288, 339, 378.5)

mbticks <- NULL
for (rc in 1:nrow(archref)) {
  ticks <- as.data.frame(list("REFX" = archref$REFX[rc],
                              "MBTS" = seq.int(0, round_any(archref$END_Mb[rc], 0.5, f = floor), by = 0.5)))
  ticks$XEND <- ifelse(grepl("[.]5", as.character(ticks$MBTS)), 2.5, 5)
  mbticks <- rbind(mbticks, ticks)
}

archqtl <- archqtl[order(archqtl$CHR, archqtl$QTLMIN_Mb), ]
archqtl$QTLX <- c(27.5, 39, 27.5, # 3 
                  67, 67, 78.5, 90, # 5
                  118, 118, 129.5, 141, 152.5, 164, 175.5, 187, 198.5, 118, # 11
                  226.5, 238, 249.5, 261, 272.5, 226.5, 226.5, 238, 249.5, 261, 226.5, 238, 226.5, 238, # 12
                  300.5, 312, 323.5, 300.5, 300.5, 300.5, # 13
                  351.5, 363, 351.5, # 17
                  391, 391, 402.5, 391) # 19
archqtl$TRT_YR <- word(archqtl$TRT, 2, 3, sep = "_")


physgrad <- NULL
for (chr in unique(archqtl$CHR)) {
  pullchr <- droplevels(archqtl[which(archqtl$CHR == chr), ])
  chrend <- round(archref$END_Mb[which(archref$CHR == chr)], 2)
  atphys <- 0
  splits <- NULL
  while (atphys != chrend) {
    if (atphys == 0) {
      atphys <- round(pullchr$QTLMIN_Mb[which.min(pullchr$QTLMIN_Mb)], 2)
      splits <- append(splits, atphys)
    } else {
      pos <- round(c(pullchr$QTLMIN_Mb, pullchr$QTLMAX_Mb, chrend), 2)
      pos <- sort(pos)
      atphys <- pos[which(pos > atphys)][1]
      splits <- append(splits, atphys)
    }
  }
  chrbks <- as.data.frame(list("CHR" = pullchr$CHR[1], 
                               "XMIN" = rep(archref$REFX[which(archref$CHR == chr)], length(splits)), 
                               "XMAX" = rep(archref$REFX[which(archref$CHR == chr)] + 5, length(splits)), 
                               "YMIN" = c(0, (splits[-length(splits)] + 0.01)), 
                               "YMAX" = c(splits)))
  physgrad <- rbind(physgrad, chrbks)
}

physgrad$DENS <- c(0, 1, 2, 1, 0, 1, 0,
                   0, 1, 0, 2, 1, 0, 
                   0, 1, 0, 1, 5, 5, 2, 0, 1, 0,
                   0, 1, 4, 5, 2, 0, 1, 0, 3, 3, 0, 2, 0, 2, 0,
                   0, 3, 1, 0, 1, 0, 1, 0, 1, 0, 
                   0, 1, 2, 0, 1, 0,
                   0, 1, 0, 1, 2, 1, 0, 1, 0)

denscale <- as.data.frame(list("CHR" = rep("SCALE", 6), 
                               "XMIN" = c(200, (200 + 33.3), (200 + (33.3 * 2)), (200 + (33.3 * 3)), 
                                          (200 + (33.3 * 4)), (200 + (33.3 * 5))), 
                               "XMAX" = c((200 + 33.3), (200 + (33.3 * 2)), (200 + (33.3 * 3)),
                                          (200 + (33.3 * 4)), (200 + (33.3 * 5)), 400), 
                               "YMIN" = rep(19.2, 6), "YMAX" = rep(19.8, 6), 
                               "DENS" = c(0, 1, 2, 3, 4, 5))) 

physgrad <- rbind(physgrad, denscale)

physgrad$COLTEST <- as.factor(ifelse(physgrad$DENS == 0, NA, physgrad$DENS))

ggplot( ) + 
  # setting plotting dimensions # 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 410), 
                     position = "top", breaks = (archref$REFX + 2.25), 
                     labels = archref$CHR) +
  scale_y_reverse(expand = c(0.01, 0), limits = c(20, -0.01),
                  breaks = seq.int(0, 20, by = 1),
                  labels = seq.int(0, 20, by = 1)) +
  # qtl density heatmap #
  geom_rect(aes(xmin = physgrad$XMIN, xmax = physgrad$XMAX,
                ymin = physgrad$YMIN, ymax = physgrad$YMAX,
                linewidth = 1, alpha = 0.55, fill = physgrad$COLTEST)) +
  scale_fill_viridis(discrete = T, begin = .4) +
  # physical Mb scale for chromosomes #
  geom_segment(aes(x = c(archref$REFX, archref$REFX + 5), xend =  c(archref$REFX, archref$REFX + 5),
                   y = c(archref$START, archref$START), yend = c(archref$END_Mb, archref$END_Mb)),
               linewidth = 1) + 
  geom_segment(aes(x = mbticks$REFX, xend = (mbticks$REFX + mbticks$XEND),
                   y = mbticks$MBTS, yend = mbticks$MBTS),
               linewidth = 1, alpha = 0.5) + 
  # plot qtl: sex inner border with species fill # 
  geom_rect(aes(xmin = (archqtl$QTLX - 3.75), xmax = (archqtl$QTLX + 3.75),
                ymin = (archqtl$QTLMIN_Mb), ymax = (archqtl$QTLMAX_Mb)),
            linewidth = 1, alpha = 0.33, 
            fill = ifelse(archqtl$PAR == "04BN051", "#006cf3",
                          ifelse(archqtl$PAR %in% c("P63", "P295", "P294"), "#00b036",
                                 ifelse(archqtl$PAR == "P336", "#29d3b9",
                                        ifelse(archqtl$PAR == "04FF016", "#a39e91",
                                               ifelse(archqtl$PAR %in% c("Jorr", "07MBG5027"), "#ff2b36",
                                                      ifelse(archqtl$PAR %in% c("94006", "94001"), "#c764d9", "white")))))),
            color = ifelse(archqtl$PAR %in% c("94001", "P63", "Jorr", "04FF016", "04BN051"), "cyan4",
                           ifelse(archqtl$PAR %in% c("94006", "P294", "P295", "07MBG5027", "P336"), "deeppink3", "black"))) +
  # plot qtl: family outer border #
  geom_rect(aes(xmin = (archqtl$QTLX - 5), xmax = (archqtl$QTLX + 5),
                ymin = archqtl$QTLMIN_Mb, ymax = archqtl$QTLMAX_Mb),
            linewidth = 1, fill = NA,
            color = ifelse(archqtl$FAM == "358", "#006cf3", 
                           ifelse(archqtl$FAM %in% c("400", "440", "443"), "#00b036",
                                  ifelse(archqtl$FAM == "426", "#29d3b9",
                                         ifelse(archqtl$FAM == "438", "#a39e91",
                                                ifelse(archqtl$FAM %in% c("407", "421"), "#ff2b36",
                                                       ifelse(archqtl$FAM == "CON", "#c764d9", "white"))))))) +
  # plot peak marker # 
  geom_segment(aes(x = (archqtl$QTLX - 2.5), xend = (archqtl$QTLX + 2.5),
                   y = (archqtl$QTLPEAK_Mb), yend = (archqtl$QTLPEAK_Mb)),
               linewidth = 1.5, col = "black") +
  # qtl text labels #
  geom_text(aes(x = archqtl$QTLX[-c(4, 23:31, 32, 35, 36, 37, 38)] - 0.5,
                y = archqtl$QTLMAX_Mb[-c(4, 23:31, 32, 35, 36, 37, 38)] + 0.2,
                label = archqtl$TRT_YR[-c(4, 23:31, 32, 35, 36, 37, 38)]),
            size = 3.75, fontface = "bold", angle = 90, hjust = "right") +
  # creating legend: qtl denstiy scale #
  geom_rect(aes(xmin = 200, xmax = 400, ymin = 19.2, ymax = 19.8), linewidth = 0.5, color = "black", fill = NA) +
  # creating legend: parent species boxes #
  geom_rect(aes(xmin = c(200, 200, (200 + 66.7), (200 + 66.7), (200 + (66.7 * 2)), (200 + (66.7 * 2))),
                xmax = c(222.2, 222.2, (222.2 + 66.7), (222.2 + 66.7), (222.2 + (66.7 * 2)), (222.2 + (66.7 * 2))),
                ymin = rep(c(18, 18.6), 3),
                ymax = rep(c(18.4, 19), 3)),
            linewidth = 0.5, color = "black", alpha = 0.33,
            fill = c("#29d3b9", "#00b036", "#a39e91", "#ff2b36", "#c764d9", "#006cf3")) +
  # creating legend: family boxes #
  geom_rect(aes(xmin = c(200, 200, 200, 
                         (200 + 66.7), (200 + 66.7), (200 + 66.7), 
                         (200 + (66.7 * 2)), (200 + (66.7 * 2)), (200 + (66.7 * 2))),
                xmax = c(222.2, 222.2, 222.2, 
                         (222.2 + 66.7), (222.2 + 66.7), (222.2 + 66.7), 
                         (222.2 + (66.7 * 2)), (222.2 + (66.7 * 2)), (222.2 + (66.7 * 2))),
                ymin = rep(c(16.2, 16.8, 17.4), 3),
                ymax = rep(c(16.6, 17.2, 17.8), 3)),
            linewidth = 1.75, fill = NA,
            color = c("#a39e91", "#00b036", "#ff2b36", "#006cf3",
                      "#c764d9",
                      "#29d3b9", "#00b036", "#00b036", "#ff2b36")) +
  # creating legend: sex inner border #
  geom_rect(aes(xmin = c((200 + 1.5), (200 + 66.7 + 1.5), 200, (200 + 66.7)),
                xmax = c((222.2 - 1.5), (222.2 + 66.7 - 1.5), 222.2, (222.2 + 66.7)),
                ymin = c(15.6, 15.6, 15.6, 15.6),
                ymax = c(16, 16, 16, 16)),
            linewidth = 1.75, fill = NA,
            color = c("deeppink3", "cyan4", "black", "black")) +
  # creating legend: text labels #
  geom_text(aes(x = c(141, (222.2 + 3), (222.2 + 66.7 + 3),
                      (222.2 + 3), (222.2 + 66.7 + 3), (222.2 + (66.7 * 2) + 3),
                      141, (222.2 + 3), (222.2 + 66.7 + 3), (222.2 + (66.7 * 2) + 3),
                      (222.2 + 3), (222.2 + 66.7 + 3), (222.2 + (66.7 * 2) + 3),
                      (222.2 + 3), (222.2 + 66.7 + 3), (222.2 + (66.7 * 2) + 3),
                      141,
                      (222.2 + 3), (222.2 + 66.7 + 3), (222.2 + (66.7 * 2) + 3),
                      141, 216.7, 250, 283.3, 316.6, 349.9, 383.2),
                y = c(15.75, 15.75, 15.75,
                      16.35, 16.35, 16.35,
                      16.95, 16.95, 16.95, 16.95,
                      17.55, 17.55, 17.55,
                      18.15, 18.15, 18.15,
                      18.45,
                      18.75, 18.75, 18.75,
                      19.5, 19.5, 19.5, 19.5, 19.5, 19.5, 19.5),
                label = c("Parent Sex", "Female", "Male",
                          "S.p x S.k", "S.p x S.u", "S.s(4) x S.p",
                          "Family", "S.p x S.s", "S.p Cons", "S.s(5) x S.p",
                          "S.p x S.v", "S.i x S.p", "S.v x S.p",
                          "S. inte.", "S. kori.", "S. purp.",
                          "Par. Species", "S. such.", "S. vimen", "S. uden",
                          "QTL Density", "0", "1", "2", "3", "4", "5")),
            size = 4, fontface = "bold", hjust = "left") +
  labs(x = "Chromosomes with Overlapping Leaf Architecture QTL" , y = "Physical Distance in 94006 v5.1 Reference Genome (Mb)") + 
  THEME

rm(chrbks, mbticks, physgrad, pullchr, archqtl, archref, ticks, atphys, chr, chrend, pos, rc, splits)

###### YIELD COMP ######

yldqtl <- droplevels(qtl[which(qtl$GRP == "YLDCOMP"), ])
yldqtl <- droplevels(yldqtl[which(yldqtl$CHR %in% c("CHR04", "CHR09", "CHR11", 
                                                    "CHR12", "CHR14", 
                                                    "CHR17", "CHR19")), ])

yldref <- ref
yldref <- droplevels(yldref[which(yldref$CHR %in% yldqtl$CHR), ])
yldref <- droplevels(yldref[which(yldref$CHR %in% c("CHR04", "CHR09", "CHR11", 
                                                    "CHR12", "CHR14", 
                                                    "CHR17", "CHR19")), ])

yldref$REFX <- c(15, 77.5, 117, 179.5, 242, 281.5, 321)

yldqtl$TRT_YR <- word(yldqtl$TRT, 2, 3, sep = "_")

mbticks <- NULL
for (rc in 1:nrow(yldref)) {
  ticks <- as.data.frame(list("REFX" = yldref$REFX[rc],
                              "MBTS" = seq.int(0, round_any(yldref$END_Mb[rc], 0.5, f = floor), by = 0.5)))
  ticks$XEND <- ifelse(grepl("[.]5", as.character(ticks$MBTS)), 2.5, 5)
  mbticks <- rbind(mbticks, ticks)
}

yldqtl <- yldqtl[order(yldqtl$CHR, yldqtl$QTLMIN_Mb), ]
yldqtl$QTLX <- c(27.5, 39, 50.5, 62, 27.5, 39, 50.5, 27.5, # 4 
                 90, 90, 101.5, 90, # 9
                 129.5, 141, 152.5, 164, # 11
                 192, 203.5, 203.5, 215, 215, 192, 203.5, 215, 226.5, 192, # 12
                 254.5, 266, 254.5, 254.5, # 14
                 294, 305.5, 294, # 17
                 333.5, 345, 356.5, 368, 379.5, 391, 402.5, 414, 425.5, 437, 448.5, 460, 471.5, 483, 494.5, 506, 517.5, 529, 540.5, 333.5) # 19

physgrad <- NULL
for (chr in unique(yldqtl$CHR)) {
  pullchr <- droplevels(yldqtl[which(yldqtl$CHR == chr), ])
  chrend <- round(yldref$END_Mb[which(yldref$CHR == chr)], 2)
  atphys <- 0
  splits <- NULL
  while (atphys != chrend) {
    if (atphys == 0) {
      atphys <- round(pullchr$QTLMIN_Mb[which.min(pullchr$QTLMIN_Mb)], 2)
      splits <- append(splits, atphys)
    } else {
      pos <- round(c(pullchr$QTLMIN_Mb, pullchr$QTLMAX_Mb, chrend), 2)
      pos <- sort(pos)
      atphys <- pos[which(pos > atphys)][1]
      splits <- append(splits, atphys)
    }
  }
  chrbks <- as.data.frame(list("CHR" = pullchr$CHR[1], 
                               "XMIN" = rep(yldref$REFX[which(yldref$CHR == chr)], length(splits)), 
                               "XMAX" = rep(yldref$REFX[which(yldref$CHR == chr)] + 5, length(splits)), 
                               "YMIN" = c(0, (splits[-length(splits)] + 0.01)), 
                               "YMAX" = c(splits)))
  physgrad <- rbind(physgrad, chrbks)
}

physgrad$DENS <- c(0, 1, 3, 4, 2, 0, 3, 0, 1, 0, 
                   0, 1, 0, 1, 2, 1, 0, 1, 0, 
                   0, 2, 1, 3, 2, 1, 0,
                   0, 1, 2, 1, 2, 1, 2, 1, 2, 0, 1, 2, 4, 0, 1, 0,
                   0, 1, 2, 0, 1, 0, 1, 0, 
                   0, 1, 2, 0, 1, 0,
                   0, 3, 2, 5, 10, 11, 12, 13, 14, 13, 10, 7, 10, 6, 6, 3, 2, 1, 1, 0, 1, 0)

denscale <- as.data.frame(list("CHR" = rep("SCALE", 15), 
                               "XMIN" = c(409, (409 + 26.2), (409 + (26.2 * 2)), (409 + (26.2 * 3)), (409 + (26.2 * 4)),
                                          409, (409 + 26.2), (409 + (26.2 * 2)), (409 + (26.2 * 3)), (409 + (26.2 * 4)),
                                          409, (409 + 26.2), (409 + (26.2 * 2)), (409 + (26.2 * 3)), (409 + (26.2 * 4))),
                               "XMAX" = c((409 + 26.2), (409 + (26.2 * 2)), (409 + (26.2 * 3)), (409 + (26.2 * 4)), 540, 
                                          (409 + 26.2), (409 + (26.2 * 2)), (409 + (26.2 * 3)), (409 + (26.2 * 4)), 540, 
                                          (409 + 26.2), (409 + (26.2 * 2)), (409 + (26.2 * 3)), (409 + (26.2 * 4)), 540),
                               "YMIN" = c(rep(15.2, 5), rep(14.4, 5), rep(13.6, 5)),
                               "YMAX" = c(rep(15.8, 5), rep(15, 5), rep(14.2, 5)), 
                               "DENS" = c(0, 1, 2, 3, 4, 
                                          5, 6, 7, 8, 9, 
                                          10, 11, 12, 13, 14))) 

physgrad <- rbind(physgrad, denscale)

physgrad$COLTEST <- as.factor(ifelse(physgrad$DENS == 0, NA, physgrad$DENS))

ggplot( ) + 
  # setting plotting dimensions # 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 550), 
                     position = "top", breaks = (yldref$REFX + 2.25), 
                     labels = yldref$CHR) +
  scale_y_reverse(expand = c(0.01, 0), limits = c(16, -0.01),
                  breaks = seq.int(0, 16, by = 1),
                  labels = seq.int(0, 16, by = 1)) +
  # qtl density heatmap #
  geom_rect(aes(xmin = physgrad$XMIN, xmax = physgrad$XMAX,
                ymin = physgrad$YMIN, ymax = physgrad$YMAX,
                linewidth = 1, alpha = 0.55, fill = physgrad$COLTEST)) +
  scale_fill_viridis(discrete = T, begin = .4) +
  # physical Mb scale for chromosomes #
  geom_segment(aes(x = c(yldref$REFX, yldref$REFX + 5), xend =  c(yldref$REFX, yldref$REFX + 5),
                   y = c(yldref$START, yldref$START), yend = c(yldref$END_Mb, yldref$END_Mb)),
               linewidth = 1) + 
  geom_segment(aes(x = mbticks$REFX, xend = (mbticks$REFX + mbticks$XEND),
                   y = mbticks$MBTS, yend = mbticks$MBTS),
               linewidth = 1, alpha = 0.5) + 
  # plot qtl: sex inner border with species fill # 
  geom_rect(aes(xmin = (yldqtl$QTLX - 3.75), xmax = (yldqtl$QTLX + 3.75),
                ymin = (yldqtl$QTLMIN_Mb), ymax = (yldqtl$QTLMAX_Mb)),
            linewidth = 1, alpha = 0.33, 
            fill = ifelse(yldqtl$PAR == "04BN051", "#006cf3",
                          ifelse(yldqtl$PAR %in% c("P63", "P295", "P294"), "#00b036",
                                 ifelse(yldqtl$PAR == "P336", "#29d3b9",
                                        ifelse(yldqtl$PAR == "04FF016", "#a39e91",
                                               ifelse(yldqtl$PAR %in% c("Jorr", "07MBG5027"), "#ff2b36",
                                                      ifelse(yldqtl$PAR %in% c("94006", "94001"), "#c764d9", "white")))))),
            color = ifelse(yldqtl$PAR %in% c("94001", "P63", "Jorr", "04FF016", "04BN051"), "cyan4",
                           ifelse(yldqtl$PAR %in% c("94006", "P294", "P295", "07MBG5027", "P336"), "deeppink3", "black"))) +
  # plot qtl: family outer border #
  geom_rect(aes(xmin = (yldqtl$QTLX - 5), xmax = (yldqtl$QTLX + 5),
                ymin = yldqtl$QTLMIN_Mb, ymax = yldqtl$QTLMAX_Mb),
            linewidth = 1, fill = NA, 
            color = ifelse(yldqtl$FAM == "358", "#006cf3", 
                           ifelse(yldqtl$FAM %in% c("400", "440", "443"), "#00b036",
                                  ifelse(yldqtl$FAM == "426", "#29d3b9",
                                         ifelse(yldqtl$FAM == "438", "#a39e91",
                                                ifelse(yldqtl$FAM %in% c("407", "421"), "#ff2b36",
                                                       ifelse(yldqtl$FAM == "CON", "#c764d9", "white"))))))) +
  # plot peak marker # 
  geom_segment(aes(x = (yldqtl$QTLX - 2.5), xend = (yldqtl$QTLX + 2.5),
                   y = (yldqtl$QTLPEAK_Mb), yend = (yldqtl$QTLPEAK_Mb)),
               linewidth = 1.5, col = "black") +
  # qtl text labels #
  geom_text(aes(x = yldqtl$QTLX[-c(2, 3, 9, 10, 18:22, 29, 31)] - 0.5,
                y = yldqtl$QTLMAX_Mb[-c(2, 3, 9, 10, 18:22, 29, 31)] + 0.2,
                label = yldqtl$TRT_YR[-c(2, 3, 9, 10, 18:22, 29, 31)]),
            size = 3.75, fontface = "bold", angle = 90, hjust = "right") +
  # creating legend: qtl denstiy scale #
  geom_rect(aes(xmin = c(409, 409, 409), xmax = c(540, 540, 540), 
                ymin = c(15.2, 14.4, 13.6), ymax = c(15.8, 15, 14.2)), 
            linewidth = 0.5, color = "black", fill = NA) +
  # creating legend: parent species boxes #
  geom_rect(aes(xmin = c(409, 409, 409, 
                         (409 + 65.5), (409 + 65.5), (409 + 65.5)),
                xmax = c(430.8, 430.8, 430.8, 
                         (430.8 + 65.5), (430.8 + 65.5), (430.8 + 65.5)),
                ymin = rep(c(11.2, 12, 12.8), 2),
                ymax = rep(c(11.8, 12.6, 13.4), 2)),
            linewidth = 0.5, color = "black", alpha = 0.33,
            fill = c("#29d3b9", "#00b036", "#a39e91", "#ff2b36", "#c764d9", "#006cf3")) +
  # creating legend: family boxes #
  geom_rect(aes(xmin = c(409, 409, 409, 409, 409,
                         (409 + 65.5), (409 + 65.5), (409 + 65.5), (409 + 65.5)),
                xmax = c(430.8, 430.8, 430.8, 430.8, 430.8,
                         (430.8 + 65.5), (430.8 + 65.5), (430.8 + 65.5), (430.8 + 65.5)),
                ymin = c(7.2, 8, 8.8, 9.6, 10.4, 7.2, 8, 8.8, 9.6),
                ymax = c(7.8, 8.6, 9.4, 10.2, 11, 7.8, 8.6, 9.4, 10.2)),
            linewidth = 1.75, fill = NA,
            color = c("#a39e91", "#00b036", "#ff2b36", "#006cf3",
                      "#c764d9",
                      "#29d3b9", "#00b036", "#00b036", "#ff2b36")) +
  # creating legend: sex inner border #
  geom_rect(aes(xmin = c((409 + 1.5), (409 + 65.5 + 1.5), 409, (409 + 65.5)),
                xmax = c((430.8 - 1.5), (430.8 + 65.5 - 1.5), 430.8, (430.8 + 65.5)),
                ymin = c(6.4, 6.4, 6.4, 6.4),
                ymax = c(7, 7, 7, 7)),
            linewidth = 1.75, fill = NA,
            color = c("deeppink3", "cyan4", "black", "black")) +
  # creating legend: text labels #
  geom_text(aes(x = c(345, (430.8 + 3), (430.8 + 65.5 + 3),
                      (430.8 + 3), (430.8 + 65.5 + 3),
                      (430.8 + 3), (430.8 + 65.5 + 3),
                      345, (430.8 + 3), (430.8 + 65.5 + 3), 
                      (430.8 + 3), (430.8 + 65.5 + 3),
                      (430.8 + 3), 
                      (430.8 + 3), (430.8 + 65.5 + 3),
                      345, (430.8 + 3), (430.8 + 65.5 + 3),
                      (430.8 + 3), (430.8 + 65.5 + 3),
                      417.1, 443.3, 469.5, 495.7, 521.9,
                      345, 419.6, 445.8, 472.0, 498.2, 524.1,
                      419.6, 445.8, 472.0, 498.2, 524.1),
                y = c(6.65, 6.65, 6.65, 
                      7.45, 7.45, 
                      8.25, 8.25, 
                      9.05, 9.05, 9.05,
                      9.85, 9.85,
                      10.65, 
                      11.5, 11.5,
                      12.25, 12.25, 12.25,
                      13.05, 13.05,
                      13.85, 13.85, 13.85, 13.85, 13.85, 
                      14.65, 14.65, 14.65, 14.65, 14.65, 14.65,
                      15.5, 15.5, 15.5, 15.5, 15.5),
                label = c("Parent Sex", "Female", "Male",
                          "S.p x S.k", "S.i x S.p", 
                          "S.p x S.s", "S.s(4) x S.p", 
                          "Family", "S.p x S.v",  "S.s(5) x S.p",
                          "S.p x S.u", "S.v x S.p", "S.p Cons",
                          "S. inte.", "S. vimen", 
                          "Par. Species",  "S. such.", "S. purp.",
                          "S. kori.", "S. uden",
                          "10", "11", "12", "13", "14",
                          "QTL Density", "5", "6", "7", "8", "9", 
                          "0", "1", "2", "3", "4")),
            size = 3.75, fontface = "bold", hjust = "left") +
  labs(x = "Chromosomes with Overlapping Yield Component QTL" , y = "Physical Distance in 94006 v5.1 Reference Genome (Mb)") + 
  THEME

rm(chrbks, mbticks, physgrad, pullchr, yldqtl, yldref, ticks, atphys, chr, chrend, pos, rc, splits)

###### FINDING HIGH QUALITY OVERLAP REGIONS ###### 

### WITHIN ### 

qtl <- read.csv("QTLCurated.csv", stringsAsFactors = T)
qtl <- qtl[ , c(1:5, 9:10, 12:14)]
qtl[ , 8:10] <- (qtl[ , 8:10] * 1000000)

ref <- read.table("SpRefChrBounds.txt", colClasses = c("character", "numeric", "numeric"), sep = "\t")
colnames(ref) <- c("CHR", "START", "END")

#qtl <- droplevels(qtl[which(qtl$GRP == "LEAFRUST"), ])
#qtl <- droplevels(qtl[which(qtl$GRP == "HERB"), ])
#qtl <- droplevels(qtl[which(qtl$GRP == "LEAFARCH"), ])
qtl <- droplevels(qtl[which(qtl$GRP == "YLDCOMP"), ])

chr <- "CHR19"
pullchr <- droplevels(qtl[which(qtl$CHR == chr), ])
chrend <- ref$END[which(ref$CHR == chr)]
atphys <- 0
splits <- NULL
while (atphys != chrend) {
  if (atphys == 0) {
    atphys <- pullchr$QTLMIN_Mb[which.min(pullchr$QTLMIN_Mb)]
    splits <- append(splits, atphys)
  } else {
    pos <- c(pullchr$QTLMIN_Mb, pullchr$QTLMAX_Mb, chrend)
    pos <- sort(pos)
    atphys <- pos[which(pos > atphys)][1]
    splits <- append(splits, atphys)
  }
}

chrbks <- as.data.frame(list("CHR" = pullchr$CHR[1], 
                             "YMIN" = c(0, (splits[-length(splits)])), 
                             "YMAX" = splits))

chrbks$BRKDENS <- unlist(lapply(1:nrow(chrbks), function(x) {
  sum((pullchr$QTLMIN_Mb == chrbks$YMIN[x] & pullchr$QTLMAX_Mb > chrbks$YMAX[x]) | 
        (pullchr$QTLMIN_Mb < chrbks$YMIN[x] & pullchr$QTLMAX_Mb == chrbks$YMAX[x]) | 
        (pullchr$QTLMIN_Mb == chrbks$YMIN[x] & pullchr$QTLMAX_Mb == chrbks$YMAX[x]) |
        (pullchr$QTLMIN_Mb < chrbks$YMIN[x] & pullchr$QTLMAX_Mb > chrbks$YMAX[x]))
  
}))

### BETWEEN ### 

qtl <- read.csv("QTLCurated.csv", stringsAsFactors = T)
qtl <- qtl[ , c(1:5, 12:14)]
qtl[ , 6:8] <- (qtl[ , 6:8] * 1000000)

ref <- read.table("SpRefChrBounds.txt", colClasses = c("character", "numeric", "numeric"), sep = "\t")
colnames(ref) <- c("CHR", "START", "END")

chr <- "CHR19"
pullchr <- droplevels(qtl[which(qtl$CHR == chr), ])
chrend <- ref$END[which(ref$CHR == chr)]
atphys <- 0
splits <- NULL
while (atphys != chrend) {
  if (atphys == 0) {
    atphys <- pullchr$QTLMIN_Mb[which.min(pullchr$QTLMIN_Mb)]
    splits <- append(splits, atphys)
  } else {
    pos <- c(pullchr$QTLMIN_Mb, pullchr$QTLMAX_Mb, chrend)
    pos <- sort(pos)
    atphys <- pos[which(pos > atphys)][1]
    splits <- append(splits, atphys)
  }
}

chrbks <- as.data.frame(list("CHR" = pullchr$CHR[1], 
                             "YMIN" = c(0, (splits[-length(splits)])), 
                             "YMAX" = splits))

chrbks$BRKDENS <- unlist(lapply(1:nrow(chrbks), function(x) {
  sum((pullchr$QTLMIN_Mb == chrbks$YMIN[x] & pullchr$QTLMAX_Mb > chrbks$YMAX[x]) | 
        (pullchr$QTLMIN_Mb < chrbks$YMIN[x] & pullchr$QTLMAX_Mb == chrbks$YMAX[x]) | 
        (pullchr$QTLMIN_Mb == chrbks$YMIN[x] & pullchr$QTLMAX_Mb == chrbks$YMAX[x]) |
        (pullchr$QTLMIN_Mb < chrbks$YMIN[x] & pullchr$QTLMAX_Mb > chrbks$YMAX[x]))
        
}))

###### CREATING REFERENCE GENE LIST WITH POSITIONS ###### 

# read in annotation file and rename columns 
annot <- read.delim("Spurpurea_519_v5.1.annotation_info.txt",
                    header = T, sep = "\t")[ , c(2, 10:13)]
colnames(annot) <- c("LOCUS", "GO", "BHARABI", "ARABISYM", "ARABIDEF")

# remove duplicate rows; isolate genes on chrs and 15W
annot <- annot[!duplicated(annot), ]
annot <- annot[!duplicated(annot$LOCUS), ]
annot <- droplevels(annot[which(!grepl("Sapur.T", annot$LOCUS)), ])
annot <- droplevels(annot[which(!grepl("Sapur.15Z", annot$LOCUS)), ])

# read in gene gff3; rename columns, 
genes <- read.delim("Spurpurea_519_v5.1_no_15Z.gene.gff3",
                    header = F, sep = "\t", skip = 3)[ , c(1, 3:5, 9)]
colnames(genes) <- c("CHR", "TYPE", "START", "STOP", "INFO")

# isolate genes on chromosomes
genes <- droplevels(genes[which(genes$TYPE == "gene"), ])
genes <- droplevels(genes[which(grepl("Chr", genes$CHR)), ]) # 31831

# pull gene name from info column; modify CHR column to match QTL df
genes$LOCUS <- str_replace(word(genes$INFO, 2, sep = "[;]"), "Name=", "")
genes$CHR <- str_replace(genes$CHR, "Chr", "CHR")
genes$CHR <- ifelse(genes$CHR == "CHR15W", "CHR15", genes$CHR)

# remove unneccesary columns and merge with annotation file 
genes <- genes[ , -c(2, 5)]

generef <- merge(genes, annot, by = "LOCUS", all = T)

write.table(generef, "94006GeneList.txt", quote = F, 
            sep = "\t", row.names = F)

###### INTO THE OVERLAP ###### 

ovrlp <- read.csv("HCPQTLOVERLAP.csv", stringsAsFactors = T)
generef <- read.delim("94006GeneList.txt", header = T, stringsAsFactors = T)

ovrlp <- ovrlp[ , -c(5:7)]

ovrlp$OVLPSTR <- ovrlp$OVLPSTR*1e6
ovrlp$OVLPSTP <- ovrlp$OVLPSTP*1e6
ovrlp$PEAKSTR <- ovrlp$PEAKSTR*1e6
ovrlp$PEAKSTP <- ovrlp$PEAKSTP*1e6

ovrlp$GENECT <- unlist(lapply(1:nrow(ovrlp), function(x) {
  pullchr <- droplevels(generef[which(generef$CHR == as.character(ovrlp$CHR)[x]), ])
  sum((pullchr$START <= ovrlp$OVLPSTR[x] & 
         pullchr$STOP > ovrlp$OVLPSTR[x] & pullchr$STOP < ovrlp$OVLPSTP[x]) |
      (pullchr$START > ovrlp$OVLPSTR[x] & pullchr$START < ovrlp$OVLPSTP[x] &
         pullchr$STOP > ovrlp$OVLPSTP[x]) |
      (pullchr$START > ovrlp$OVLPSTR[x] & pullchr$STOP < ovrlp$OVLPSTP[x]) | 
      (pullchr$START < ovrlp$OVLPSTR[x] & pullchr$STOP > ovrlp$OVLPSTP[x]))
}))

ovrlp$PEAKCT <- unlist(lapply(1:nrow(ovrlp), function(x) {
  pullchr <- droplevels(generef[which(generef$CHR == as.character(ovrlp$CHR)[x]), ])
  sum((pullchr$START <= ovrlp$PEAKSTR[x] & 
         pullchr$STOP > ovrlp$PEAKSTR[x] & pullchr$STOP < ovrlp$PEAKSTP[x]) |
        (pullchr$START > ovrlp$PEAKSTR[x] & pullchr$START < ovrlp$PEAKSTP[x] &
           pullchr$STOP > ovrlp$PEAKSTP[x]) |
        (pullchr$START > ovrlp$PEAKSTR[x] & pullchr$STOP < ovrlp$PEAKSTP[x]) | 
        (pullchr$START < ovrlp$PEAKSTR[x] & pullchr$STOP > ovrlp$PEAKSTP[x]))
}))

ovrlp$INT <- ifelse(ovrlp$PEAKCT == 0, "OVLP",
                    ifelse(ovrlp$GENECT == ovrlp$PEAKCT, "PKMR", 
                           ifelse(ovrlp$GENECT < ovrlp$PEAKCT, "OVLP", "PKMR")))

top5 <- lapply(1:nrow(ovrlp), function(x) {
  pullchr <- droplevels(generef[which(generef$CHR == as.character(ovrlp$CHR)[x]), ])
  pullqtl <- droplevels(pullchr[which((pullchr$START <= ovrlp$OVLPSTR[x] & 
                                         pullchr$STOP > ovrlp$OVLPSTR[x] & pullchr$STOP < ovrlp$OVLPSTP[x]) |
                                        (pullchr$START > ovrlp$OVLPSTR[x] & pullchr$START < ovrlp$OVLPSTP[x] &
                                           pullchr$STOP > ovrlp$OVLPSTP[x]) |
                                        (pullchr$START > ovrlp$OVLPSTR[x] & pullchr$STOP < ovrlp$OVLPSTP[x]) | 
                                        (pullchr$START < ovrlp$OVLPSTR[x] & pullchr$STOP > ovrlp$OVLPSTP[x])), ])
  
  defline <- as.data.frame(table(pullqtl$ARABIDEF))
  defline <- defline[order(-defline$Freq), ]
  if ("" %in% defline$Var1) { 
    unannot <- defline$Freq[which(defline$Var1 == "")]
  } else {
    unannot <- 0
  }
  
  defline <- droplevels(defline[which(defline$Var1 != ""), ])

  if (nrow(defline) >= 5) {
    c(as.character(defline$Var1[1]), defline$Freq[1], 
      as.character(defline$Var1[2]), defline$Freq[2],
      as.character(defline$Var1[3]), defline$Freq[3],
      as.character(defline$Var1[4]), defline$Freq[4],
      as.character(defline$Var1[5]), defline$Freq[5], 
      unannot) 
  } else {
    if (nrow(defline) == 4) {
      c(as.character(defline$Var1[1]), defline$Freq[1], 
        as.character(defline$Var1[2]), defline$Freq[2],
        as.character(defline$Var1[3]), defline$Freq[3],
        as.character(defline$Var1[4]), defline$Freq[4],
        NA, NA, 
        unannot)
    } else { 
      if (nrow(defline) == 3) {
        c(as.character(defline$Var1[1]), defline$Freq[1], 
          as.character(defline$Var1[2]), defline$Freq[2],
          as.character(defline$Var1[3]), defline$Freq[3],
          NA, NA, NA, NA,
          unannot)
      } else {
        if (nrow(defline) == 2) {
          c(as.character(defline$Var1[1]), defline$Freq[1], 
            as.character(defline$Var1[2]), defline$Freq[2],
            NA, NA, NA, NA, NA, NA,
            unannot)
        } else { 
          if (nrow(defline) == 1) {
            c(as.character(defline$Var1[1]), defline$Freq[1], 
              NA, NA, NA, NA, NA, NA, NA, NA,
              unannot)
            } 
        }
      }
    }
  }
})

top5 <- do.call("rbind", top5)
colnames(top5) <- c("FIRST", "FSTCT",
                    "SECOND", "SNDCT",
                    "THIRD", "TRDCT", 
                    "FOURTH", "FRHCT", 
                    "FIFTH", "FTHCT", 
                    "UNANCT")

ovrlp <- cbind(ovrlp, top5)

write.csv(ovrlp, "QTLMapping/GenesNSummary.csv", quote = F, 
            row.names = F)
