ptm <- proc.time()
require(ape)
require(phangorn)
require(phytools)
require(caper)

filename <- read.csv('Settings.csv')
filename <- filename$runName
simlist <- (list())
estlist <- (list())
folders1 <- grep("Model1", dir(), value = T)
folders2 <- grep("Model2", dir(), value = T)
folders3 <- grep("Model3", dir(), value = T)
folders4 <- grep("Model4", dir(), value = T)
folders5 <- grep("Model5", dir(), value = T)
M1FPDtest <-  vector()
M2FPDtest <-  vector()
M3FPDtest <-  vector()
M4FPDtest <-  vector()
M1FPDVvalue <- vector()
M2FPDVvalue <- vector()
M3FPDVvalue <- vector()
M4FPDVvalue <- vector() 
M1simFPDstn <- vector()
M1estFPDstn <- vector()
M1RankSimFPDstn <- vector()
M1RankEstFPDstn <- vector()
M2simFPDstn <- vector()
M2estFPDstn <- vector()
M2RankSimFPDstn <- vector()
M2RankEstFPDstn <- vector()
M3simFPDstn <- vector()
M3estFPDstn <- vector()
M3RankSimFPDstn <- vector()
M3RankEstFPDstn <- vector()
M4simFPDstn <- vector()
M4estFPDstn <- vector()
M4RankSimFPDstn <- vector()
M4RankEstFPDstn <- vector()
n <- 0

#Model1
for(i in 1: length(folders1)){
  #print(getwd())
  setwd (folders1[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders1[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  FPDSimDist <- vector()
  FPDEstDist <- vector()
  tipsample <- list()
  for(i in 1:100){
    relatedness <- vector()
    #we choose our starting tip and name it x
    x <- sample(simtree$tip.label, 1)
    #we make a vector of all the tip labels other than x
    remaindertips <- simtree$tip.label 
    remaindertips <- remaindertips[!remaindertips==x]
    for(i in 1:length(remaindertips)){
      tipsub <- c(x, (remaindertips[i]))
      PairPD <- pd.calc(simtree, tip.subset = tipsub)
      relatedness <- c(relatedness, PairPD)
    }
    tipswpd <- data.frame(I(remaindertips), relatedness)
    tipswpd <- tipswpd[order(tipswpd$relatedness),]
    tips <- tipswpd$remaindertips
    tips <- tips[1:19]
    tips <- c(x, tips)
    tipsample <- rbind(tipsample, tips)
    FPDSim <- pd.calc(simtree, tip.subset = tips)
    FPDSimStn <- FPDSim/TBLsimtree
    FPDEst <- pd.calc(esttree, tip.subset = tips)
    FPDEstStn <- FPDEst/TBLesttree
    FPDSimDist <- c(FPDSimDist, FPDSimStn)
    FPDEstDist <- c(FPDEstDist, FPDEstStn)
    print("we are the champions my friends")
  }
  print(FPDSimDist)
  print(FPDEstDist)
  save(tipsample, file = 'sampletips.Rdata')
  FPDWilcox <- wilcox.test(FPDSimDist, FPDEstDist, paired = F)
  FPDttest <- t.test(FPDSimDist, FPDEstDist, paired = F)
  SWTestSim <- shapiro.test(FPDSimDist)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(FPDEstDist)
  SWTestEstP <- SWTestEst$p.value
  p <- FPDWilcox$p.value
  M1FPDtest <-  c(M1FPDtest, p)
  V <- FPDWilcox$statistic
  M1FPDVvalue <- c(M1FPDVvalue, V) 
  RankedFPDSimDist <- rank(FPDSimDist)
  RankedFPDEstDist <- rank(FPDEstDist)
  M1simFPDstn <- c(M1simFPDstn, FPDSimDist)
  M1estFPDstn <- c(M1estFPDstn, FPDEstDist)
  M1RankSimFPDstn <- c(M1RankSimFPDstn, RankedFPDSimDist)
  M1RankEstFPDstn <- c(M1RankEstFPDstn, RankedFPDEstDist)
  n <- n+1
  print(n)
  setwd("..")
}
setwd("results/FPD/Clade/Sample20")
save(M1FPDtest, file = "M1FPDpvalueCd.Rdata")
save(M1FPDVvalue, file = "M1FPDVvalueCd.Rdata")

png("Model1_FPDCd_p.values_hist.png")
hist(M1FPDtest, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

setwd("../../../..")

#Model2
for(i in 1: length(folders2)){
  #print(getwd())
  setwd (folders2[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders2[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  FPDSimDist <- vector()
  FPDEstDist <- vector()
  tipsample <- list()
  for(i in 1:100){
    relatedness <- vector()
    
    #we choose our starting tip and name it x
    x <- sample(simtree$tip.label, 1)
    #we make a vector of all the tip labels other than x
    remaindertips <- simtree$tip.label 
    remaindertips <- remaindertips[!remaindertips==x]
    for(i in 1:length(remaindertips)){
      tipsub <- c(x, (remaindertips[i]))
      PairPD <- pd.calc(simtree, tip.subset = tipsub)
      relatedness <- c(relatedness, PairPD)
    }
    tipswpd <- data.frame(I(remaindertips), relatedness)
    tipswpd <- tipswpd[order(tipswpd$relatedness),]
    tips <- tipswpd$remaindertips
    tips <- tips[1:19]
    tips <- c(x, tips)
    tipsample <- rbind(tipsample, tips)
    FPDSim <- pd.calc(simtree, tip.subset = tips)
    FPDSimStn <- FPDSim/TBLsimtree
    FPDEst <- pd.calc(esttree, tip.subset = tips)
    FPDEstStn <- FPDEst/TBLesttree
    FPDSimDist <- c(FPDSimDist, FPDSimStn)
    FPDEstDist <- c(FPDEstDist, FPDEstStn)
    print("we are the champions my friends")
  }
  print(FPDSimDist)
  print(FPDEstDist)
  save(tipsample, file = 'sampletips.Rdata')
  FPDWilcox <- wilcox.test(FPDSimDist, FPDEstDist, paired = F)
  FPDttest <- t.test(FPDSimDist, FPDEstDist, paired = F)
  SWTestSim <- shapiro.test(FPDSimDist)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(FPDEstDist)
  SWTestEstP <- SWTestEst$p.value
  p <- FPDWilcox$p.value
  M2FPDtest <-  c(M2FPDtest, p)
  V <- FPDWilcox$statistic
  M2FPDVvalue <- c(M2FPDVvalue, V) 
  RankedFPDSimDist <- rank(FPDSimDist)
  RankedFPDEstDist <- rank(FPDEstDist)
  M2simFPDstn <- c(M2simFPDstn, FPDSimDist)
  M2estFPDstn <- c(M2estFPDstn, FPDEstDist)
  M2RankSimFPDstn <- c(M2RankSimFPDstn, RankedFPDSimDist)
  M2RankEstFPDstn <- c(M2RankEstFPDstn, RankedFPDEstDist)
  n <- n+1
  print(n)
  setwd("..")
}
setwd("results/FPD/Clade/Sample20")
save(M2FPDtest, file = "M2FPDpvalueCd.Rdata")
save(M2FPDVvalue, file = "M2FPDVvalueCd.Rdata")  

png("Model2_FPDCd_p.values_hist.png")
hist(M2FPDtest, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

setwd("../../../..")

#Model3
for(i in 1: length(folders3)){
  #print(getwd())
  setwd (folders3[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders3[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  FPDSimDist <- vector()
  FPDEstDist <- vector()
  tipsample <- list()
  for(i in 1:100){
    relatedness <- vector()
    
    #we choose our starting tip and name it x
    x <- sample(simtree$tip.label, 1)
    #we make a vector of all the tip labels other than x
    remaindertips <- simtree$tip.label 
    remaindertips <- remaindertips[!remaindertips==x]
    for(i in 1:length(remaindertips)){
      tipsub <- c(x, (remaindertips[i]))
      PairPD <- pd.calc(simtree, tip.subset = tipsub)
      relatedness <- c(relatedness, PairPD)
    }
    tipswpd <- data.frame(I(remaindertips), relatedness)
    tipswpd <- tipswpd[order(tipswpd$relatedness),]
    tips <- tipswpd$remaindertips
    tips <- tips[1:19]
    tips <- c(x, tips)
    tipsample <- rbind(tipsample, tips)
    FPDSim <- pd.calc(simtree, tip.subset = tips)
    FPDSimStn <- FPDSim/TBLsimtree
    FPDEst <- pd.calc(esttree, tip.subset = tips)
    FPDEstStn <- FPDEst/TBLesttree
    FPDSimDist <- c(FPDSimDist, FPDSimStn)
    FPDEstDist <- c(FPDEstDist, FPDEstStn)
    print("we are the champions my friends")
  }
  
  print(FPDSimDist)
  print(FPDEstDist)
  save(tipsample, file = 'sampletips.Rdata')
  FPDWilcox <- wilcox.test(FPDSimDist, FPDEstDist, paired = F)
  FPDttest <- t.test(FPDSimDist, FPDEstDist, paired = F)
  SWTestSim <- shapiro.test(FPDSimDist)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(FPDEstDist)
  SWTestEstP <- SWTestEst$p.value
  p <- FPDWilcox$p.value
  M3FPDtest <-  c(M3FPDtest, p)
  V <- FPDWilcox$statistic
  M3FPDVvalue <- c(M3FPDVvalue, V) 
  RankedFPDSimDist <- rank(FPDSimDist)
  RankedFPDEstDist <- rank(FPDEstDist)
  M3simFPDstn <- c(M3simFPDstn, FPDSimDist)
  M3estFPDstn <- c(M3estFPDstn, FPDEstDist)
  M3RankSimFPDstn <- c(M3RankSimFPDstn, RankedFPDSimDist)
  M3RankEstFPDstn <- c(M3RankEstFPDstn, RankedFPDEstDist)
  n <- n+1
  print(n)
  setwd("..")
}
setwd("results/FPD/Clade/Sample20")
save(M3FPDtest, file = "M3FPDpvalueCd.Rdata")
save(M3FPDVvalue, file = "M3FPDVvalueCd.Rdata")

png("Model3_FPDCd_p.values_hist.png")
hist(M3FPDtest, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

setwd("../../../..")

#Model4
for(i in 1: length(folders4)){
  #print(getwd())
  setwd (folders4[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders4[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  FPDSimDist <- vector()
  FPDEstDist <- vector()
  tipsample <- list()
  for(i in 1:100){
    relatedness <- vector()
    
    #we choose our starting tip and name it x
    x <- sample(simtree$tip.label, 1)
    #we make a vector of all the tip labels other than x
    remaindertips <- simtree$tip.label 
    remaindertips <- remaindertips[!remaindertips==x]
    for(i in 1:length(remaindertips)){
      tipsub <- c(x, (remaindertips[i]))
      PairPD <- pd.calc(simtree, tip.subset = tipsub)
      relatedness <- c(relatedness, PairPD)
    }
    tipswpd <- data.frame(I(remaindertips), relatedness)
    tipswpd <- tipswpd[order(tipswpd$relatedness),]
    tips <- tipswpd$remaindertips
    tips <- tips[1:19]
    tips <- c(x, tips)
    tipsample <- rbind(tipsample, tips)
    FPDSim <- pd.calc(simtree, tip.subset = tips)
    FPDSimStn <- FPDSim/TBLsimtree
    FPDEst <- pd.calc(esttree, tip.subset = tips)
    FPDEstStn <- FPDEst/TBLesttree
    FPDSimDist <- c(FPDSimDist, FPDSimStn)
    FPDEstDist <- c(FPDEstDist, FPDEstStn)
    print("we are the champions my friends")
  }
  
  print(FPDSimDist)
  print(FPDEstDist)
  save(tipsample, file = 'sampletips.Rdata')
  FPDWilcox <- wilcox.test(FPDSimDist, FPDEstDist, paired = F)
  FPDttest <- t.test(FPDSimDist, FPDEstDist, paired = F)
  SWTestSim <- shapiro.test(FPDSimDist)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(FPDEstDist)
  SWTestEstP <- SWTestEst$p.value
  p <- FPDWilcox$p.value
  M4FPDtest <-  c(M4FPDtest, p)
  V <- FPDWilcox$statistic
  M4FPDVvalue <- c(M4FPDVvalue, V) 
  RankedFPDSimDist <- rank(FPDSimDist)
  RankedFPDEstDist <- rank(FPDEstDist)
  M4simFPDstn <- c(M4simFPDstn, FPDSimDist)
  M4estFPDstn <- c(M4estFPDstn, FPDEstDist)
  M4RankSimFPDstn <- c(M4RankSimFPDstn, RankedFPDSimDist)
  M4RankEstFPDstn <- c(M4RankEstFPDstn, RankedFPDEstDist)
  n <- n+1
  print(n)
  setwd("..")
}
setwd("results/FPD/Clade/Sample20")
save(M4FPDtest, file = "M4FPDpvalueCd.Rdata")
save(M4FPDVvalue, file = "M4FPDVvalueCd.Rdata")

png("Model4_FPDCd_p.values_hist.png")
hist(M4FPDtest, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

png("Model1_FPD_scores.png")
plot(M1simFPDstn, M1estFPDstn, xlab = "Simulated FPD scores", ylab = "Reconstructed FPD scores", xlim = range(.0, .3), ylim = range(0.0, .3), col = rgb(0, 0, 0, .2))
dev.off()
png("Model2_FPD_scores.png")
plot(M2simFPDstn, M2estFPDstn, xlab = "Simulated FPD scores", ylab = "Reconstructed FPD scores", xlim = range(.0, .3), ylim = range(0.0, .3), col = rgb(0, 0, 0, .2))
dev.off()
png("Model3_FPD_scores.png")
plot(M3simFPDstn, M3estFPDstn, xlab = "Simulated FPD scores", ylab = "Reconstructed FPD scores", xlim = range(.0, .3), ylim = range(0.0, .3), col = rgb(0, 0, 0, .2))
dev.off()
png("Model4_FPD_scores.png")
plot(M4simFPDstn, M4estFPDstn, xlab = "Simulated FPD scores", ylab = "Reconstructed FPD scores", xlim = range(.0, .3), ylim = range(0.0, .3), col = rgb(0, 0, 0, .2))
dev.off()
png("Model1_FPD_scores_ranked.png")
plot(M1RankSimFPDstn, M1RankEstFPDstn, xlab = "Ranked simulated FPD scores", ylab = "Ranked reconstructed FPD scores", xlim = range(0, 160), ylim = range(0, 160), col = rgb(0, 0, 0, .2))
dev.off()
png("Model2_FPD_scores_ranked.png")
plot(M2RankSimFPDstn, M2RankEstFPDstn, xlab = "Ranked simulated FPD scores", ylab = "Ranked reconstructed FPD scores", xlim = range(0, 160), ylim = range(0, 160), col = rgb(0, 0, 0, .2))
dev.off()
png("Model3_FPD_scores_ranked.png")
plot(M3RankSimFPDstn, M3RankEstFPDstn, xlab = "Ranked simulated FPD scores", ylab = "Ranked reconstructed FPD scores", xlim = range(0, 160), ylim = range(0, 160), col = rgb(0, 0, 0, .2))
dev.off()
png("Model4_FPD_scores_ranked.png")
plot(M4RankSimFPDstn, M4RankEstFPDstn, xlab = "Ranked simulated FPD scores", ylab = "Ranked reconstructed FPD scores", xlim = range(0, 160), ylim = range(0, 160), col = rgb(0, 0, 0, .2))
dev.off()

save(M1simFPDstn, file = "Model1_FPD_Sim_values.Rdata")
save(M1estFPDstn, file = "Model1_FPD_Est_values.Rdata")
save(M2simFPDstn, file = "Model2_FPD_Sim_values.Rdata")
save(M2estFPDstn, file = "Model2_FPD_Est_values.Rdata")
save(M3simFPDstn, file = "Model3_FPD_Sim_values.Rdata")
save(M3estFPDstn, file = "Model3_FPD_Est_values.Rdata")
save(M4simFPDstn, file = "Model4_FPD_Sim_values.Rdata")
save(M4estFPDstn, file = "Model4_FPD_Est_values.Rdata")
save(M1RankSimFPDstn, file = "Model1_FPD_SimRanked_values.Rdata")
save(M1RankEstFPDstn, file = "Model1_FPD_EstRanked_values.Rdata")
save(M2RankSimFPDstn, file = "Model2_FPD_SimRanked_values.Rdata")
save(M2RankEstFPDstn, file = "Model2_FPD_EstRanked_values.Rdata")
save(M3RankSimFPDstn, file = "Model3_FPD_SimRanked_values.Rdata")
save(M3RankEstFPDstn, file = "Model3_FPD_EstRanked_values.Rdata")
save(M4RankSimFPDstn, file = "Model4_FPD_SimRanked_values.Rdata")
save(M4RankEstFPDstn, file = "Model4_FPD_EstRanked_values.Rdata")

setwd("../../../..")

proc.time() - ptm