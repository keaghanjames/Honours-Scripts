M1EDpvalues <- vector()
M2EDpvalues <- vector()
M3EDpvalues <- vector()
M4EDpvalues <- vector()
M5EDpvalues <- vector()
M1EDVvalues <- vector()
M2EDVvalues <- vector()
M3EDVvalues <- vector()
M4EDVvalues <- vector()
M5EDVvalues <- vector()
M1EDConfInt <- vector()
M2EDConfInt <- vector()
M3EDConfInt <- vector()
M4EDConfInt <- vector()
M1simEDstn <- vector()
M1estEDstn <- vector()
M2simEDstn <- vector()
M2estEDstn <- vector()
M3simEDstn <- vector()
M3estEDstn <- vector()
M4simEDstn <- vector()
M4estEDstn <- vector()
M1RankSimEDstn <- vector()
M1RankEstEDstn <- vector()
M2RankSimEDstn <- vector()
M2RankEstEDstn <- vector()
M3RankSimEDstn <- vector()
M3RankEstEDstn <- vector()
M4RankSimEDstn <- vector()
M4RankEstEDstn <- vector()
n <- 0

filename <- read.csv('Settings.csv')
filename <- filename$runName
simlist <- (list())
estlist <- (list())
folders1 <- grep("Model1", dir(), value = T)
folders2 <- grep("Model2", dir(), value = T)
folders3 <- grep("Model3", dir(), value = T)
folders4 <- grep("Model4", dir(), value = T)
folders5 <- grep("Model5", dir(), value = T)


#Model1
for(i in 1: length(folders1)){
  #print(getwd())
  setwd (folders1[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders1[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simEDstn <- (simED$spp$ED)/TBLsimtree
  estED <- ed.calc(esttree)
  estEDstn <- (estED$spp$ED)/TBLesttree
  RankSimEDstn <- rank(simEDstn)
  RankEstEDstn <- rank(estEDstn)
  M1RankSimEDstn <- c(M1RankSimEDstn, RankSimEDstn)
  M1RankEstEDstn <- c(M1RankEstEDstn, RankEstEDstn)
  M1simEDstn <- c(M1simEDstn, simEDstn)
  M1estEDstn <- c(M1estEDstn, estEDstn)
  EDWilcox <- wilcox.test(simEDstn, estEDstn, paired = T, conf.int = T)
  EDttest <- t.test(simEDstn, estEDstn, paired = T)
  EDWilcox
  EDttest
  SWTestSim <- shapiro.test(simEDstn)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(estEDstn)
  SWTestEstP <- SWTestEst$p.value
  p <- EDWilcox$p.value
  M1EDpvalues <- c(M1EDpvalues, p)
  M1EDVvalues <- c(M1EDVvalues, (EDWilcox$statistic))
  M1EDConfInt <- c(M1EDConfInt, (EDWilcox$conf.int))
  #print(getwd())
  setwd("..")
  #print(getwd())
  print("chaching")
  n <- n+1
  print(n)
}

M2simED <- vector()
M2estED <- vector()

#Model2
for(i in 1: length(folders2)){
  #print(getwd())
  setwd (folders2[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders2[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simEDstn <- (simED$spp$ED)/TBLsimtree
  simED <- simED$spp$ED
  estED <- ed.calc(esttree)
  estEDstn <- (estED$spp$ED)/TBLesttree
  estED <- estED$spp$ED
  RankSimEDstn <- rank(simEDstn)
  RankEstEDstn <- rank(estEDstn)
  M2RankSimEDstn <- c(M2RankSimEDstn, RankSimEDstn)
  M2RankEstEDstn <- c(M2RankEstEDstn, RankEstEDstn)
  M2simEDstn <- c(M2simEDstn, simEDstn)
  M2estEDstn <- c(M2estEDstn, estEDstn)
  M2simED <- c(M2simED, simED)
  M2estED <- c(M2estED, estED)
  EDWilcox <- wilcox.test(simEDstn, estEDstn, paired = T, conf.int = T)
  EDttest <- t.test(simEDstn, estEDstn, paired = T)
  EDWilcox
  EDttest
  SWTestSim <- shapiro.test(simEDstn)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(estEDstn)
  SWTestEstP <- SWTestEst$p.value
  p <- EDWilcox$p.value
  M2EDpvalues <- c(M2EDpvalues, p)
  M2EDVvalues <- c(M2EDVvalues, (EDWilcox$statistic))
  M2EDConfInt <- c(M2EDConfInt, (EDWilcox$conf.int))
  #print(getwd())
  setwd("..")
  #print(getwd())
  print("chaching")
  n <- n+1
  print(n)
}

#Model3
for(i in 1: length(folders3)){
  #print(getwd())
  setwd (folders3[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders3[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simEDstn <- (simED$spp$ED)/TBLsimtree
  estED <- ed.calc(esttree)
  estEDstn <- (estED$spp$ED)/TBLesttree
  RankSimEDstn <- rank(simEDstn)
  RankEstEDstn <- rank(estEDstn)
  M3RankSimEDstn <- c(M3RankSimEDstn, RankSimEDstn)
  M3RankEstEDstn <- c(M3RankEstEDstn, RankEstEDstn)
  M3simEDstn <- c(M3simEDstn, simEDstn)
  M3estEDstn <- c(M3estEDstn, estEDstn)
  EDWilcox <- wilcox.test(simEDstn, estEDstn, paired = T)
  EDttest <- t.test(simEDstn, estEDstn, paired = T)
  EDWilcox
  EDttest
  SWTestSim <- shapiro.test(simEDstn)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(estEDstn)
  SWTestEstP <- SWTestEst$p.value
  p <- EDWilcox$p.value
  M3EDpvalues <- c(M3EDpvalues, p)
  M3EDVvalues <- c(M3EDVvalues, (EDWilcox$statistic))
  M3EDConfInt <- c(M3EDConfInt, (EDWilcox$conf.int))
  #print(getwd())
  setwd("..")
  #print(getwd())
  print("chaching")
  n <- n+1
  print(n)
}

#Model4
for(i in 1: length(folders4)){
  #print(getwd())
  setwd (folders4[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders4[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simEDstn <- (simED$spp$ED)/TBLsimtree
  estED <- ed.calc(esttree)
  estEDstn <- (estED$spp$ED)/TBLesttree
  RankSimEDstn <- rank(simEDstn)
  RankEstEDstn <- rank(estEDstn)
  M4RankSimEDstn <- c(M4RankSimEDstn, RankSimEDstn)
  M4RankEstEDstn <- c(M4RankEstEDstn, RankEstEDstn)
  M4simEDstn <- c(M4simEDstn, simEDstn)
  M4estEDstn <- c(M4estEDstn, estEDstn)
  EDWilcox <- wilcox.test(simEDstn, estEDstn, paired = T, conf.int = T)
  EDttest <- t.test(simEDstn, estEDstn, paired = T)
  EDWilcox
  EDttest
  SWTestSim <- shapiro.test(simEDstn)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(estEDstn)
  SWTestEstP <- SWTestEst$p.value
  p <- EDWilcox$p.value
  M4EDpvalues <- c(M4EDpvalues, p)
  M4EDVvalues <- c(M4EDVvalues, (EDWilcox$statistic))
  M4EDConfInt <- c(M4EDConfInt, (EDWilcox$conf.int))
  #print(getwd())
  setwd("..")
  #print(getwd())
  print("chaching")
  n <- n+1
  print(n)
}


setwd("results/ED")
save(M1EDpvalues, file = "Model1_ED_p.values.Rdata")
save(M2EDpvalues, file = "Model2_ED_p.values.Rdata")
save(M3EDpvalues, file = "Model3_ED_p.values.Rdata")
save(M4EDpvalues, file = "Model4_ED_p.values.Rdata")
save(M1EDVvalues, file = "Model1_ED_Vvalues.Rdata")
save(M2EDVvalues, file = "Model2_ED_Vvalues.Rdata")
save(M3EDVvalues, file = "Model3_ED_Vvalues.Rdata")
save(M4EDVvalues, file = "Model4_ED_Vvalues.Rdata")
save(M1simEDstn, file = "Model1_ED_Sim_values.Rdata")
save(M1estEDstn, file = "Model1_ED_Est_values.Rdata")
save(M2simEDstn, file = "Model2_ED_Sim_values.Rdata")
save(M2estEDstn, file = "Model2_ED_Est_values.Rdata")
save(M3simEDstn, file = "Model3_ED_Sim_values.Rdata")
save(M3estEDstn, file = "Model3_ED_Est_values.Rdata")
save(M4simEDstn, file = "Model4_ED_Sim_values.Rdata")
save(M4estEDstn, file = "Model4_ED_Est_values.Rdata")
save(M1RankSimEDstn, file = "Model1_ED_SimRanked_values.Rdata")
save(M1RankEstEDstn, file = "Model1_ED_EstRanked_values.Rdata")
save(M2RankSimEDstn, file = "Model2_ED_SimRanked_values.Rdata")
save(M2RankEstEDstn, file = "Model2_ED_EstRanked_values.Rdata")
save(M3RankSimEDstn, file = "Model3_ED_SimRanked_values.Rdata")
save(M3RankEstEDstn, file = "Model3_ED_EstRanked_values.Rdata")
save(M4RankSimEDstn, file = "Model4_ED_SimRanked_values.Rdata")
save(M4RankEstEDstn, file = "Model4_ED_EstRanked_values.Rdata")

#plots a histogram of p.value frequency, groupted by breaks of 0.05 for each of the Models
#Model1
png("Model1_ED_p.values_hist.png")
hist(M1EDpvalues, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

#Model2
png("Model2_ED_p.values_hist.png")
hist(M2EDpvalues, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

#Model3
png("Model3_ED_p.values_hist.png")
hist(M3EDpvalues, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

#Model4
png("Model4_ED_p.values_hist.png")
hist(M4EDpvalues, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

#creates scatter plots of the simulated and reconstructed ED scores and ranked ED scores
png("Model1_ED_scores.png")
plot(M1simEDstn, M1estEDstn, xlab = "Simulated ED scores", ylab = "Reconstructed ED scores", xlim = range(0, 0.04), ylim = range(0, 0.04), col = rgb(0, 0, 0, .2))
dev.off()
png("Model2_ED_scores.png")
plot(M2simEDstn, M2estEDstn, xlab = "Simulated ED scores", ylab = "Reconstructed ED scores", xlim = range(0, 0.04), ylim = range(0, 0.04), col = rgb(0, 0, 0, .2))
dev.off()
png("Model3_ED_scores.png")
plot(M3simEDstn, M3estEDstn, xlab = "Simulated ED scores", ylab = "Reconstructed ED scores", xlim = range(0, 0.04), ylim = range(0, 0.04), col = rgb(0, 0, 0, .2))
dev.off()
png("Model4_ED_scores.png")
plot(M4simEDstn, M4estEDstn, xlab = "Simulated ED scores", ylab = "Reconstructed ED scores", xlim = range(0, 0.04), ylim = range(0, 0.04), col = rgb(0, 0, 0, .2))
dev.off()
png("Model1_ED_scores_ranked.png")
plot(M1RankSimEDstn, M1RankEstEDstn, xlab = "Ranked simulated ED scores", ylab = "Ranked reconstructed ED scores", xlim = range(0, 210), ylim = range(0, 210), col = rgb(0, 0, 0, .2))
dev.off()
png("Model2_ED_scores_ranked.png")
plot(M2RankSimEDstn, M2RankEstEDstn, xlab = "Ranked simulated ED scores", ylab = "Ranked reconstructed ED scores", xlim = range(0, 210), ylim = range(0, 210), col = rgb(0, 0, 0, .2))
dev.off()
png("Model3_ED_scores_ranked.png")
plot(M3RankSimEDstn, M3RankEstEDstn, xlab = "Ranked simulated ED scores", ylab = "Ranked reconstructed ED scores", xlim = range(0, 210), ylim = range(0, 210), col = rgb(0, 0, 0, .2))
dev.off()
png("Model4_ED_scores_ranked.png")
plot(M4RankSimEDstn, M4RankEstEDstn, xlab = "Ranked simulated ED scores", ylab = "Ranked reconstructed ED scores", xlim = range(0, 210), ylim = range(0, 210), col = rgb(0, 0, 0, .2))
dev.off()

setwd("../..")