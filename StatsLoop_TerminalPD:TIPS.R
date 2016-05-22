#should you be standerdising by the total length of the tree or the total length the tips

M1TPDpvalues <- vector()
M2TPDpvalues <- vector()
M3TPDpvalues <- vector()
M4TPDpvalues <- vector()
M1TPDVvalues <- vector()
M2TPDVvalues <- vector()
M3TPDVvalues <- vector()
M4TPDVvalues <- vector()
M1TPDdif <- vector()
M2TPDdif <- vector()
M3TPDdif <- vector()
M4TPDdif <- vector()

filename <- read.csv('Settings.csv')
filename <- filename$runName
simlist <- (list())
estlist <- (list())
folders1 <- grep("Model1", dir(), value = T)
folders2 <- grep("Model2", dir(), value = T)
folders3 <- grep("Model3", dir(), value = T)
folders4 <- grep("Model4", dir(), value = T)

#Model1
for(i in 1: length(folders1)){
  #print(getwd())
  setwd (folders1[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders1[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree, method = "TIP")
  TBLesttree <- pd.calc(esttree, method = "TIP")
  simterms <- simtree$edge[, 2] <= Ntip(simtree)
  simTPD <- simtree$edge.length[simterms]
  simTPDstn <- (simTPD)/TBLsimtree
  estterms <- esttree$edge[, 2] <= Ntip(esttree)
  estTPD <- esttree$edge.length[estterms]
  estTPDstn <- (estTPD)/TBLesttree
  TPDdif <- estTPDstn-simTPDstn
  TPDdif <- median(TPDdif)
  TPDWilcox <- wilcox.test(simTPDstn, estTPDstn, paired = T)
  TPDttest <- t.test(simTPDstn, estTPDstn, paired = T)
  TPDWilcox
  TPDttest
  SWTestSim <- shapiro.test(simTPDstn)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(estTPDstn)
  SWTestEstP <- SWTestEst$p.value
  p <- TPDWilcox$p.value
  M1TPDpvalues <- c(M1TPDpvalues, p)
  M1TPDVvalues <- c(M1TPDVvalues, (TPDWilcox$statistic))
  M1TPDdif <- c(M1TPDdif, TPDdif)
  #print(getwd())
  setwd("..")
  #print(getwd())
  print("win")
}

#Model2
for(i in 1: length(folders2)){
  #print(getwd())
  setwd (folders2[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders2[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree, method = "TIP")
  TBLesttree <- pd.calc(esttree, method = "TIP")
  simterms <- simtree$edge[, 2] <= Ntip(simtree)
  simTPD <- simtree$edge.length[simterms]
  simTPDstn <- (simTPD)/TBLsimtree
  estterms <- esttree$edge[, 2] <= Ntip(esttree)
  estTPD <- esttree$edge.length[estterms]
  estTPDstn <- (estTPD)/TBLesttree
  TPDdif <- estTPDstn-simTPDstn
  TPDdif <- median(TPDdif)
  TPDWilcox <- wilcox.test(simTPDstn, estTPDstn, paired = T)
  TPDttest <- t.test(simTPDstn, estTPDstn, paired = T)
  TPDWilcox
  TPDttest
  SWTestSim <- shapiro.test(simTPDstn)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(estTPDstn)
  SWTestEstP <- SWTestEst$p.value
  p <- TPDWilcox$p.value
  M2TPDpvalues <- c(M2TPDpvalues, p)
  M2TPDVvalues <- c(M2TPDVvalues, (TPDWilcox$statistic))
  M2TPDdif <- c(M2TPDdif, TPDdif)
  #print(getwd())
  setwd("..")
  #print(getwd())
  print("win")
}

#Model3
for(i in 1: length(folders3)){
  #print(getwd())
  setwd (folders3[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders3[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree, method = "TIP")
  TBLesttree <- pd.calc(esttree, method = "TIP")
  simterms <- simtree$edge[, 2] <= Ntip(simtree)
  simTPD <- simtree$edge.length[simterms]
  simTPDstn <- (simTPD)/TBLsimtree
  estterms <- esttree$edge[, 2] <= Ntip(esttree)
  estTPD <- esttree$edge.length[estterms]
  estTPDstn <- (estTPD)/TBLesttree
  TPDdif <- estTPDstn-simTPDstn
  TPDdif <- median(TPDdif)
  TPDWilcox <- wilcox.test(simTPDstn, estTPDstn, paired = T)
  TPDttest <- t.test(simTPDstn, estTPDstn, paired = T)
  TPDWilcox
  TPDttest
  SWTestSim <- shapiro.test(simTPDstn)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(estTPDstn)
  SWTestEstP <- SWTestEst$p.value
  p <- TPDWilcox$p.value
  M3TPDpvalues <- c(M3TPDpvalues, p)
  M3TPDVvalues <- c(M3TPDVvalues, (TPDWilcox$statistic))
  M3TPDdif <- c(M3TPDdif, TPDdif)
  #print(getwd())
  setwd("..")
  #print(getwd())
  print("win")
}

#Model4
for(i in 1: length(folders4)){
  #print(getwd())
  setwd (folders4[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders4[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  TBLsimtree <- pd.calc(simtree, method = "TIP")
  TBLesttree <- pd.calc(esttree, method = "TIP")
  simterms <- simtree$edge[, 2] <= Ntip(simtree)
  simTPD <- simtree$edge.length[simterms]
  simTPDstn <- (simTPD)/TBLsimtree
  estterms <- esttree$edge[, 2] <= Ntip(esttree)
  estTPD <- esttree$edge.length[estterms]
  estTPDstn <- (estTPD)/TBLesttree
  TPDdif <- estTPDstn-simTPDstn
  TPDdif <- median(TPDdif)
  TPDWilcox <- wilcox.test(simTPDstn, estTPDstn, paired = T)
  TPDttest <- t.test(simTPDstn, estTPDstn, paired = T)
  TPDWilcox
  TPDttest
  SWTestSim <- shapiro.test(simTPDstn)
  SWTestSimP <- SWTestSim$p.value
  SWTestEst <- shapiro.test(estTPDstn)
  SWTestEstP <- SWTestEst$p.value
  p <- TPDWilcox$p.value
  M4TPDpvalues <- c(M4TPDpvalues, p)
  M4TPDVvalues <- c(M4TPDVvalues, (TPDWilcox$statistic))
  M4TPDdif <- c(M1TPDdif, TPDdif)
  #print(getwd())
  setwd("..")
  #print(getwd())
  print("win")
}

setwd("results/TPD/DividedByTIPS")
save(M1TPDpvalues, file = "Model1_TPD_p.values.Rdata")
save(M2TPDpvalues, file = "Model2_TPD_p.values.Rdata")
save(M3TPDpvalues, file = "Model3_TPD_p.values.Rdata")
save(M4TPDpvalues, file = "Model4_TPD_p.values.Rdata")
save(M1TPDVvalues, file = "Model1_TPD_Vvalues.Rdata")
save(M2TPDVvalues, file = "Model2_TPD_Vvalues.Rdata")
save(M3TPDVvalues, file = "Model3_TPD_Vvalues.Rdata")
save(M4TPDVvalues, file = "Model4_TPD_Vvalues.Rdata")
save(M1TPDdif, file = "Model1_TPD_median_differences.Rdata")
save(M2TPDdif, file = "Model1_TPD_median_differences.Rdata")
save(M3TPDdif, file = "Model1_TPD_median_differences.Rdata")
save(M4TPDdif, file = "Model1_TPD_median_differences.Rdata")

#plots a histogram of p.value frequency, groupted by breaks of 0.05 for each of the Models
#Model1
png("Model1_TPD_p.values_hist.png")
hist(M1TPDpvalues, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

#Model2
png("Model2_TPD_p.values_hist.png")
hist(M2TPDpvalues, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

#Model3
png("Model3_TPD_p.values_hist.png")
hist(M3TPDpvalues, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

#Model4
png("Model4_TPD_p.values_hist.png")
hist(M4TPDpvalues, col = 'lightblue', breaks = seq(0, 1, 0.05), main = NULL, xlim = c(0,1), ylim = c(0,100), xlab = paste("p.values"), labels = T)
abline(v=.05,col=2,lty=1)
text(0.08, y = 90, labels = "0.05", col = 'red')
dev.off()

setwd("../../..")