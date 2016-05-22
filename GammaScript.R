M1SimGamaStat <- vector()
M1EstGamaStat <- vector()
M2SimGamaStat <- vector()
M2EstGamaStat <- vector()
M3SimGamaStat <- vector()
M3EstGamaStat <- vector()
M4SimGamaStat <- vector()
M4EstGamaStat <- vector()
M5SimGamaStat <- vector()
M5EstGamaStat <- vector()
MeanSim <- vector()
SDSim <- vector()
MeanEst <- vector()
SDEst <- vector()

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
  SG <- gammaStat(simtree)
  M1SimGamaStat <- c(M1SimGamaStat, SG)
  EG <- gammaStat(esttree)
  M1EstGamaStat <- c(M1EstGamaStat, EG)
  setwd("..")
  print("!")
}
MeanSim <- c(MeanSim, (mean(M1SimGamaStat)))
SDSim <- c(SDSim, (sd(M1SimGamaStat)))
MeanEst <- c(MeanEst, (mean(M1EstGamaStat)))
SDEst <- c(SDEst, (sd(M1EstGamaStat)))

setwd("results/Gamma")
M1ttest <- t.test(M1SimGamaStat, M1EstGamaStat, paired = T)
save(M1ttest, file = "M1ttest.Rdata")
print(M1ttest)

png("Model1_Sim_Gamma_hist.png")
hist(M1SimGamaStat, col = 'lightblue')
dev.off()

png("Model1_Est_Gamma_hist.png")
hist(M1EstGamaStat, col = 'lightblue')
dev.off()

setwd("../..")

#Model2

for(i in 1: length(folders2)){
  #print(getwd())
  setwd (folders2[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders2[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  SG <- gammaStat(simtree)
  M2SimGamaStat <- c(M2SimGamaStat, SG)
  EG <- gammaStat(esttree)
  M2EstGamaStat <- c(M2EstGamaStat, EG)
  setwd("..")
  print("!")
}

setwd("results/Gamma")
M2ttest <- t.test(M2SimGamaStat, M2EstGamaStat, paired = T)
save(M2ttest, file = "M2ttest.Rdata")
print(M2ttest)

MeanSim <- c(MeanSim, (mean(M2SimGamaStat)))
SDSim <- c(SDSim, (sd(M2SimGamaStat)))
MeanEst <- c(MeanEst, (mean(M2EstGamaStat)))
SDEst <- c(SDEst, (sd(M2EstGamaStat)))

png("Model2_Sim_Gamma_hist.png")
hist(M2SimGamaStat, col = 'lightblue')
dev.off()

png("Model2_Est_Gamma_hist.png")
hist(M2EstGamaStat, col = 'lightblue')
dev.off()

setwd("../..")

#Model3

for(i in 1: length(folders3)){
  #print(getwd())
  setwd (folders3[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders3[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  SG <- gammaStat(simtree)
  M3SimGamaStat <- c(M3SimGamaStat, SG)
  EG <- gammaStat(esttree)
  M3EstGamaStat <- c(M3EstGamaStat, EG)
  setwd("..")
  print("!")
}

setwd("results/Gamma")
M3ttest <- t.test(M3SimGamaStat, M3EstGamaStat, paired = T)
save(M3ttest, file = "M3ttest.Rdata")
print(M3ttest)

MeanSim <- c(MeanSim, (mean(M3SimGamaStat)))
SDSim <- c(SDSim, (sd(M3SimGamaStat)))
MeanEst <- c(MeanEst, (mean(M3EstGamaStat)))
SDEst <- c(SDEst, (sd(M3EstGamaStat)))

png("Model3_Sim_Gamma_hist.png")
hist(M3SimGamaStat, col = 'lightblue')
dev.off()

png("Model3_Est_Gamma_hist.png")
hist(M3EstGamaStat, col = 'lightblue')
dev.off()

setwd("../..")

#Model4

for(i in 1: length(folders4)){
  #print(getwd())
  setwd (folders4[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders4[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  SG <- gammaStat(simtree)
  M4SimGamaStat <- c(M4SimGamaStat, SG)
  EG <- gammaStat(esttree)
  M4EstGamaStat <- c(M4EstGamaStat, EG)
  setwd("..")
  print("!")
}

setwd("results/Gamma")
M4ttest <- t.test(M4SimGamaStat, M4EstGamaStat, paired = T)
save(M4ttest, file = "M4ttest.Rdata")
print(M4ttest)

MeanSim <- c(MeanSim, (mean(M4SimGamaStat)))
SDSim <- c(SDSim, (sd(M4SimGamaStat)))
MeanEst <- c(MeanEst, (mean(M4EstGamaStat)))
SDEst <- c(SDEst, (sd(M4EstGamaStat)))

png("Model4_Sim_Gamma_hist.png")
hist(M4SimGamaStat, col = 'lightblue')
dev.off()

png("Model4_Est_Gamma_hist.png")
hist(M4EstGamaStat, col = 'lightblue')
dev.off()

setwd("../..")

#Model5

for(i in 1: length(folders5)){
  #print(getwd())
  setwd (folders5[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders5[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  SG <- gammaStat(simtree)
  M5SimGamaStat <- c(M5SimGamaStat, SG)
  EG <- gammaStat(esttree)
  M5EstGamaStat <- c(M5EstGamaStat, EG)
  setwd("..")
  print("!")
}

setwd("results/Gamma")
M5ttest <- t.test(M5SimGamaStat, M5EstGamaStat, paired = T)
save(M5ttest, file = "M5ttest.Rdata")
print(M5ttest)

png("Model5_Sim_Gamma_hist.png")
hist(M5SimGamaStat, col = 'lightblue')
dev.off()

png("Model5_Est_Gamma_hist.png")
hist(M5EstGamaStat, col = 'lightblue')
dev.off()

GammaTable <- data.frame(MeanSim, SDSim, MeanEst, SDEst)
write.csv(GammaTable, file = "GammaTable.csv")

setwd("../..")