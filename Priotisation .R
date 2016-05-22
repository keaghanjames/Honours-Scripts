filename <- read.csv('Settings.csv')
filename <- filename$runName
simlist <- (list())
estlist <- (list())
folders1 <- grep("Model1", dir(), value = T)
folders2 <- grep("Model2", dir(), value = T)
folders3 <- grep("Model3", dir(), value = T)
folders4 <- grep("Model4", dir(), value = T)
folders5 <- grep("Model5", dir(), value = T)
M1EDrank <- vector()
M2EDrank <- vector()
M3EDrank <- vector()
M4EDrank <- vector()
n <- 0

#model1
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
  
  simEDscores <- simED$spp
  simEDscores <- simEDscores[sort.list(simEDscores[,2], decreasing = T), ]
  simEDscores <- simEDscores[, 1]
  simEDscores <- simEDscores[1:20]
  estEDscores <- estED$spp
  estEDscores <- estEDscores[sort.list(estEDscores[,2], decreasing = T), ]
  estEDscores <- estEDscores[, 1]
  estEDscores <- estEDscores[1:20]
  x <- setdiff(simEDscores, estEDscores)
  x <- length(x)
  M1EDrank <- c(M1EDrank, x)
  n <- n+1
  print(n)
  setwd("..")
}

#model2
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
  estED <- ed.calc(esttree)
  estEDstn <- (estED$spp$ED)/TBLesttree

simEDscores <- simED$spp
simEDscores <- simEDscores[sort.list(simEDscores[,2], decreasing = T), ]
simEDscores <- simEDscores[, 1]
simEDscores <- simEDscores[1:20]
estEDscores <- estED$spp
estEDscores <- estEDscores[sort.list(estEDscores[,2], decreasing = T), ]
estEDscores <- estEDscores[, 1]
estEDscores <- estEDscores[1:20]
x <- setdiff(simEDscores, estEDscores)
x <- length(x)
M2EDrank <- c(M2EDrank, x)
n <- n+1
print(n)
setwd("..")
}

#model3
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
  
  #having calculated the ED scores I then creates two matrices, one with the 
  simEDscores <- simED$spp
  simEDscores <- simEDscores[sort.list(simEDscores[,2], decreasing = T), ]
  simEDscores <- simEDscores[, 1]
  simEDscores <- simEDscores[1:20]
  estEDscores <- estED$spp
  estEDscores <- estEDscores[sort.list(estEDscores[,2], decreasing = T), ]
  estEDscores <- estEDscores[, 1]
  estEDscores <- estEDscores[1:20]
  x <- setdiff(simEDscores, estEDscores)
  x <- length(x)
  M3EDrank <- c(M2EDrank, x)
  n <- n+1
  print(n)
  setwd("..")
}

#model4
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
  
  #having calculated the ED scores I then creates two matrices, one with the 
  simEDscores <- simED$spp
  simEDscores <- simEDscores[sort.list(simEDscores[,2], decreasing = T), ]
  simEDscores <- simEDscores[, 1]
  simEDscores <- simEDscores[1:20]
  estEDscores <- estED$spp
  estEDscores <- estEDscores[sort.list(estEDscores[,2], decreasing = T), ]
  estEDscores <- estEDscores[, 1]
  estEDscores <- estEDscores[1:20]
  x <- setdiff(simEDscores, estEDscores)
  x <- length(x)
  M4EDrank <- c(M2EDrank, x)
  n <- n+1
  print(n)
  setwd("..")
}

