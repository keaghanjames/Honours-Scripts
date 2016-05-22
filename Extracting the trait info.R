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
#NA the correlation test will not work because there is no variance in the substitution rate between tips. 

#Model2
M2subrates <- vector()
M2EDproportion <- vector()

for(i in 1: length(folders2)){
  #print(getwd())
  setwd (folders2[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders2[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  load(paste0(folders2[i], ".Rdata"))
  
simtermstraces <- tree.sim$timephyloFULL$edge[, 2] <= Ntip(tree.sim$timephyloFULL)
tiptraces <- tree.sim$traces[simtermstraces,]
tip.labels.full <- tree.sim$timephyloFULL$tip.label
tip.labels <- tree.sim$timephylo$tip.label
tracesFULL <- cbind(tip.labels.full, tiptraces)
rownames(tracesFULL) <- tracesFULL$tip.labels.full
tracesCLIPED <- tracesFULL[tip.labels,]

#V4 = trait, V5 = substitution rate, V6 = speciation rate

#we create an object which contains only the trait values
subrates <- tracesCLIPED$V5

#scale and centre subrates
subrates <- scale(subrates, center = T, scale = T)
#we calculate ED scores for all of the tips on the simulated and estimated trees

TBLsimtree <- pd.calc(simtree)
TBLesttree <- pd.calc(esttree)
simED <- ed.calc(simtree)
simEDstn <- (simED$spp$ED)/TBLsimtree
estED <- ed.calc(esttree)
estEDstn <- (estED$spp$ED)/TBLesttree

#we create and object which contains the difference between the PD value of interest on the simulated and estimated trees
EDproportion <- (estEDstn/simEDstn)

#we scale and centre ED proportion
EDproportion <- scale(EDproportion, center = T, scale = T)

#we add the scaled subrates and the scaled EDproportion to their respective vectors
M2subrates <- c(M2subrates, subrates)
M2EDproportion <- c(M2EDproportion, EDproportion)

setwd("..")
print("win")
}

#Model3
M3subrates <- vector()
M3EDproportion <- vector()

for(i in 1: length(folders3)){
  #print(getwd())
  setwd (folders3[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders3[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  load(paste0(folders3[i], ".Rdata"))
  
  simtermstraces <- tree.sim$timephyloFULL$edge[, 2] <= Ntip(tree.sim$timephyloFULL)
  tiptraces <- tree.sim$traces[simtermstraces,]
  tip.labels.full <- tree.sim$timephyloFULL$tip.label
  tip.labels <- tree.sim$timephylo$tip.label
  tracesFULL <- cbind(tip.labels.full, tiptraces)
  rownames(tracesFULL) <- tracesFULL$tip.labels.full
  tracesCLIPED <- tracesFULL[tip.labels,]
  
  #V4 = trait, V5 = substitution rate, V6 = speciation rate
  
  #we create an object which contains only the trait values
  subrates <- tracesCLIPED$V5
  
  #scale and centre subrates
  subrates <- scale(subrates, center = T, scale = T)
  #we calculate ED scores for all of the tips on the simulated and estimated trees
  
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simEDstn <- (simED$spp$ED)/TBLsimtree
  estED <- ed.calc(esttree)
  estEDstn <- (estED$spp$ED)/TBLesttree
  
  #we create and object which contains the difference between the PD value of interest on the simulated and estimated trees
  EDproportion <- (estEDstn/simEDstn)
  
  #we scale and centre ED proportion
  EDproportion <- scale(EDproportion, center = T, scale = T)
  
  #we add the scaled subrates and the scaled EDproportion to their respective vectors
  M3subrates <- c(M3subrates, subrates)
  M3EDproportion <- c(M3EDproportion, EDproportion)
  
  setwd("..")
  print("win")
}

#Model4
M4subrates <- vector()
M4EDproportion <- vector()

for(i in 1: length(folders4)){
  #print(getwd())
  setwd (folders4[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders4[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  load(paste0(folders4[i], ".Rdata"))
  
  simtermstraces <- tree.sim$timephyloFULL$edge[, 2] <= Ntip(tree.sim$timephyloFULL)
  tiptraces <- tree.sim$traces[simtermstraces,]
  tip.labels.full <- tree.sim$timephyloFULL$tip.label
  tip.labels <- tree.sim$timephylo$tip.label
  tracesFULL <- cbind(tip.labels.full, tiptraces)
  rownames(tracesFULL) <- tracesFULL$tip.labels.full
  tracesCLIPED <- tracesFULL[tip.labels,]
  
  #V4 = trait, V5 = substitution rate, V6 = speciation rate
  
  #we create an object which contains only the trait values
  subrates <- tracesCLIPED$V5
  
  #scale and centre subrates
  subrates <- scale(subrates, center = T, scale = T)
  #we calculate ED scores for all of the tips on the simulated and estimated trees
  
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simEDstn <- (simED$spp$ED)/TBLsimtree
  estED <- ed.calc(esttree)
  estEDstn <- (estED$spp$ED)/TBLesttree
  
  #we create and object which contains the difference between the PD value of interest on the simulated and estimated trees
  EDproportion <- (estEDstn/simEDstn)
  
  #we scale and centre ED proportion
  EDproportion <- scale(EDproportion, center = T, scale = T)
  
  #we add the scaled subrates and the scaled EDproportion to their respective vectors
  M4subrates <- c(M4subrates, subrates)
  M4EDproportion <- c(M4EDproportion, EDproportion)
  
  setwd("..")
  print("win")
}

#Model5
M5subrates <- vector()
M5EDproportion <- vector()

for(i in 1: length(folders5)){
  #print(getwd())
  setwd (folders5[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders5[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  load(paste0(folders5[i], ".Rdata"))
  
  simtermstraces <- tree.sim$timephyloFULL$edge[, 2] <= Ntip(tree.sim$timephyloFULL)
  tiptraces <- tree.sim$traces[simtermstraces,]
  tip.labels.full <- tree.sim$timephyloFULL$tip.label
  tip.labels <- tree.sim$timephylo$tip.label
  tracesFULL <- cbind(tip.labels.full, tiptraces)
  rownames(tracesFULL) <- tracesFULL$tip.labels.full
  tracesCLIPED <- tracesFULL[tip.labels,]
  
  #V4 = trait, V5 = substitution rate, V6 = speciation rate
  
  #we create an object which contains only the trait values
  subrates <- tracesCLIPED$V5
  
  #scale and centre subrates
  subrates <- scale(subrates, center = T, scale = T)
  #we calculate ED scores for all of the tips on the simulated and estimated trees
  
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simEDstn <- (simED$spp$ED)/TBLsimtree
  estED <- ed.calc(esttree)
  estEDstn <- (estED$spp$ED)/TBLesttree
  
  #we create and object which contains the difference between the PD value of interest on the simulated and estimated trees
  EDproportion <- (estEDstn/simEDstn)
  
  #we scale and centre ED proportion
  EDproportion <- scale(EDproportion, center = T, scale = T)
  
  #we add the scaled subrates and the scaled EDproportion to their respective vectors
  M5subrates <- c(M5subrates, subrates)
  M5EDproportion <- c(M5EDproportion, EDproportion)
  
  setwd("..")
  print("win")
}