require(phangorn)
require(ape)
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

#Model 1

M1RhoDif <- vector()
M1pDif <- vector()
M1RhoProp <- vector()
M1pProp <- vector()
M1EDMedianDif <- vector()
M1EDMedianProp <- vector()
M1EDMedianTV <- vector()
M1RTV <- vector()

for(i in 1: length(folders1)){
  #print(getwd())
  setwd (folders1[i])
  #print(getwd())
  simtree <- read.tree(paste0(folders1[i], "_sim.tre"))
  esttree <- read.nexus("BEASTestimated.tree")
  load(paste0(folders1[i], ".Rdata"))
  
  simtermstraces <- tree.sim$timephyloFULL$edge[, 2] <= Ntip(tree.sim$timephyloFULL)
  tiptraces <- tree.sim$traces[simtermstraces,]
  tip.labels.full <- tree.sim$timephyloFULL$tip.label
  tip.labels <- tree.sim$timephylo$tip.label
  tracesFULL <- cbind(tip.labels.full, tiptraces)
  rownames(tracesFULL) <- tracesFULL$tip.labels.full
  tracesCLIPED <- tracesFULL[tip.labels,]
  traitvalues <- (tracesCLIPED$V4)
  traitvaluesMedian <- median(traitvalues)
  Rangetraitvalues <- range(traitvalues)
  Rangetraitvalues <- diff(Rangetraitvalues)
  M1RTV <- c(M1RTV, Rangetraitvalues)
  M1EDMedianTV <- c(M1EDMedianTV, median(traitvalues))
  
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simED <- simED$spp$ED
  estED <- ed.calc(esttree)
  estED <- estED$spp$ED
  EDproportion <- (estED/simED)
  EDdifference <- (estED - simED)
  M1EDMedianDif <- c(M1EDMedianDif, median(EDdifference))
  M1EDMedianProp <- c(M1EDMedianProp, median(EDproportion))
  
  SpearDif <- cor.test(traitvalues, EDdifference, method = "spearman", alternative = "greater")
  SpearProp <- cor.test(traitvalues, EDproportion, method = "spearman", alternative = "greater")
  rhodif <- SpearDif$estimate
  pdif <- SpearDif$p.value
  rhoprop <- SpearProp$estimate
  pprop <- SpearProp$p.value
  M1RhoDif <- c(M1RhoDif, rhodif)
  M1pDif <- c(M1pDif, pdif)
  M1RhoProp <- c(M1RhoProp, rhoprop)
  M1pProp <- c(M1pProp, pprop)
  print("win")
  setwd("..")
}

setwd("results/ED/Spearman")

M1EDMedianDif <- abs(M1EDMedianDif)
M1EDMedianProp <- abs(M1EDMedianProp)
M1EDMedianTV <- abs(M1EDMedianTV)
M1RhoTrVaEDDif <- cor.test(M1EDMedianTV, M1EDMedianDif, method = "spearman", alternative = "greater")
save(M1RhoTrVaEDDif, file = 'M1Rho_TrVv_EDDif.Rdata')
M1RhoTrVaEDProp <- cor.test(M1EDMedianTV, M1EDMedianProp, method = "spearman", alternative = "greater")
save(M1RhoTrVaEDProp, file = 'M1Rho_TrVv_EDProp.Rdata')
M1RhoRangeEDDif <- cor.test(M1RTV, M1EDMedianDif, method = "spearman", alternative = "greater")
save(M1RhoRangeEDDif, file = 'M1Rho_Range_EDDif.Rdata')
png("M1Rho_Range_EDDif.png")
plot(M1RTV, M1EDMedianDif, xlab = "Trait Value Range", ylab = 'Median ED Error')
dev.off()
M1RhoRangeEDProp <- cor.test(M1RTV, M1EDMedianProp, method = "spearman", alternative = "greater")
save(M1RhoRangeEDProp, file = 'M1Rho_Range_EDProp.Rdata')
png("M1Rho_Range_EDProp.png")
plot(M1RTV, M1EDMedianProp, xlab = "Trait Value Range", ylab = 'Median ED Error')
dev.off()
setwd("../../..")

#Model 2

M2RhoDif <- vector()
M2pDif <- vector()
M2RhoProp <- vector()
M2pProp <- vector()
M2EDMedianDif <- vector()
M2EDMedianProp <- vector()
M2EDMedianTV <- vector()
M2RTV <- vector()

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
traitvalues <- (tracesCLIPED$V4)
traitvaluesMedian <- median(traitvalues)
Rangetraitvalues <- range(traitvalues)
Rangetraitvalues <- diff(Rangetraitvalues)
M2RTV <- c(M2RTV, Rangetraitvalues)
M2EDMedianTV <- c(M2EDMedianTV, median(traitvalues))

TBLsimtree <- pd.calc(simtree)
TBLesttree <- pd.calc(esttree)
simED <- ed.calc(simtree)
simED <- simED$spp$ED
estED <- ed.calc(esttree)
estED <- estED$spp$ED
EDproportion <- (estED/simED)
EDdifference <- (estED - simED)
M2EDMedianDif <- c(M2EDMedianDif, median(EDdifference))
M2EDMedianProp <- c(M2EDMedianProp, median(EDproportion))

SpearDif <- cor.test(traitvalues, EDdifference, method = "spearman", alternative = "greater")
SpearProp <- cor.test(traitvalues, EDproportion, method = "spearman", alternative = "greater")
rhodif <- SpearDif$estimate
pdif <- SpearDif$p.value
rhoprop <- SpearProp$estimate
pprop <- SpearProp$p.value
M2RhoDif <- c(M2RhoDif, rhodif)
M2pDif <- c(M2pDif, pdif)
M2RhoProp <- c(M2RhoProp, rhoprop)
M2pProp <- c(M2pProp, pprop)
print("win")
setwd("..")
}

setwd("results/ED/Spearman")

M2EDMedianDif <- abs(M2EDMedianDif)
M2EDMedianProp <- abs(M2EDMedianProp)
M2EDMedianTV <- abs(M2EDMedianTV)
M2RhoTrVaEDDif <- cor.test(M2EDMedianTV, M2EDMedianDif, method = "spearman", alternative = "greater")
save(M2RhoTrVaEDDif, file = 'M2Rho_TrVv_EDDif.Rdata')
M2RhoTrVaEDProp <- cor.test(M2EDMedianTV, M2EDMedianProp, method = "spearman", alternative = "greater")
save(M2RhoTrVaEDProp, file = 'M2Rho_TrVv_EDProp.Rdata')
M2RhoRangeEDDif <- cor.test(M2RTV, M2EDMedianDif, method = "spearman", alternative = "greater")
save(M2RhoRangeEDDif, file = 'M2Rho_Range_EDDif.Rdata')
png("M2Rho_Range_EDDif.png")
plot(M2RTV, M2EDMedianDif, xlab = "Trait Value Range", ylab = 'Median ED Error')
dev.off()
M2RhoRangeEDProp <- cor.test(M2RTV, M2EDMedianProp, method = "spearman", alternative = "greater")
save(M2RhoRangeEDProp, file = 'M2Rho_Range_EDProp.Rdata')
png("M2Rho_Range_EDProp.png")
plot(M2RTV, M2EDMedianProp, xlab = "Trait Value Range", ylab = 'Median ED Error')
dev.off()
setwd("../../..")

#Model3

M3RhoDif <- vector()
M3pDif <- vector()
M3RhoProp <- vector()
M3pProp <- vector()
M3EDMedianDif <- vector()
M3EDMedianProp <- vector()
M3EDMedianTV <- vector()
M3RTV <- vector()

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
  traitvalues <- (tracesCLIPED$V4)
  traitvaluesMedian <- median(traitvalues)
  Rangetraitvalues <- range(traitvalues)
  Rangetraitvalues <- diff(Rangetraitvalues)
  M3RTV <- c(M3RTV, Rangetraitvalues)
  M3EDMedianTV <- c(M3EDMedianTV, median(traitvalues))
  
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simED <- simED$spp$ED
  estED <- ed.calc(esttree)
  estED <- estED$spp$ED
  EDproportion <- (estED/simED)
  EDdifference <- (estED - simED)
  M3EDMedianDif <- c(M3EDMedianDif, median(EDdifference))
  M3EDMedianProp <- c(M3EDMedianProp, median(EDproportion))
  
  SpearDif <- cor.test(traitvalues, EDdifference, method = "spearman", alternative = "greater")
  SpearProp <- cor.test(traitvalues, EDproportion, method = "spearman", alternative = "greater")
  rhodif <- SpearDif$estimate
  pdif <- SpearDif$p.value
  rhoprop <- SpearProp$estimate
  pprop <- SpearProp$p.value
  M3RhoDif <- c(M3RhoDif, rhodif)
  M3pDif <- c(M3pDif, pdif)
  M3RhoProp <- c(M3RhoProp, rhoprop)
  M3pProp <- c(M3pProp, pprop)
  print("win")
  setwd("..")
}

setwd("results/ED/Spearman")

M3EDMedianDif <- abs(M3EDMedianDif)
M3EDMedianProp <- abs(M3EDMedianProp)
M3EDMedianTV <- abs(M3EDMedianTV)
M3RhoTrVaEDDif <- cor.test(M3EDMedianTV, M3EDMedianDif, method = "spearman", alternative = "greater")
save(M3RhoTrVaEDDif, file = 'M3Rho_TrVv_EDDif.Rdata')
M3RhoTrVaEDProp <- cor.test(M3EDMedianTV, M3EDMedianProp, method = "spearman", alternative = "greater")
save(M3RhoTrVaEDProp, file = 'M3Rho_TrVv_EDProp.Rdata')
M3RhoRangeEDDif <- cor.test(M3RTV, M3EDMedianDif, method = "spearman", alternative = "greater")
save(M3RhoRangeEDDif, file = 'M3Rho_Range_EDDif.Rdata')
png("M3Rho_Range_EDDif.png")
plot(M3RTV, M3EDMedianDif, xlab = "Trait Value Range", ylab = 'Median ED Error')
dev.off()
M3RhoRangeEDProp <- cor.test(M3RTV, M3EDMedianProp, method = "spearman", alternative = "greater")
save(M3RhoRangeEDProp, file = 'M3Rho_Range_EDProp.Rdata')
png("M3Rho_Range_EDProp.png")
plot(M3RTV, M3EDMedianProp, xlab = "Trait Value Range", ylab = 'Median ED Error')
dev.off()
setwd("../../..")

#Model4
M4RhoDif <- vector()
M4pDif <- vector()
M4RhoProp <- vector()
M4pProp <- vector()
M4EDMedianDif <- vector()
M4EDMedianProp <- vector()
M4EDMedianTV <- vector()
M4RTV <- vector()

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
  traitvalues <- (tracesCLIPED$V4)
  traitvaluesMedian <- median(traitvalues)
  Rangetraitvalues <- range(traitvalues)
  Rangetraitvalues <- diff(Rangetraitvalues)
  M4RTV <- c(M4RTV, Rangetraitvalues)
  M4EDMedianTV <- c(M4EDMedianTV, median(traitvalues))
  
  TBLsimtree <- pd.calc(simtree)
  TBLesttree <- pd.calc(esttree)
  simED <- ed.calc(simtree)
  simED <- simED$spp$ED
  estED <- ed.calc(esttree)
  estED <- estED$spp$ED
  EDproportion <- (estED/simED)
  EDdifference <- (estED - simED)
  M4EDMedianDif <- c(M4EDMedianDif, median(EDdifference))
  M4EDMedianProp <- c(M4EDMedianProp, median(EDproportion))
  
  SpearDif <- cor.test(traitvalues, EDdifference, method = "spearman", alternative = "greater")
  SpearProp <- cor.test(traitvalues, EDproportion, method = "spearman", alternative = "greater")
  rhodif <- SpearDif$estimate
  pdif <- SpearDif$p.value
  rhoprop <- SpearProp$estimate
  pprop <- SpearProp$p.value
  M4RhoDif <- c(M4RhoDif, rhodif)
  M4pDif <- c(M4pDif, pdif)
  M4RhoProp <- c(M4RhoProp, rhoprop)
  M4pProp <- c(M4pProp, pprop)
  print("win")
  setwd("..")
}

setwd("results/ED/Spearman")

M4EDMedianDif <- abs(M4EDMedianDif)
M4EDMedianProp <- abs(M4EDMedianProp)
M4EDMedianTV <- abs(M4EDMedianTV)
M4RhoTrVaEDDif <- cor.test(M4EDMedianTV, M4EDMedianDif, method = "spearman", alternative = "greater")
save(M4RhoTrVaEDDif, file = 'M4Rho_TrVv_EDDif.Rdata')
M4RhoTrVaEDProp <- cor.test(M4EDMedianTV, M4EDMedianProp, method = "spearman", alternative = "greater")
save(M4RhoTrVaEDProp, file = 'M4Rho_TrVv_EDProp.Rdata')
M4RhoRangeEDDif <- cor.test(M4RTV, M4EDMedianDif, method = "spearman", alternative = "greater")
save(M4RhoRangeEDDif, file = 'M4Rho_Range_EDDif.Rdata')
png("M4Rho_Range_EDDif.png")
plot(M4RTV, M4EDMedianDif, xlab = "Trait Value Range", ylab = 'Median ED Error')
dev.off()
M4RhoRangeEDProp <- cor.test(M4RTV, M4EDMedianProp, method = "spearman", alternative = "greater")
save(M4RhoRangeEDProp, file = 'M4Rho_Range_EDProp.Rdata')
png("M4Rho_Range_EDProp.png")
plot(M4RTV, M4EDMedianProp, xlab = "Trait Value Range", ylab = 'Median ED Error')
dev.off()
setwd("../../..")