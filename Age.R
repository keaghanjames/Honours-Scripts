age <- vector()
tips <- vector()
M1age <- vector()
M2age <- vector()
M3age <- vector()
M4age <- vector()
M1tips <- vector()
M2tips <- vector()
M3tips <- vector()
M4tips <- vector()

filename <- read.csv('Settings.csv')
filename <- filename$runName
simlist <- (list())
estlist <- (list())
folders1 <- grep("Model1", dir(), value = T)
folders2 <- grep("Model2", dir(), value = T)
folders3 <- grep("Model3", dir(), value = T)
folders4 <- grep("Model4", dir(), value = T)
folders5 <- grep("Model5", dir(), value = T)

for(i in 1: length(folders1)){
  #print(getwd())
  setwd(folders1[i])
  #print(getwd())
  load(paste0(folders1[i], ".Rdata"))
  x <- tail(tree.sim$infotab)
  x <- x[6,1]
  age <- c(age, x)
  M1age <- c(M1age, x)
  simtree <- read.tree(paste0(folders1[i], "_sim.tre"))
  y <- length(simtree$tip.label)
  M1tips <- c(M1tips, y)
  setwd("..")
  print(":)")
}

for(i in 1: length(folders2)){
  #print(getwd())
  setwd(folders2[i])
  #print(getwd())
  load(paste0(folders2[i], ".Rdata"))
  x <- tail(tree.sim$infotab)
  x <- x[6,1]
  age <- c(age, x)
  M2age <- c(M2age, x)
  simtree<- read.tree(paste0(folders2[i], "_sim.tre"))
  y <- length(simtree$tip.label)
  M2tips <- c(M2tips, y)
  setwd("..")
  print(":)")}

for(i in 1: length(folders3)){
  #print(getwd())
  setwd(folders3[i])
  #print(getwd())
  load(paste0(folders3[i], ".Rdata"))
  x <- tail(tree.sim$infotab)
  x <- x[6,1]
  age <- c(age, x)
  M3age <- c(M3age, x)
  simtree<- read.tree(paste0(folders3[i], "_sim.tre"))
  y <- length(simtree$tip.label)
  M3tips <- c(M3tips, y)
  setwd("..")
  print(":)")}

for(i in 1: length(folders4)){
  #print(getwd())
  setwd(folders4[i])
  #print(getwd())
  load(paste0(folders4[i], ".Rdata"))
  x <- tail(tree.sim$infotab)
  x <- x[6,1]
  age <- c(age, x)
  M4age <- c(M4age, x)
  simtree<- read.tree(paste0(folders4[i], "_sim.tre"))
  y <- length(simtree$tip.label)
  M4tips <- c(M4tips, y)
  setwd("..")
  print(":)")}
tips <- c(M1tips, M2tips, M3tips, M4tips)
median(age)