filename <- read.csv('Settings.csv')
filename <- filename$runName
simlist <- (list())
estlist <- (list())
folders1 <- grep("Model1", dir(), value = T)
folders2 <- grep("Model2", dir(), value = T)
folders3 <- grep("Model3", dir(), value = T)
folders4 <- grep("Model4", dir(), value = T)
folders5 <- grep("Model5", dir(), value = T)
M1steps <- vector()
M1subs <- vector()
M2steps <- vector()
M2subs <- vector()
M3steps <- vector()
M3subs <- vector()
M4steps <- vector()
M4subs <- vector()
n <- 0

for(i in 1:length(folders1)){
  setwd (folders1[i])
  load(paste0(folders1[i], ".Rdata"))
  steps <- tree.sim$infotab[, 1]
  subrates <- tree.sim$infotab[, 2]
  stepsbysubs <- data.frame(steps, subrates)
  stepsbysubs <- stepsbysubs
  M1steps <- c(M1steps, steps)
  M1subs <- c(M1subs, subrates)
  setwd("..")
  n <- n+1
  print(n)
}
M1stepsbysubs <- data.frame(M1steps, M1subs)
M1stepsbysubs <- M1stepsbysubs[order(M1steps),]
M1stepsbysubs <- M1stepsbysubs[-(1:100),]


for(i in 1:length(folders2)){
  setwd (folders2[i])
  load(paste0(folders2[i], ".Rdata"))
  steps <- tree.sim$infotab[, 1]
  subrates <- tree.sim$infotab[, 2]
  stepsbysubs <- data.frame(steps, subrates)
  stepsbysubs <- stepsbysubs
  M2steps <- c(M2steps, steps)
  M2subs <- c(M2subs, subrates)
  setwd("..")
  n <- n+1
  print(n)
}
M2stepsbysubs <- data.frame(M2steps, M2subs)
M2stepsbysubs <- M2stepsbysubs[order(M2steps),]
M2stepsbysubs <- M2stepsbysubs[-(1:100),]

for(i in 1:length(folders3)){
  setwd (folders3[i])
  load(paste0(folders3[i], ".Rdata"))
  steps <- tree.sim$infotab[, 1]
  subrates <- tree.sim$infotab[, 2]
  stepsbysubs <- data.frame(steps, subrates)
  stepsbysubs <- stepsbysubs
  M3steps <- c(M3steps, steps)
  M3subs <- c(M3subs, subrates)
  setwd("..")
  n <- n+1
  print(n)
}
M3stepsbysubs <- data.frame(M3steps, M3subs)
M3stepsbysubs <- M3stepsbysubs[order(M3steps),]
M3stepsbysubs <- M3stepsbysubs[-(1:100),]

for(i in 1:length(folders4)){
  setwd (folders4[i])
  load(paste0(folders4[i], ".Rdata"))
  steps <- tree.sim$infotab[, 1]
  subrates <- tree.sim$infotab[, 2]
  stepsbysubs <- data.frame(steps, subrates)
  stepsbysubs <- stepsbysubs
  M4steps <- c(M4steps, steps)
  M4subs <- c(M4subs, subrates)
  setwd("..")
  n <- n+1
  print(n)
}
M4stepsbysubs <- data.frame(M4steps, M4subs)
M4stepsbysubs <- M4stepsbysubs[order(M4steps),]
M4stepsbysubs <- M4stepsbysubs[-(1:100),]

mean(M1stepsbysubs$M1subs)
mean(M2stepsbysubs$M2subs)
mean(M3stepsbysubs$M3subs)
mean(M4stepsbysubs$M4subs)
sd(M1stepsbysubs$M1subs)
sd(M2stepsbysubs$M2subs)
sd(M3stepsbysubs$M3subs)
sd(M4stepsbysubs$M4subs)
range(M1stepsbysubs$M1subs)
range(M2stepsbysubs$M2subs)
range(M3stepsbysubs$M3subs)
range(M4stepsbysubs$M4subs)