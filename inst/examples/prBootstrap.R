## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(0,1,length=nsteps+1)
h <- 1/nsteps

## Ages
ages <- outer(Days,Ages,FUN="+")
## Age-length and length-weight conversions
ls <- vonBertalanffyAL(ages,t0=0.0667,K=0.5,Linf=500)
ws <- powerLW(ls,9E-10,3.32)

## Constant intra-annual natural mortality
ms <- matrix(1,nsteps+1,length(Ages))
ms <- ms/mean(trapz(ms,h))
Ms <- ctrapz(ms,h)
Msf <- final(Ms)
M <- 0.2

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsIB(mnR=0.3,vrR=0.01,Msf)

## Simulate new mean and variance for 17 surveys
prBootstrap(prRecruitsIB,ps,17,Msf)

## Distribution of 2000 samples
bs <- replicate(2000,unlist(prBootstrap(prRecruitsIB,ps,17,Msf)))
opar <- par(mfrow=c(1,2))
hist(bs[1,],50,xlab="Mean R",main="")
abline(v=median(bs[1,]),col="red")
hist(bs[2,],50,xlab="Variance R",main="")
abline(v=median(bs[2,]),col="red")
par(opar)

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsG(mnR=0.3,vrR=0.01,Msf)

## Simulate new mean and variance for 17 surveys
prBootstrap(prRecruitsG,ps,17,Msf)

## Distribution of 2000 samples
bs <- replicate(2000,unlist(prBootstrap(prRecruitsG,ps,17,Msf)))
opar <- par(mfrow=c(1,2))
hist(bs[1,],50,xlab="Mean R",main="")
abline(v=median(bs[1,]),col="red")
hist(bs[2,],50,xlab="Variance R",main="")
abline(v=median(bs[2,]),col="red")
par(opar)

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsLN(mnR=0.3,vrR=0.01,Msf)

## Simulate new mean and variance for 17 surveys
prBootstrap(prRecruitsLN,ps,17,Msf)

## Distribution of 2000 samples
bs <- replicate(2000,unlist(prBootstrap(prRecruitsLN,ps,17,Msf)))
opar <- par(mfrow=c(1,2))
hist(bs[1,],50,xlab="Mean R",main="")
abline(v=median(bs[1,]),col="red")
hist(bs[2,],50,xlab="Variance R",main="")
abline(v=median(bs[2,]),col="red")
par(opar)

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsGYM(mnR=0.3,vrR=0.01,Msf)

## Simulate new mean and variance for 17 surveys
prBootstrap(prRecruitsGYM,ps,17,Msf)

## Distribution of 2000 samples
bs <- replicate(2000,unlist(prBootstrap(prRecruitsGYM,ps,17,Msf)))
opar <- par(mfrow=c(1,2))
hist(bs[1,],50,xlab="Mean R",main="")
abline(v=median(bs[1,]),col="red")
hist(bs[2,],50,xlab="Variance R",main="")
abline(v=median(bs[2,]),col="red")
par(opar)
