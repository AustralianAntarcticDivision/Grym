## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(from=0, to=1, length=nsteps+1)
h <- 1/nsteps

## Ages
ages <- outer(X=Days, Y=Ages, FUN="+")
## Age-length and length-weight conversions
ls <- vonBertalanffyAL(A=ages, t0=0.0667, K=0.5, Linf=500)
ws <- powerLW(L=ls, a=9E-10, b=3.32)

## Constant intra-annual natural mortality
ms <- matrix(data=1, nrow=nsteps+1, ncol=length(x=Ages))
ms <- ms/mean(x=trapz(fs=ms, h=h))
Ms <- ctrapz(fs=ms, h=h)
Msf <- final(P=Ms)
M <- 0.2

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsIB(mnR=0.3, vrR=0.01, Msf=Msf)

## Simulate new mean and variance for 17 surveys
prBootstrap(prRec=prRecruitsIB, ps=ps, n=17, Msf=Msf)

## Distribution of 2000 samples
bs <- replicate(n=2000, expr=unlist(x=prBootstrap(prRec=prRecruitsIB, ps=ps, n=17, Msf=Msf)))
opar <- par(mfrow=c(1,2))
hist(x=bs[1,], breaks=50, xlab="Mean R", main="")
abline(v=median(bs[1,]), col="red")
hist(x=bs[2,], breaks=50, xlab="Variance R", main="")
abline(v=median(x=bs[2,]), col="red")
par(opar)

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsG(mnR=0.3, vrR=0.01, Msf=Msf)

## Simulate new mean and variance for 17 surveys
prBootstrap(prRec=prRecruitsG, ps=ps, n=17, Msf=Msf)

## Distribution of 2000 samples
bs <- replicate(n=2000, expr=unlist(x=prBootstrap(prRec=prRecruitsG, ps=ps, n=17, Msf=Msf)))
opar <- par(mfrow=c(1,2))
hist(x=bs[1,], breaks=50, xlab="Mean R", main="")
abline(v=median(x=bs[1,]), col="red")
hist(x=bs[2,], breaks=50, xlab="Variance R", main="")
abline(v=median(x=bs[2,]), col="red")
par(opar)

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsLN(mnR=0.3, vrR=0.01, Msf=Msf)

## Simulate new mean and variance for 17 surveys
prBootstrap(prRec=prRecruitsLN, ps=ps, n=17, Msf=Msf)

## Distribution of 2000 samples
bs <- replicate(n=2000, expr=unlist(x=prBootstrap(prRec=prRecruitsLN, ps=ps, n=17, Msf=Msf)))
opar <- par(mfrow=c(1,2))
hist(x=bs[1,], breaks=50, xlab="Mean R", main="")
abline(v=median(x=bs[1,]), col="red")
hist(x=bs[2,], breaks=50, xlab="Variance R", main="")
abline(v=median(x=bs[2,]),col="red")
par(opar)

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsGYM(mnR=0.3, vrR=0.01, Msf=Msf)

## Simulate new mean and variance for 17 surveys
prBootstrap(prRec=prRecruitsGYM, ps=ps, n=17, Msf=Msf)

## Distribution of 2000 samples
bs <- replicate(n=2000,expr=unlist(prBootstrap(prRec=prRecruitsGYM, ps=ps, n=17, Msf=Msf)))
opar <- par(mfrow=c(1,2))
hist(x=bs[1,], breaks=50, xlab="Mean R", main="")
abline(v=median(x=bs[1,]), col="red")
hist(x=bs[2,], breaks=50, xlab="Variance R", main="")
abline(v=median(x=bs[2,]), col="red")
par(opar)
