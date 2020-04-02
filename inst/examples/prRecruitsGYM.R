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
ps <- prRecruitParsGYM(Msf=Msf, mnR=0.3, vrR=0.01)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileGYM(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsGYM(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsGYM(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsGYM(n=10, ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 1]/rowSums(x=N)
## Estimate mean and variance of R
mean(x=R)
var(x=R)

## Histogram of simulated recruitments
rs <- prRecruitsGYM(n=100000, ps=ps, mnA=100)
hist(x=rs[rs < 350], seq(from=0, to=350, by=5), freq=FALSE, xlab="Recruits", main="GYM")

## Compute model parameters that give mean R = 0.3 and var R = 0.01, 
## when R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older
ps <- prRecruitParsGYM(mnR=0.3, vrR=0.01, Msf=Msf, r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileGYM(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsGYM(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsGYM(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=100000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsGYM(n=10, ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 3]/rowSums(x=N[, 3:7])
## Estimate mean and variance of R
mean(x=R)
var(x=R)
