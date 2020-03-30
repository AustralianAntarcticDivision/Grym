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
ps <- prRecruitParsGYM(mnR=0.3,vrR=0.01,Msf)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileGYM(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsGYM(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsGYM(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsGYM(10,ps,mnA=100),Msf,ps$M)
R <- N[,1]/rowSums(N)
## Estimate mean and variance of R
mean(R)
var(R)

## Histogram of simulated recruitments
rs <- prRecruitsGYM(100000,ps,mnA=100)
hist(rs[rs < 350],seq(0,350,5),freq=FALSE,xlab="Recruits",main="GYM")

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older
ps <- prRecruitParsGYM(mnR=0.3,vrR=0.01,Msf,r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileGYM(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsGYM(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsGYM(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,100000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsGYM(10,ps,mnA=100),Msf,ps$M)
R <- N[,3]/rowSums(N[,3:7])
## Estimate mean and variance of R
mean(R)
var(R)
