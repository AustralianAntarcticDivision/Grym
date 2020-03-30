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

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsIB(mnR=0.2,vrR=0.01,Msf)
## Estimated natural mortality scaling
ps$M

R <- matrix(prRecruitsIB(10000,ps),1000,10)

## Maturity at age
gs <- rampOgive(ls,410,40)

## Increments in the spawning period
spawn <- 201:281

## Virgin spawning stock abundance
spawningB0S(R,ws,gs,Ms,ps$M,spawn)
