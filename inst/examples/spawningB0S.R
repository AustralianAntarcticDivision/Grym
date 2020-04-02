## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(from=0,to=1,length=nsteps+1)
h <- 1/nsteps

## Ages
ages <- outer(X=Days,Y=Ages,FUN="+")
## Age-length and length-weight conversions
ls <- vonBertalanffyAL(A=ages,t0=0.0667,K=0.5,Linf=500)
ws <- powerLW(L=ls, a=9E-10, b=3.32)

## Constant intra-annual natural mortality
ms <- matrix(data=1,nrow=nsteps+1,ncol=length(Ages))
ms <- ms/mean(x=trapz(fs=ms,h=h))
Ms <- ctrapz(fs=ms,h=h)
Msf <- final(P=Ms)

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## when R is the proportion of recruits in the population
ps <- prRecruitParsIB(Msf=Msf,mnR=0.2,vrR=0.01)
## Estimated natural mortality scaling
ps$M

R <- matrix(data=prRecruitsIB(n=10000,ps=ps),nrow=1000,ncol=10)

## Maturity at age
gs <- rampOgive(x=ls,x50=410,xrange=40)

## Increments in the spawning period
spawn <- 201:281

## Virgin spawning stock abundance
spawningB0S(R=R,gs=gs,ws=ws,Ms=Ms,M=ps$M,spawn=spawn)
