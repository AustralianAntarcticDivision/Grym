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

## Deterministic initial age structure assuming no fishing and no plus
## group
N0 <- ageStructureD(M*Msf)
N0

## Annual projection from N0
pr <- project(ws,M*Ms,Nref=N0)
## Initial and final abundance
pr$N[c(1,nsteps+1),]

## Deterministic initial age structure with no plus class
ageStructureD(M*Msf)
## Stochastic initial age structure with constant recruitment
ageStructureS(rep(1,20),Msf,M)

## Deterministic initial age structure assuming the plus class
## aggregates 50 additional age classes
ageStructureD(M*Msf,plus=20)
## Stochastic initial age structure with constant recruitment with
## plus class
ageStructureS(rep(1,7+20),Msf,M,plus=TRUE)

## Deterministic initial age structure assuming no upper limit on age
ageStructureD(M*Msf,plus=Inf)

## Scale to recruitment
N0 <- ageStructureD(M*Msf,R=8000)
N0

## Scale age structure to biomass estimated from a five day mid-year
## survey (days 150 to 154)
Bsurvey <- 15000 
N0 <- ageStructureD(M*Msf)
pr <- project(ws,M*Ms,0,0,Nref=N0,nref=1,Bref=Bsurvey,bref=151:155)
N0 <- initial(pr$N)
N0
