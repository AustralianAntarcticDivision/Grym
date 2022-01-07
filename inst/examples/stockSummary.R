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
M <- 0.1

## Normalized within year distribution of fishing effort
fwy <- double(length=nsteps+1)
fwy[81:220] <- 1
fwy <- fwy/trapz(fs=fwy, h=h)
## Length based selectivity
ss <- rampOgive(x=ls, x50=420, xrange=30)

## Intra-annual fishing mortality - effort by selectivity
fs <- fwy*ss
Fs <- ctrapz(fs=fs, h=h)
Fsf <- final(P=Fs)
F <- 0.2

## Deterministic initial age structure
N0 <- ageStructureD(MMsf=M*Msf, FFsf=F*Fsf)
N0

## Annual projection from N0
pr <- project(ws=ws, MMs=M*Ms, FFs=F*Fs, Ffs=F*fs, Nref=N0, yield=2)

## Initial abundance
initial(P=pr$N)
## Final abundance
final(P=pr$N)

monitor <- 101:150
## Mean total abundance in the monitoring period
meanStock(P=pr$N, period=monitor)
## Mean total biomass in the monitoring period
meanStock(P=pr$B, period=monitor)

## Maturity at age
gs <- rampOgive(x=ls, x50=410, xrange=40)

## Increments in the spawning period
spawn <- 201:281

## Spawning stock abundance
spawningStock(P=pr$N, gs=gs, period=spawn)
## Spawning stock biomass
spawningStock(P=pr$B, gs=gs, period=spawn)

## Monitoring period within fishing season
monitor <- 101:151
## Exploitable biomass
exploitableStock(P=pr$B, ss=ss, period=monitor)
## Weighted vulnerable biomass
vulnerableStock(P=pr$B, ss=ss, fwy=fwy, period=monitor)

## Monitoring period spans fishing season
monitor <- 41:101
## Exploitable biomass
exploitableStock(P=pr$B, ss=ss, period=monitor)
## Vulnerable biomass weights by effort
vulnerableStock(P=pr$B, ss=ss, fwy=fwy, period=monitor)

## No effort in monitoring period
monitor <- 11:31
## Exploitable biomass
exploitableStock(P=pr$B, ss=ss, period=monitor)
## Vulnerable biomass weights by effort
vulnerableStock(P=pr$B, ss=ss, fwy=fwy, period=monitor)
