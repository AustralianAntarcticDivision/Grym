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
M <- 0.1

## Normalized within year distribution of fishing effort
fwy <- double(nsteps+1)
fwy[81:220] <- 1
fwy <- fwy/trapz(fwy,h)
## Length based selectivity
ss <- rampOgive(ls,420,30)

## Intra-annual fishing mortality - effort by selectivity
fs <- fwy*ss
Fs <- ctrapz(fs,h)
Fsf <- final(Fs)
F <- 0.2

## Deterministic initial age structure
N0 <- ageStructureD(M*Msf,F*Fsf)
N0

## Annual projection from N0
pr <- project(ws,M*Ms,F*Fs,F*fs,Nref=N0,yield=2)

## Initial abundance
initial(pr$N)
## Final abundance
final(pr$N)

## Maturity at age
gs <- rampOgive(ls,410,40)

## Increments in the spawning period
spawn <- 201:281

## Spawning stock abundance
spawningStock(pr$N,gs,spawn)
## Spawning stock biomass
spawningStock(pr$B,gs,spawn)

## Monitoring period within fishing season
monitor <- 101:151
## Exploitable biomass
exploitableStock(pr$B,ss,monitor)
## Weighted vulnerable biomass
vulnerableStock(pr$B,ss,fwy,monitor)

## Monitoring period spans fishing season
monitor <- 41:101
## Exploitable biomass
exploitableStock(pr$B,ss,monitor)
## Vulnerable biomass weights by effort
vulnerableStock(pr$B,ss,fwy,monitor)

## No effort in monitoring period
monitor <- 11:31
## Exploitable biomass
exploitableStock(pr$B,ss,monitor)
## Vulnerable biomass weights by effort
vulnerableStock(pr$B,ss,fwy,monitor)
