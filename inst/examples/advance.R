## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(from=0,to=1,length=nsteps+1)
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

## Define fishing mortalities - mid-year fishing with length based selectivity 
fs <- double(length=nsteps+1)
fs[81:220] <- 1
fs <- fs/trapz(fs=fs, h=h)
fs <- fs*rampOgive(x=ls, x50=420, xrange=30)
Fs <- ctrapz(fs=fs, h=h)
Fsf <- final(P=Fs)
F <- 0.2

## Deterministic initial age structure
N0 <- ageStructureD(MMsf=M*Msf, FFsf=F*Fsf)
N0

## Annual projection from N0
pr <- project(ws=ws, MMs=M*Ms, FFs=F*Fs, Ffs=F*fs, Nref=N0)

## Final abundance
pr$N[nsteps+1,]

## Initial abundance for next projection assuming no plus group
advance(N=pr$N, R=1, plus=FALSE)

## Initial abundance for next projection assuming plus group
advance(N=pr$N, R=1, plus=TRUE)

## Ten year projection from virgin stock collating final annual
## abundance, biomass and yield
Years <- 1:10
Nf <- Bf <- Yf <- matrix(data=NA, nrow=length(x=Years), ncol=length(x=Ages)) #produce identical matrices Nf, Bf and Yf

## Deterministic initial age structure assuming no fishing and 8000
## new recruits each year
N0 <- ageStructureD(MMsf=M*Msf, R=8000)

## Project forward with constant fishing mortality
for(k in 1:nrow(x=Nf)) {
  pr <- project(ws=ws, MMs=M*Ms, FFs=F*Fs, Ffs=F*fs, Nref=N0, yield=1)
  N0 <- advance(N=pr$N, R=8000)
  Nf[k,] <- final(P=pr$N)
  Bf[k,] <- final(P=pr$B)
  Yf[k,] <- pr$Y
}

## Plot annual abundance, biomass and yield
opar <- par(mfrow=c(2,2),mar=c(5,4,2,2)+0.1)
pal <- hcl(h=seq(from=15, to=375, length=8)[1:7], l=65, c=100)
matplot(x=Years, y=Nf, xlab="Year", ylab="Abundance", type="l", lty=1, col=pal)
matplot(x=Years, y=Bf, xlab="Year", ylab="Biomass", type="l", lty=1, col=pal)
matplot(x=Years, y=Yf, xlab="Year", ylab="Yield", type="l", lty=1, col=pal)
par(opar)

## Ten year projection from virgin stock collating final annual
## abundance, biomass and yield
Years <- 1:10
Nf <- Bf <- Yf <- matrix(data=NA, nrow=length(x=Years), ncol=length(x=Ages))

## Deterministic initial age structure assuming no fishing, 8000
## new recruits each year and a plus class
N0 <- ageStructureD(MMsf=M*Msf, R=8000, plus=50)

## Project forward with constant fishing mortality assuming plus class
for(k in 1:nrow(x=Nf)) {
  pr <- project(ws=ws, MMs=M*Ms, FFs=F*Fs, Ffs=F*fs, Nref=N0, yield=1)
  N0 <- advance(N=pr$N, R=8000, plus=TRUE)
  Nf[k,] <- final(P=pr$N)
  Bf[k,] <- final(P=pr$B)
  Yf[k,] <- pr$Y
}

## Plot annual abundance, biomass and yield
opar <- par(mfrow=c(2,2), mar=c(5,4,2,2)+0.1)
pal <- hcl(h=seq(from=15, to=375, length=8)[1:7], l=65, c=100)
matplot(x=Years, y=Nf, xlab="Year", ylab="Abundance", type="l", lty=1, col=pal)
matplot(x=Years, y=Bf, xlab="Year", ylab="Biomass", type="l", lty=1, col=pal)
matplot(x=Years, y=Yf, xlab="Year", ylab="Yield", type="l", lty=1, col=pal)
par(opar)
