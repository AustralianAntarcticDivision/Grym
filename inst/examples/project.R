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

## Final abundance
pr$N[nsteps+1,]
## Total annual yield
sum(pr$Y[nsteps+1,])

## Plot initial abundance, projected abundance, biomass and yield
opar <- par(mfrow=c(2,2), mar=c(5,4,2,2)+0.1)
pal <- hcl(h=seq(from=15, to=375, length=8)[1:7], l=65, c=100)
plot(x=Ages, y=N0, xlab="Age", ylab="Abundance")
matplot(x=Days, y=pr$N, xlab="Day", ylab="Abundance", type="l", lty=1, col=pal)
matplot(x=Days, y=pr$B, xlab="Day", ylab="Biomass", type="l", lty=1, col=pal)
matplot(x=Days, y=pr$Y, xlab="Day", ylab="Yield", type="l", lty=1, col=pal)
par(opar)

## Annual projection given relative abundance and total
## biomass estimated from a five day mid-year survey (days 150 to 154)
Nsurvey <- c(1.0,0.67,0.46,0.33,0.25,0.18,0.13) 
Bsurvey <- 15000 
pr <- project(ws=ws, MMs=M*Ms, FFs=F*Fs, Ffs=F*fs, Nref=Nsurvey, 
              nref=151:155, Bref=Bsurvey, yield=2)

## Total annual yield
sum(pr$Y[nsteps+1,])

## Plot initial abundance, projected abundance, biomass and yield
opar <- par(mfrow=c(2,2), mar=c(5,4,2,2)+0.1)
pal <- hcl(h=seq(from=15, to=375, length=8)[1:7],l=65,c=100)
plot(x=Ages, y=pr$N[1,], xlab="Age", ylab="Abundance")
matplot(x=Days, y=pr$N, xlab="Day", ylab="Abundance", type="l" ,lty=1, col=pal)
matplot(x=Days, y=pr$B, xlab="Day", ylab="Biomass", type="l", lty=1, col=pal)
matplot(x=Days, y=pr$Y, xlab="Day", ylab="Yield", type="l", lty=1, col=pal)
par(opar)

## Ten year projection from virgin stock collating final annual
## abundance, biomass and yield
Years <- 1:10
Nf <- Bf <- Yf <- matrix(data=NA, nrow=length(Years), ncol= length(Ages))

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
opar <- par(mfrow=c(2,2), mar=c(5,4,2,2)+0.1)
pal <- hcl(h=seq(from=15, to=375, length=8)[1:7], l=65, c=100)
matplot(x=Years, y=Nf, xlab="Year", ylab="Abundance", type="l", lty=1, col=pal)
matplot(x=Years, y=Bf, xlab="Year", ylab="Biomass", type="l", lty=1, col=pal)
matplot(x=Years, y=Yf, xlab="Year", ylab="Yield", type="l", lty=1, col=pal)
par(opar) 
