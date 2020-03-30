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

## Ten year projection from virgin stock collating final annual
## abundance, biomass and yield
Years <- 1:10
Nf <- Bf <- Yf <- matrix(NA,length(Years),length(Ages))

## Deterministic initial age structure assuming no fishing
N0 <- ageStructureD(M*Msf,R=8000)

## Project first year
pr <- project(ws,M*Ms,F*Fs,F*fs,Nref=N0,yield=1)
N0 <- advance(pr$N,R=8000)
Nf[1,] <- final(pr$N)
Bf[1,] <- final(pr$B)
Yf[1,] <- pr$Y

## Compute projections for subsequent years by rescaling
for(k in 2:nrow(Nf)) {
  pr <- rescaleProjection(pr,Nref=N0)
  N0 <- advance(pr$N,R=8000)
  Nf[k,] <- final(pr$N)
  Bf[k,] <- final(pr$B)
  Yf[k,] <- pr$Y
}

## Plot annual abundance, biomass and yield
opar <- par(mfrow=c(2,2),mar=c(5,4,2,2)+0.1)
pal <- hcl(h=seq(15,375,length=8)[1:7],l=65,c=100)
matplot(Years,Nf,xlab="Year",ylab="Abundance",type="l",lty=1,col=pal)
matplot(Years,Bf,xlab="Year",ylab="Biomass",type="l",lty=1,col=pal)
matplot(Years,Yf,xlab="Year",ylab="Yield",type="l",lty=1,col=pal)
par(opar)
