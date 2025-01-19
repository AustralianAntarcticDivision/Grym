library(Grym)
library(microbenchmark)
library(ggplot2)

source("./R/Grym.R")

prRecruit <-function(n,mn,vr) {
  q <- rbeta(n,mn*(mn*(1+mn)+vr)/vr,mn*(1+mn)/vr+2)
  q/(1-q)/mn
}

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

##### ctrapz
Fs <- ctrapz(fs=fs, h=h)
Fa <- ctrapz_alt(fs=fs, h=h)
identical(Fs,Fa)

microbenchmark(
  original = ctrapz(fs=fs, h=h),
  alt = ctrapz_alt(fs=fs, h=h),
  times = 10000,
  control = list(warmup = 50)
)




##### ageStructureS

ps <- c(0.7684538, 0.4684104, 1.1680451, 0.2993452)
pr_recs <- prRecruit(length(Msf),ps[3],ps[4])

N0 <- ageStructureS(pr_recs,Msf,M)
Na <- ageStructureS_alt(pr_recs,Msf,M)
identical(N0,Na)

microbenchmark(
  original = ageStructureS(pr_recs,Msf,M),
  alt = ageStructureS_alt(pr_recs,Msf,M),
  times = 10000,
  control = list(warmup = 50)
)

##### project
pr <- project(ws=ws, MMs=M*Ms, FFs=F*Fs, Ffs=F*fs, Nref=N0, yield=2)
pa <- project_alt(ws=ws, MMs=M*Ms, FFs=F*Fs, Ffs=F*fs, Nref=N0, yield=2)
identical(pr,pa)

bnmk <- microbenchmark(
  original = project(ws=ws, MMs=M*Ms, FFs=F*Fs, Ffs=F*fs, Nref=N0, yield=2),
  alt = project_alt(ws=ws, MMs=M*Ms, FFs=F*Fs, Ffs=F*fs, Nref=N0, yield=2),
  times = 10000,
  control = list(warmup = 50)
                )

bnmk

##### rescaleProjection
pq <- rescaleProjection(pr=pr, Nref=N0)
pa <- rescaleProjection_alt(pr=pr, Nref=N0)
identical(pq,pa)

microbenchmark(
  original = rescaleProjection(pr=pr, Nref=N0),
  alt = rescaleProjection_alt(pr=pr, Nref=N0),
  times = 10000,
  control = list(warmup = 50)
)


##### projectC
pr <- projectC(ws=ws, MMs=M*Ms, Fs=Fs, fs=fs, Catch=0.2811238, Nref=N0, yield=2)
pa <- projectC_alt(ws=ws, MMs=M*Ms, Fs=Fs, fs=fs, Catch=0.2811238, Nref=N0, yield=2)
identical(pr,pa)

microbenchmark(
  original = projectC(ws=ws, MMs=M*Ms, Fs=Fs, fs=fs, Catch=0.2811238, Nref=N0, yield=2),
  alt = projectC_alt(ws=ws, MMs=M*Ms, Fs=Fs, fs=fs, Catch=0.2811238, Nref=N0, yield=2),
  Uniroot_alt = projectC_alt2(ws=ws, MMs=M*Ms, Fs=Fs, fs=fs, Catch=0.2811238, Nref=N0, yield=2),
  times = 10000,
  control = list(warmup = 50)
)







