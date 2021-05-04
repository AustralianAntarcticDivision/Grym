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
fwy1 <- fwy2 <- fwy
## Fleet 1 doesn't fish from day 150 to 219
fwy1[151:220] <- 0
## Fleet 2 doesn't fish from day 80 to 149
fwy2[81:150] <- 0

## Length based selectivity
ss <- rampOgive(x=ls, x50=420, xrange=30)

## Intra-annual fishing mortality - effort by selectivity
fs <- array(0,c(nsteps+1,length(Ages),2))
fs[,,1] <- fwy1*ss
fs[,,2] <- fwy2*ss
Fs <- ctrapz(fs=fs, h=h)
Fsf <- Fs[nrow(Fs),,] #final(P=Fs)
F <- c(0.2,0.2)

## Deterministic initial age structure
N0 <- ageStructureD(MMsf=M*Msf, FFsf=Fsf%*%F)
N0

## Scale by F
FFs <- scale3(Fs,F,sum=TRUE)
Ffs <- scale3(fs,F,sum=FALSE)
## Annual projection from N0
pr <- projectMF(ws=ws, MMs=M*Ms, FFs=FFs, Ffs=Ffs, Nref=N0, yield=1)
pr$Y

colSums(pr$Y)
