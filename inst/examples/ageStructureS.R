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
ms <- matrix(data=1, nrow=nsteps+1, ncol=length(Ages))
ms <- ms/mean(x=trapz(fs=ms, h=h))
Ms <- ctrapz(fs=ms, h=h)
Msf <- final(P=Ms)
M <- 0.2

## Stochastic initial age structure with log Normal recruitment
ageStructureS(R=rlnorm(n=8, meanlog=log(x=10), sdlog=0.2), Msf=Msf, M=M)

## Plot 100 replicates
Rs <- replicate(n=100, 
                expr=ageStructureS(R=rlnorm(n=8, meanlog=log(x=10), sdlog=0.2), Msf=Msf, M=M))
matplot(x=Rs,type="l", lty=1, col=rgb(0,0,0,0.1), xlab="Age", ylab="Abundance")
lines(x=ageStructureD(MMsf=M*Msf, R=exp(x=log(x=10)+0.1^2/2)), col="firebrick")

Rs <- matrix(data=0, nrow=length(x=Ages), ncol=100)
Rs[,1] <- ageStructureS(R=rlnorm(n=8, meanlog=log(x=10), sdlog=0.2), Msf=Msf, M=M)
for(k in 2:ncol(x=Rs))
  Rs[,k] <- ageStructureS(R=rlnorm(n=8, meanlog=log(x=10), sdlog=0.2),Msf=Msf, M=M, N0=Rs[,k-1])
matplot(x=Rs, type="l", lty=1, col=rgb(0,0,0,0.1), xlab="Age", ylab="Abundance")
lines(x=ageStructureD(MMsf=M*Msf, R=exp(x=log(x=10)+0.1^2/2)), col="firebrick")

