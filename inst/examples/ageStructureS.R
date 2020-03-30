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

## Stochastic initial age structure with log Normal recruitment
ageStructureS(rlnorm(8,log(10),0.2),Msf,M)

## Plot 100 replicates
Rs <- replicate(100,ageStructureS(rlnorm(8,log(10),0.2),Msf,M))
matplot(Rs,type="l",lty=1,col=rgb(0,0,0,0.1),xlab="Age",ylab="Abundance")
lines(ageStructureD(M*Msf,R=exp(log(10)+0.1^2/2)),col="firebrick")

Rs <- matrix(0,length(Ages),100)
Rs[,1] <- ageStructureS(rlnorm(8,log(10),0.2),Msf,M)
for(k in 2:ncol(Rs))
  Rs[,k] <- ageStructureS(rlnorm(8,log(10),0.2),Msf,M,N0=Rs[,k-1])
matplot(Rs,type="l",lty=1,col=rgb(0,0,0,0.1),xlab="Age",ylab="Abundance")
lines(ageStructureD(M*Msf,R=exp(log(10)+0.1^2/2)),col="firebrick")
