age <- 2:8

## Convert age to length
len <- vonBertalanffyAL(age,t0=0.0667,K=0.5,Linf=500)
len

## Convert length to age
age <- vonBertalanffyLA(len,t0=0.0667,K=0.5,Linf=500)
age

## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(0,1,length=nsteps+1)
h <- 1/nsteps

## Ages
ages <- outer(Days,Ages,FUN="+")

## Age-length conversion
ls <- vonBertalanffyAL(ages,t0=0.0667,K=0.5,Linf=500)
matplot(Days,ls,type="l",lty=1,xlab="Day",ylab="Length",main="Length at Age")

## Age-length conversion - growth occurs in middle of the year
ls <- vonBertalanffyRAL(ages,t0=0.0667,K=0.5,Linf=500,f0=0.3,f1=0.7)
matplot(Days,ls,type="l",lty=1,xlab="Day",ylab="Length",main="Length at Age")
