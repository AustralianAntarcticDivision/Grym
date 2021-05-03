age <- 2:8

## Convert age to length
len <- vonBertalanffyAL(A=age, t0=0.0667, K=0.5, Linf=500)
len

## Convert length to age
age <- vonBertalanffyLA(L=len, t0=0.0667, K=0.5, Linf=500)
age

## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(from=0, to=1, length=nsteps+1)
h <- 1/nsteps

## Ages
ages <- outer(X=Days, Y=Ages, FUN="+")

## Age-length conversion
ls <- vonBertalanffyAL(A=ages, t0=0.0667, K=0.5, Linf=500)
matplot(x=Days, y=ls, type="l", lty=1, xlab="Day", ylab="Length", main="Length at Age")

## Age-length conversion - growth occurs in middle of the year
ls <- vonBertalanffyRAL(A=ages, t0=0.0667, K=0.5, Linf=500, f0=0.3, f1=0.7)
matplot(x=Days, y=ls, type="l", lty=1, xlab="Day", ylab="Length", main="Length at Age")

## Age-length conversion - growth occurs across the year boundary
ls <- vonBertalanffyRAL(A=ages, t0=0.0667, K=0.5, Linf=500, f0=0.8, f1=1.2)
matplot(x=Days, y=ls, type="l", lty=1, xlab="Day", ylab="Length", main="Length at Age")
