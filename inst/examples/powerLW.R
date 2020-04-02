len <- seq(from=150, to=450, by=50)

## Convert length to weight
wgt <- powerLW(L=len, a=9.0E-10, b=3.32)
wgt

## Convert weigth to length
len <- powerWL(W=wgt, a=9.0E-10, b=3.32)
len

## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(from=0, to=1, length=nsteps+1)
h <- 1/nsteps

## Ages
ages <- outer(X=Days, Y=Ages, FUN="+")

## Age-length conversion
ls <- vonBertalanffyAL(A=ages, t0=0.0667, K=0.5, Linf=500)

## Length-weight conversion
ws <- powerLW(L=ls, a=9E-10, b=3.32)
matplot(x=Days, y=ws, type="l",lty=1,xlab="Day",ylab="Weight",main="Weight at Age")
