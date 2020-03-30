len <- seq(150,450,50)

## Convert length to weight
wgt <- powerLW(len,9.0E-10,3.32)
wgt

## Convert weigth to length
len <- powerWL(wgt,9.0E-10,3.32)
len

## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(0,1,length=nsteps+1)
h <- 1/nsteps

## Ages
ages <- outer(Days,Ages,FUN="+")

## Age-length conversion
ls <- vonBertalanffyAL(ages,t0=0.0667,K=0.5,Linf=500)

## Length-weight conversion
ws <- powerLW(ls,9E-10,3.32)
matplot(Days,ws,type="l",lty=1,xlab="Day",ylab="Weight",main="Weight at Age")
