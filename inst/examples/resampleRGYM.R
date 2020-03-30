## Simulate the varianceof R
vs <- replicate(10000,resampleRGYM(mnR=0.5,vrR=0.17,n=17)$vrR)
## The simulated variance is often larger than it is possible for the
## variance of a proprtion to be
hist(vs,50,xlab="Variance R",main="Resampled Variances")
