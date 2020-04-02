## Simulate the varianceof R
vs <- replicate(n=10000, expr=resampleRGYM(mnR=0.5, vrR=0.17, n=17)$vrR)
## The simulated variance is often larger than it is possible for the
## variance of a proprtion to be
hist(x=vs, breaks=50, xlab="Variance R", main="Resampled Variances")
