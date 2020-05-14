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

## Compute model parameters that give mean R = 0.3 and var R = 0.1, 
## where R is the proportion of recruits in the population
ps <- prRecruitParsIB(Msf=Msf, mnR=0.3, vrR=0.01)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileIB(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsIB(10, ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsIB(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsIB(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 1]/rowSums(x=N)
## Estimate mean and variance of R
mean(R)
var(R)

## Histogram of simulated recruitments
rs <- prRecruitsIB(n=100000, ps=ps, mnA=100)
hist(x=rs[rs < 350], breaks=seq(from=0, to=350, by=5), freq=FALSE, 
     xlab="Recruits", main="Inverse Beta")

## Compute model parameters that give mean R = 0.3 and var R = 0.01, 
## where R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older
ps <- prRecruitParsIB(Msf=Msf, mnR=0.3, vrR=0.01, r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileIB(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsIB(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsIB(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsIB(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 3]/rowSums(x=N[, 3:7])
## Estimate mean and variance of R
mean(x=R)
var(x=R)

## Compute model parameters that give mean R = 0.3 and var R = 0.01, 
## where R is the proportion of recruits in the population
ps <- prRecruitParsG(Msf=Msf, mnR=0.3, vrR=0.01)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileG(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsG(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsG(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsG(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 1]/rowSums(x=N)
## Estimate mean and variance of R
mean(x=R)
var(x=R)

## Histogram of simulated recruitments
rs <- prRecruitsG(100000, ps, mnA=100)
hist(x=rs[rs < 350], breaks=seq(from=0, to=350, by=5), freq=FALSE, 
     xlab="Recruits", main="Gamma")

## Compute model parameters that give mean R = 0.3 and var R = 0.01, 
## where R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older
ps <- prRecruitParsG(Msf=Msf, mnR=0.3, vrR=0.01, r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileG(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsG(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsG(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsG(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 3]/rowSums(x=N[, 3:7])
## Estimate mean and variance of R
mean(x=R)
var(x=R)

## Compute model parameters that give mean R = 0.3 and var R = 0.01, 
## where R is the proportion of recruits in the population
ps <- prRecruitParsLN(Msf=Msf, mnR=0.3, vrR=0.01)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileLN(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsLN(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsLN(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsLN(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 1]/rowSums(x=N)
## Estimate mean and variance of R
mean(x=R)
var(x=R)

## Histogram of simulated recruitments
rs <- prRecruitsLN(n=100000, ps=ps, mnA=100)
hist(x=rs[rs < 350], breaks=seq(from=0, to=350, by=5), freq=FALSE, 
     xlab="Recruits", main="Log Normal")

## Compute model parameters that give mean R = 0.3 and var R = 0.01, 
## where R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older
ps <- prRecruitParsLN(Msf=Msf, mnR=0.3, vrR=0.01, r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileLN(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsLN(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsLN(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsLN(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 3]/rowSums(x=N[, 3:7])
## Estimate mean and variance of R
mean(x=R)
var(x=R)

## Compute model parameters that give mean R = 0.3 and var R = 0.2, 
## where R is the proportion of recruits in the population,  assuming
## zero recruitment occurs with probability 0.4
ps <- prRecruitParsDIB(Msf=Msf, mnR=0.3, vrR=0.2, p=0.4)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDIB(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsDIB(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsDIB(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsDIB(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 1]/rowSums(x=N)
## Estimate mean and variance of R
mean(x=R, na.rm=TRUE)
var(x=R, na.rm=TRUE)

## Histogram of simulated recruitments
rs <- prRecruitsDIB(n=100000, ps=ps, mnA=100)
hist(x=rs[rs >0 & rs < 600], breaks=seq(from=0, to=600, by=5), freq=FALSE, 
     xlab="Recruits", main="Delta Inverse Beta")

## Fraction of zero recruitments
mean(x=rs==0)

## Compute model parameters that give mean R = 0.3 and var R = 0.2, 
## where R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older,  assuming
## recruitment occurs with probability 0.4
ps <- prRecruitParsDIB(Msf=Msf, mnR=0.3, vrR=0.2, p=0.4, r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDIB(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsDIB(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsDIB(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsDIB(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 3]/rowSums(x=N[, 3:7])
## Estimate mean and variance of R
mean(x=R, na.rm=TRUE)
var(x=R, na.rm=TRUE)

## Compute model parameters that give mean R = 0.3 and var R = 0.2, 
## where R is the proportion of recruits in the population,  assuming
## zero recruitment occurs with probability 0.4
ps <- prRecruitParsDG(Msf=Msf, mnR=0.3, vrR=0.2, p=0.4)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDG(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsDG(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsDG(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsDG(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 1]/rowSums(x=N)
## Estimate mean and variance of R
mean(x=R, na.rm=TRUE)
var(x=R, na.rm=TRUE)

## Histogram of simulated recruitments
rs <- prRecruitsDG(n=100000, ps=ps, mnA=100)
hist(x=rs[rs >0 & rs < 600], breaks=seq(from=0, to=600, by=5), freq=FALSE, 
     xlab="Recruits", main="Delta Gamma")

## Fraction of zero recruitments
mean(x=rs==0)

## Compute model parameters that give mean R = 0.3 and var R = 0.2, 
## where R is the proportion of individual in age class 3 as a
## fraction of individuals that are that age or older, assuming
## recruitment occurs with probability 0.4
ps <- prRecruitParsDG(Msf=Msf, mnR=0.3, vrR=0.2, p=0.4, r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDG(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsDG(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsDG(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsDG(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 3]/rowSums(x=N[, 3:7])
## Estimate mean and variance of R
mean(x=R, na.rm=TRUE)
var(x=R, na.rm=TRUE)

## Compute model parameters that give mean R = 0.3 and var R = 0.2, 
## where R is the proportion of recruits in the population,  assuming
## zero recruitment events occur with probability 0.4.
ps <- prRecruitParsDLN(Msf=Msf, mnR=0.3, vrR=0.1, p=0.4)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDLN(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsDLN(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsDLN(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsDLN(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 1]/rowSums(x=N)
## Estimate mean and variance of R
mean(x=R, na.rm=TRUE)
var(x=R, na.rm=TRUE)

## Histogram of simulated recruitments
rs <- prRecruitsDLN(n=100000, ps=ps, mnA=100)
hist(x=rs[rs >0 & rs < 600], breaks=seq(from=0, to=600, by=5), freq=FALSE, 
     xlab="Recruits", main="Delta Log Normal")

## Fraction of zero recruitments
mean(x=rs==0)

## Compute model parameters that give mean R = 0.3 and var R = 0.2, 
## where R is the proportion of individual in age class 3 as a
## fraction of individuals that are that age or older,  assuming
## recruitment occurs with probability 0.4
ps <- prRecruitParsDLN(Msf=Msf, mnR=0.3, vrR=0.2, p=0.4, r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDLN(q=0.5, ps=ps, mnA=100)

## Ten random recruit numbers
prRecruitsDLN(n=10, ps=ps, mnA=100)

## Mean number of recruits
mean(x=prRecruitsDLN(n=10000, ps=ps, mnA=100))

## Simulate age profiles and calculate R
N <- matrix(data=0, nrow=10000, ncol=ncol(x=Ms))
for(k in 1:nrow(x=N))
  N[k, ] <- ageStructureS(R=prRecruitsDLN(n=ncol(x=Ms), ps=ps, mnA=100), Msf=Msf, M=ps$M)
R <- N[, 3]/rowSums(x=N[, 3:7])
## Estimate mean and variance of R
mean(x=R, na.rm=TRUE)
var(x=R, na.rm=TRUE)
