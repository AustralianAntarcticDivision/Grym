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

## Compute model parameters that give mean R = 0.3 and var R = 0.1,
## where R is the proportion of recruits in the population
ps <- prRecruitParsIB(Msf,mnR=0.3,vrR=0.01)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileIB(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsIB(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsIB(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsIB(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,1]/rowSums(N)
## Estimate mean and variance of R
mean(R)
var(R)

## Histogram of simulated recruitments
rs <- prRecruitsIB(100000,ps,mnA=100)
hist(rs[rs < 350],seq(0,350,5),freq=FALSE,xlab="Recruits",main="Inverse Beta")

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## where R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older
ps <- prRecruitParsIB(Msf,mnR=0.3,vrR=0.01,r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileIB(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsIB(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsIB(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsIB(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,3]/rowSums(N[,3:7])
## Estimate mean and variance of R
mean(R)
var(R)

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## where R is the proportion of recruits in the population
ps <- prRecruitParsG(Msf,mnR=0.3,vrR=0.01)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileG(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsG(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsG(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsG(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,1]/rowSums(N)
## Estimate mean and variance of R
mean(R)
var(R)

## Histogram of simulated recruitments
rs <- prRecruitsG(100000,ps,mnA=100)
hist(rs[rs < 350],seq(0,350,5),freq=FALSE,xlab="Recruits",main="Gamma")

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## where R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older
ps <- prRecruitParsG(Msf,mnR=0.3,vrR=0.01,r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileG(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsG(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsG(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsG(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,3]/rowSums(N[,3:7])
## Estimate mean and variance of R
mean(R)
var(R)

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## where R is the proportion of recruits in the population
ps <- prRecruitParsLN(Msf,mnR=0.3,vrR=0.01)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileLN(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsLN(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsLN(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsLN(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,1]/rowSums(N)
## Estimate mean and variance of R
mean(R)
var(R)

## Histogram of simulated recruitments
rs <- prRecruitsLN(100000,ps,mnA=100)
hist(rs[rs < 350],seq(0,350,5),freq=FALSE,xlab="Recruits",main="Log Normal")

## Compute model parameters that give mean R = 0.3 and var R = 0.01,
## where R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older
ps <- prRecruitParsLN(Msf,mnR=0.3,vrR=0.01,r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileLN(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsLN(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsLN(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsLN(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,3]/rowSums(N[,3:7])
## Estimate mean and variance of R
mean(R)
var(R)

## Compute model parameters that give mean R = 0.3 and var R = 0.2,
## where R is the proportion of recruits in the population, assuming
## zero recruitment occurs with probability 0.4
ps <- prRecruitParsDIB(Msf,mnR=0.3,vrR=0.2,p=0.4)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDIB(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsDIB(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsDIB(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsDIB(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,1]/rowSums(N)
## Estimate mean and variance of R
mean(R,na.rm=TRUE)
var(R,na.rm=TRUE)

## Histogram of simulated recruitments
rs <- prRecruitsDIB(100000,ps,mnA=100)
hist(rs[rs >0 & rs < 600],seq(0,600,5),freq=FALSE,xlab="Recruits",main="Delta Inverse Beta")

## Fraction of zero recruitments
mean(rs==0)

## Compute model parameters that give mean R = 0.3 and var R = 0.2,
## where R is the proportion of individual in age class 3 as a fraction
## of individuals that are that age or older, assuming
## recruitment occurs with probability 0.4
ps <- prRecruitParsDIB(Msf,mnR=0.3,vrR=0.2,p=0.4,r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDIB(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsDIB(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsDIB(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsDIB(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,3]/rowSums(N[,3:7])
## Estimate mean and variance of R
mean(R,na.rm=TRUE)
var(R,na.rm=TRUE)

## Compute model parameters that give mean R = 0.3 and var R = 0.2,
## where R is the proportion of recruits in the population, assuming
## zero recruitment occurs with probability 0.4
ps <- prRecruitParsDG(Msf,mnR=0.3,vrR=0.2,p=0.4)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDG(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsDG(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsDG(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsDG(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,1]/rowSums(N)
## Estimate mean and variance of R
mean(R,na.rm=TRUE)
var(R,na.rm=TRUE)

## Histogram of simulated recruitments
rs <- prRecruitsDG(100000,ps,mnA=100)
hist(rs[rs >0 & rs < 600],seq(0,600,5),freq=FALSE,xlab="Recruits",main="Delta Gamma")

## Fraction of zero recruitments
mean(rs==0)

## Compute model parameters that give mean R = 0.3 and var R = 0.2,
## where R is the proportion of individual in age class 3 as a
## fraction of individuals that are that age or older, assuming
## recruitment occurs with probability 0.4
ps <- prRecruitParsDG(Msf,mnR=0.3,vrR=0.2,p=0.4,r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDG(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsDG(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsDG(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsDG(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,3]/rowSums(N[,3:7])
## Estimate mean and variance of R
mean(R,na.rm=TRUE)
var(R,na.rm=TRUE)

## Compute model parameters that give mean R = 0.3 and var R = 0.2,
## where R is the proportion of recruits in the population, assuming
## zero recruitment events occur with probability 0.4.
ps <- prRecruitParsDLN(Msf,mnR=0.3,vrR=0.1,p=0.4)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDLN(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsDLN(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsDLN(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsDLN(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,1]/rowSums(N)
## Estimate mean and variance of R
mean(R,na.rm=TRUE)
var(R,na.rm=TRUE)

## Histogram of simulated recruitments
rs <- prRecruitsDLN(100000,ps,mnA=100)
hist(rs[rs >0 & rs < 600],seq(0,600,5),freq=FALSE,xlab="Recruits",main="Delta Log Normal")

## Fraction of zero recruitments
mean(rs==0)

## Compute model parameters that give mean R = 0.3 and var R = 0.2,
## where R is the proportion of individual in age class 3 as a
## fraction of individuals that are that age or older, assuming
## recruitment occurs with probability 0.4
ps <- prRecruitParsDLN(Msf,mnR=0.3,vrR=0.2,p=0.4,r=3)
## Estimated natural mortality scaling
ps$M

## Median number of recruits
prRecruitsQuantileDLN(0.5,ps,mnA=100)

## Ten random recruit numbers
prRecruitsDLN(10,ps,mnA=100)

## Mean number of recruits
mean(prRecruitsDLN(10000,ps,mnA=100))

## Simulate age profiles and calculate R
N <- matrix(0,10000,ncol(Ms))
for(k in 1:nrow(N))
  N[k,] <- ageStructureS(prRecruitsDLN(ncol(Ms),ps,mnA=100),Msf,ps$M)
R <- N[,3]/rowSums(N[,3:7])
## Estimate mean and variance of R
mean(R,na.rm=TRUE)
var(R,na.rm=TRUE)
