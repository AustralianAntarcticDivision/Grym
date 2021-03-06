% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Grym.R
\name{ageStructureS}
\alias{ageStructureS}
\title{Stochastic Age Structure}
\usage{
ageStructureS(R, Msf, M, Fsf = 0, F = 0, plus = FALSE, N0 = 0)
}
\arguments{
\item{R}{vector of yearly recruitment numbers}

\item{Msf}{vector of final \emph{unscaled} integrated natural mortality}

\item{M}{vector of annual natural mortalities}

\item{Fsf}{vector of final \emph{unscaled} integrated fishing mortality}

\item{F}{vector of annual fishing mortalities}

\item{plus}{whether the final class is a "plus" class}

\item{N0}{an optional initial abundance}
}
\value{
\code{ageStructureS} returns a vector of the abundance of each age class.
}
\description{
Stochastically estimate an initial age structure of the population.
}
\details{
Simulate an initial age structure for a projection by projecting
random recruitment forward in time. It is necessary for the vector
of random recruits \code{R} to be at least as long as there are age
classes with non neglible numbers (including age classes
aggregated into the plus group).

The annual scalings of natural or fishing mortality may vary
inter-annually, and so the arguments \code{M} and \code{F} should be vectors
of length 1 representing a constant scaling, or of \code{length(R)-1}
representing scalings that vary inter-annually.
}
\examples{
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

}
\seealso{
\code{\link[=ageStructureD]{ageStructureD()}} for deterministic initial age structure,
\code{\link[=project]{project()}}, \code{\link[=projectC]{projectC()}} for annual projection
}
