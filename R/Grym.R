##' The von Bertalanffy growth curve and its inverse.
##'
##' These functions compute the von Bertalanffy growth curve using the
##' standard parameterization implemented in the GYM.
##'
##' The functions `vonBertalanffyAL` and `vonBertalanffyLA` use a von
##' Bertalanffy growth to map age to length and length to age. The
##' variants `vonBertalanffyRAL` and `vonBertalanffyRLA` assume that
##' growth is restricted to some fraction of each year.
##'
##' @title von Bertalanffy Growth
##' @param A vector of ages
##' @param L vector of lengths
##' @param t0 reference time parameter
##' @param K growth rate parameter
##' @param Linf asymptotic length parameter
##' @param f0 fraction of year at which growth starts
##' @param f1 fraction of year at which growth ends
##' @return `vonBertalanffyAL` returns the vector of lengths
##'   corresponding to the ages `A`, and `vonBertalanffy` returns a
##'   vector of ages corresponding to the lengths `L`.
##' @seealso [powerLW()] for length-weight modelling.
##' @example inst/examples/vonBertalanffyAL.R
##' @export
## ---- vonBertalanffyAL
vonBertalanffyAL <- function(A,t0,K,Linf) {
  Linf*(1-exp(-K*(A-t0)))
}
## ----

##' @rdname vonBertalanffyAL
##' @export
## ---- vonBertalanffyLA
vonBertalanffyLA <- function(L,t0,K,Linf) {
  t0 + log(Linf/(Linf-L))/K
}
## ----


##' @rdname vonBertalanffyAL
##' @export
## ---- vonBertalanffyRAL
vonBertalanffyRAL <- function(A,t0,K,Linf,f0=0,f1=1) {
  ## Remap interval of growth to the full year
  A <- (A-f0)%/%1 + pmin(1,pmax(0,((A-f0)%%1)/(f1-f0)))
  A[A<0] <- 0
  vonBertalanffyAL(A,t0,K,Linf)
}
## ----

##' @rdname vonBertalanffyAL
##' @export
## ---- vonBertalanffyRLA
vonBertalanffyRLA <- function(L,t0,K,Linf,f0=0,f1=1) {
  A <- vonBertalanffyLA(L,t0,K,Linf)
  A%/%1 +(f1-f0)*A%%1+f0
}
## ----



##' Power law length-weight relationship and its inverse.
##'
##' These functions compute power law length weight relationship
##' using the same parameterization implemented in the GYM.
##'
##' @title Power law length weight
##' @param L a vector or matrix of lengths
##' @param W a vector or matrix of weights
##' @param a scale parameter
##' @param b power parameter
##' @return `powerLW` returns a vector or matrix of weights
##'   corresponding to the lengths `L` and `powerLW` returns a vector
##'   or matrix of lengths corresponding to the weights `W`.
##' @seealso [vonBertalanffyAL()] for age-length modelling.
##' @example inst/examples/powerLW.R
##' @export
## ---- powerLW
powerLW <- function(L,a,b) {
  a*L^b
}
## ----

##' @rdname powerLW
##' @export
## ---- powerWL
powerWL <- function(W,a,b) {
  (W/a)^(1/b)
}
## ----



##' Ramp shaped ogive function.
##'
##' Compute a ramp shaped ogive function, parameterized so that `x50`
##' is the middle of the ramp, and `xrange` its width.
##' @title Ramp ogive
##' @param x a matrix or vector
##' @param x50 mid-point of ramp
##' @param xrange width of ramp
##' @return A matrix or vector of probabilities of the same size as x
##' @example inst/examples/rampOgive.R
##' @importFrom stats punif
##' @export
## ---- rampOgive
rampOgive <- function(x,x50,xrange) {
  punif(x,x50-xrange/2,x50+xrange/2)
}
## ----


##' Logistic shaped ogive function.
##'
##' Compute a logistic shaped ogive function, parameterized so that
##' `x50` is the middle of the ramp, and `d95` is the difference in
##' the 95th and 50th percentiles.
##' 
##' @title Logistic Ogive
##' @param x a matrix or vector
##' @param x50 50th percentile
##' @param d95 difference between the 95th and 50th percentiles
##' @return A matrix or vector of probabilities of the same size as x
##' @example inst/examples/logisticOgive.R
##' @importFrom stats plogis qlogis
##' @export
## ---- logisticOgive
logisticOgive <- function(x,x50,d95) {
  plogis(x,x50,d95/qlogis(0.95))
}
## ----


##' Double Normal Selectivity.
##'
##' Compute a double Normal or "dome shaped" selectivity function,
##' parameterized in terms of the endpoint, width and minima of the
##' ascending and descending limbs.
##'
##' @title Double Normal Selectivity
##' @param x a matrix or vector
##' @param x1 endpoint of the ascending limb
##' @param x2 endpoint of the descending limb
##' @param s1 width of the ascending limb
##' @param s2 width of the descending limb
##' @param y1 minimum selectivity the ascending limb
##' @param y2 dispersion of the descending limb
##' @return A matrix or vector of probabilities of the same size as x
##' @example inst/examples/doubleNormalSelectivity.R
##' @export
## ---- doubleNormalSelectivity
doubleNormalSelectivity <- function(x,x1,x2,s1,s2,y1,y2) {
  ifelse(x<x1,
         y1+(1-y1)*exp(-((x-x1)/s1)^2/2),
         ifelse(x>x2,
                y2+(1-y2)*exp(-((x-x2)/s2)^2/2),1))
}
## ----


##' Numerical quadrature by the composite trapezoidal rule.
##'
##' Given a vector `fs` of function values evaluated over an even
##' partition of an interval, `trapz` approximates the definite
##' integral over the interval by the composite [trapezoidal
##' rule](wikipedia.org/wiki/Trapezoidal_rule), while `ctrapz`
##' approximates the cumulative integral over successive subintervals.
##
##' If `fs` is a matrix then each column is assumed to correspond to a
##' different function and an approximation is computed for each
##' column.
##'
##' @title Trapezoidal Quadrature
##' @param fs a vector or matrix of function values
##' @param h the width of each subinterval
##' @return
##' * `trap` returns the approximation to the integral of the function
##' or functions represented by `fs`.
##' * `ctrapz` returns a vector or matrix of the same size as `fs`
##' representing the approximate cumulative integrals of the function
##' or functions represented by `fs`.
##' @example inst/examples/trapz.R
##' @export
## ---- trapz
trapz <- function(fs,h=1) {
  if(is.matrix(fs))
    h*(colSums(fs)-(fs[1L,]+fs[nrow(fs),])/2)
  else
    h*(sum(fs)-(fs[1L]+fs[length(fs)])/2)
}
## ----


##' @rdname trapz
##' @export
## ---- ctrapz
ctrapz <- function(fs,h=1) {
  if(is.matrix(fs)) {
    Fs <- matrix(0,nrow(fs),ncol(fs))
    for(j in seq_len(ncol(fs)))
      Fs[-1L,j] <- cumsum((h/2)*(fs[-1L,j]+fs[-nrow(fs),j]))
  } else {
    Fs <- c(0,cumsum((h/2)*(fs[-1L]+fs[-length(fs)])))
  }
  Fs
}
## ----


##' Compute averages over a time interval, but the interval represents
##' a single point in time, the value at that time point is returned.
##'
##' Given a vector `fs` of function values evaluated over an even
##' partition of an interval, `trapzMean` approximates the average
##' over the interval by the composite [trapezoidal
##' rule](wikipedia.org/wiki/Trapezoidal_rule).
##'
##  If `fs` is a ' matrix then each column is assumed to correspond to
##a different ' function and the average is computed columnwise.
##'
##' In each case, if `fs` does not span a time interval but represents
##' a single point in time, the value at that time point is returned.
##'
##' @title Average over a time interval.
##' @param fs a vector or matrix of function values
##' @return
##' Returns the approximation to the average of the function or functions
##' represented by `fs`.
##' @example inst/examples/trapzMeans.R
##' @export
## ---- trapzMeans
trapzMeans <- function(fs) {
  if(NROW(fs)==1L) fs else trapz(fs)/(NROW(fs)-1L)
}
## ----



##' Project the abundance, biomass and yield in each age class forward
##' over one year when natural and fishing mortalities are known.
##'
##' The `project` function calculates the annual projection of
##' abundance `N`, biomass `B` and yield `Y` that reproduces the
##' prescribed reference abundance or biomass given known natural and
##' fishing mortalities.  This function requires the natural and
##' fishing mortality be prescribed in the form of the *scaled*
##' integrated natural and fishing mortalities and the *scaled*
##' instantaneous fishing mortality.
##'
##' When no reference biomass is specified, the vector `Nref` of
##' reference abundances specifies the mean abundance of each age
##' class for the time steps specified by `nref` (`nref=1` corresponds
##' to specifying an initial abundance).  When a reference biomass
##' `Bref` is specified, `Nref` is interpreted as relative abundances
##' and the projections are scaled to the reference biomass.
##'
##' The matrices `ws`, `MMs`, `FFs`, `Ffs`, `N`, `B` and `F` are all
##' of equal size, and have one column for each age class, and one row
##' for each time step through the year.  So if a daily time step is
##' adopted, each matrix will have 366 rows with the first
##' representing the midnight on the last day of the previous year,
##' and the last representing midnight on the last day of the current
##' year.
##'
##' If `yield=0`, then yield is not calculated.  If `yield=1`, then
##' only the final yield from each age class is calculated and `Y` is
##' a vector with one element for each age class.  If `yield=2` then
##' `Y` is the yield per time step and age class and is a matrix the
##' same size as `N` and `B`.
##'
##' @title Yearly projections for known fishing mortality.
##' @param ws matrix of weight at age
##' @param MMs matrix of *scaled* integrated natural mortalities
##' @param FFs matrix of *scaled* integrated fishing mortalities
##' @param Ffs matrix of *scaled* fishing mortalities
##' @param Nref vector of reference abundances for each age class
##' @param nref vector of time steps at which the reference abundance
##'     is estimated
##' @param Bref optional reference biomass
##' @param bref vector of time steps at which the reference biomass is
##'     estimated
##' @param yield should yield be calculated for the final time step,
##'     ever time step, or not at all (see details)
##' @return `project` returns a list with elements
##' * `N` matrix of abundance
##' * `B` matrix of biomass
##' * `Y` matrix or vector of yield
##' @seealso [projectC()], [advance()] for annual projection,
##'   [ageStructureD()], [ageStructureS()] for initial age structure
##' @example inst/examples/project.R
##' @export
## ---- project
project <- function(ws,MMs,FFs=0,Ffs=0,Nref=1,nref=1,Bref=NA,bref=nref,yield=2L) {

  if(length(Nref)==1L) Nref <- rep.int(Nref,ncol(MMs))

  ## Integrate N and scale to reference abundance
  N <- exp(-MMs-FFs)
  for(j in seq_len(ncol(N))) N[,j] <- Nref[j]/mean(N[nref,j])*N[,j]

  ## Scale by weight at age
  B <- ws*N

  ## Rescale to match reference biomass
  if(!is.na(Bref)) {
    r <- Bref/sum(trapzMeans(B[bref,,drop=FALSE]))
    N <- r*N
    B <- r*B
  }

  ## Integrate yield
  Y <- switch(yield,
              trapz(Ffs*B,1/(nrow(B)-1L)),
              ctrapz(Ffs*B,1/(nrow(B)-1L)))

  list(N=N,B=B,Y=Y)
}
## ----

##' Rescale a projection to match a reference abundance and biomass.
##'
##' Rescale an existing annual projection of abundance `N`, biomass
##' `B` and yield `Y` to reproduce the prescribed reference abundance
##' or biomass.
##'
##' When no reference biomass is specified, the vector `Nref` of
##' reference abundances specifies the mean abundance of each age
##' class for the time steps specified by `nref` (`nref=1` corresponds
##' to specifying an initial abundance).  When a reference biomass
##' `Bref` is specified, `Nref` is interpreted as relative abundances
##' and the projections are scaled to the reference biomass.
##'
##' @title Rescale projection
##' @param pr a projection
##' @param Nref vector of reference abundances for each age class
##' @param nref vector of time steps at which the reference abundance is estimated
##' @param Bref optional reference biomass
##' @param bref vector of time steps at which the reference biomass is estimated
##' @return `project` returns a list with elements
##' * `N` matrix of abundance
##' * `B` matrix of biomass
##' * `Y` matrix or vector of yield
##' * `F` the fishing mortality from `pr`
##' @seealso [project()], [advance()] for annual projection,
##'   [ageStructureD()], [ageStructureS()] for initial age structure
##' @example inst/examples/rescaleProjection.R
##' @export
## ---- rescaleProjection
rescaleProjection <- function(pr,Nref=1,nref=1,Bref=NA,bref=nref) {

  N <- pr$N; B <- pr$B; Y <- pr$Y
  if(length(Nref)==1L) Nref <- rep.int(Nref,ncol(N))

  ## Rescale to match reference abundance
  for(j in seq_len(ncol(N))) {
    r <- Nref[j]/mean(pr$N[nref,j])
    N[,j] <- r*N[,j]
    B[,j] <- r*B[,j]
    if(!is.null(Y))
      if(is.matrix(Y)) Y[,j] <- r*Y[,j] else Y[j] <- r*Y[j]
  }

  ## Rescale to match reference biomass
  if(!is.na(Bref)) {
    r <- Bref/sum(trapzMeans(B[bref,,drop=FALSE]))
    N <- r*N
    B <- r*B
    if(!is.null(Y)) Y <- r*Y
  }

  list(N=N,B=B,Y=Y,F=pr$F)
}
## ----


##' Projection of abundance, biomass and yield by age class forward
##' over a year when natural mortality and total annual catch are known.
##'
##' The `projectC` function calculates the annual projection of
##' abundance `N`, biomass `B` and yield `Y` that reproduces the
##' prescribed reference abundance or biomass given known natural
##' mortality and total catch. This function requires the *scaled*
##' integrated natural mortality, and the *unscaled* integrated and
##' instantaneous fishing mortalities, and it repeatedly calls
##' `project` to determine the annual scaling of fishing mortality `F`
##' that reproduces the target catch.
##'
##' If the target catch cannot be recovered with an `F` smaller than
##' `Fmax`, then the projection is performed with `F` capped at
##' `Fmax`, and it is the responsibility of the user to check if the
##' target catch has been recovered.  As the target catch is only ever
##' approximately recovered it is advised to test for `F==Fmax`.
##'
##' When no reference biomass is specified, the vector `Nref` of
##' reference abundances specifies the mean abundance of each age
##' class for the time steps specified by `nref` (`nref=1` corresponds
##' to specifying an initial abundance).  When a reference biomass
##' `Bref` is specified, `Nref` is interpreted as relative abundances
##' and the projections are scaled to the reference biomass.
##'
##' The matrices `ws`, `MMs`, `Fs`, `fs`, `N`, `B` and `F` are all of
##' equal size, and have one column for each age class, and one row
##' for each time step through the year.  So if a daily time step is
##' adopted, each matrix will have 366 rows with the first
##' representing the midnight on the last day of the previous year,
##' and the last representing midnight on the last day of the current
##' year.
##'
##' If `yield=0`, then yield is not calculated.  If `yield=1`, then
##' only the final yield from each age class is calculated and `Y` is
##' a vector with one element for each age class.  If `yield=2` then
##' `Y` is the yield per time step and age class and is a matrix the
##' same size as `N` and `B`.
##'
##' @title Yearly projections for known catch.
##' @param ws matrix of weight at age
##' @param MMs matrix of *scaled* integrated natural mortalities
##' @param Fs matrix of *unscaled* integrated fishing mortalities
##' @param fs matrix of *unscaled* fishing mortalities
##' @param Catch target total annual catch
##' @param Nref vector of reference abundances for each age class
##' @param nref vector of time steps at which the reference abundance is estimated
##' @param Bref optional reference biomass
##' @param bref vector of time steps at which the reference biomass is estimated
##' @param yield should yield be calculated for the final time step,
##'     ever time step, or not at all (see details)
##' @param Fmax the maximum reasonable annual scaling of fishing mortality
##' @param tol the tolerance on target catch
##' @return `projectC` returns a list with elements
##' * `N` matrix of abundance
##' * `B` matrix of biomass
##' * `Y` vector or matrix of yield
##' * `F` annual fishing mortality required to reproduce the observed catch.
##' @seealso [project()], [advance()] for annual projection,
##'   [ageStructureD()], [ageStructureS()] for initial age structure
##' @example inst/examples/projectC.R
##' @importFrom stats uniroot
##' @export
## ---- projectC
projectC <- function(ws,MMs,Fs,fs,Catch,Nref,nref=1,Bref=NA,bref=nref,yield=2L,
                     Fmax=2.5,tol=1.0E-6) {

  ## Function to calculate error in modelled yield
  err <- function(F) {
    sum(project(ws,MMs,F*Fs,F*fs,Nref,nref,Bref,bref,yield=1)$Y) - Catch
  }

  if(Catch > 0) {
    errmax <- err(Fmax)
    if(errmax < 0) {
      ## Catch is unattainable
      F <- Fmax
    } else {
      ## Solve err(F) = 0
      sol <- uniroot(err,lower=0,upper=Fmax,f.lower=-Catch,f.upper=errmax,tol=tol)
      F <- sol$root
    }
  } else {
    F <- 0
  }

  ## Recalculate
  c(project(ws,MMs,F*Fs,F*fs,Nref,nref,Bref,bref,yield),F=F)
}
## ----


##' Given projected abundances, construct the initial abundances for
##' the next projection.
##'
##' At the end of each projection, the final cohort numbers advance by
##' one age class to become the initial numbers for the next
##' projection. If the final class is a "plus" class, the numbers from
##' the previous year are retained, otherwise the final age class is
##' discarded.
##'
##' @title Advance age classes.
##' @param N projected abundance
##' @param R number of recruits
##' @param plus whether the final age class is a plus class
##' @return A vector of initial abundances for the next year.
##' @seealso [project()], [projectC()] for annual projection
##' @export
## ---- advance
advance <- function(N,R=0,plus=FALSE) {
  ## Extract final numbers and advance cohorts
  N1 <- N[nrow(N),]
  N0 <- c(R,N1[-length(N1)])
  ## Adjust plus group
  if(plus) N0[length(N0)] <- N0[length(N0)] + N1[length(N0)]
  N0
}
## ----



##' Deterministically compute the equilibrium age structure of the population.
##'
##' The equilibrium age structure is determined assuming the number
##' of recruits and the natural and fishing mortalities do not vary
##' inter-annually.
##'
##' @title Deterministic Age Structure
##' @param MMsf vector of *scaled* final integrated natural mortality
##' @param FFsf vector of *scaled* final integrated fishing mortality
##' @param R the mean numbers of recruits
##' @param plus the number of additional age classes in the plus class
##'   (`plus=0` implies no plus class).
##' @return `ageStructureD` returns a vector of the abundance of each age class.
##' @seealso [ageStructureS()] for stochastic initial age structure,
##'   [project()], [projectC()] for annual projection
##' @example inst/examples/ageStructureD.R
##' @export
## ---- ageStructureD
ageStructureD <- function(MMsf,FFsf=0,R=1,plus=0) {
  ## Number of age classes
  n <- length(MMsf)
  ## Single year survival for each age class
  S <- exp(-MMsf-FFsf)
  ## Survival of plus class
  Sp <- S[n]
  ## Cumulative survival
  S <- cumprod(c(1,S[-n]))
  ## Adjust for plus class
  if(plus) S[n] <- S[n]*(1-Sp^(1+plus))/(1-Sp)
  R*S
}
## ----


##' Stochastically estimate an initial age structure of the population.
##'
##' Simulate an initial age structure for a projection by projecting
##' random recruitment forward in time. It is necessary for the vector
##' of random recruits `R` to be at least as long as there are age
##' classes with non neglible numbers (including age classes
##' aggregated into the plus group).
##'
##' The annual scalings of natural or fishing mortality may vary
##' inter-annually, and so the arguments `M` and `F` should be vectors
##' of length 1 representing a constant scaling, or of `length(R)-1`
##' representing scalings that vary inter-annually.
##'
##' @title Stochastic Age Structure
##' @param R vector of yearly recruitment numbers
##' @param Msf vector of final *unscaled* integrated natural mortality
##' @param M vector of annual natural mortalities
##' @param Fsf vector of final *unscaled* integrated fishing mortality
##' @param F vector of annual fishing mortalities
##' @param plus whether the final class is a "plus" class
##' @param N0 an optional initial abundance
##' @return `ageStructureS` returns a vector of the abundance of each age class.
##' @seealso [ageStructureD()] for deterministic initial age structure,
##'   [project()], [projectC()] for annual projection
##' @example inst/examples/ageStructureS.R
##' @export
## ---- ageStructureS
ageStructureS <- function(R,Msf,M,Fsf=0,F=0,plus=FALSE,N0=0) {

  n <- length(Msf)
  ## Initialize abundance
  N <- rep(N0,length.out=n)
  N[1L] <- R[1L]

  ## Optimize for constant mortalities
  if(length(M)==1L && length(F)==1L) {

    ## Precompute survivals
    S <- exp(-M*Msf-F*Fsf)
    for(k in seq_len(length(R)-1L)) {
      ## Project forward one year
      N <- N*S
      ## Advance age classes and recruit
      Np <- N[n]
      N <- c(R[k+1L],N[-n])
      if(plus) N[n] <- N[n]+Np
    }

  } else {

    if(length(M)==1L) M <- rep.int(M,length(R)-1L)
    if(length(F)==1L) F <- rep.int(F,length(R)-1L)

    for(k in seq_len(length(R)-1L)) {
      ## Project forward one year
      N <- N*exp(-M[k]*Msf-F[k]*Fsf)
      ## Advance age classes and recruit
      Np <- N[n]
      N <- c(R[k+1L],N[-n])
      if(plus) N[n] <- N[n]+Np
    }
  }
  N
}
## ----


##' Calculate the stock summaries for a given monitroing period.
##'
##' The `spawningStock` and `exploitableStock` functions are just
##' special cases of `meanStock`.  `vulnerableStock` differs in that
##' it also weights by the within year fishing effort.
##'
##' @title Stock summaries
##' @param gs matrix of fraction mature at age
##' @param ss matrix of selectivity at age
##' @param ws matrix of weights
##' @param P projected abundance or biomass
##' @param period the time steps corresponding to the
##'   spawning/monitoring period
##' @param fwy vector of within year distribution of effort for the
##'   year.
##' @return
##'
##' `initial` and `final` return the first and last rows the yearly
##' projected abundance or biomass.
##'
##' `meanStock` returns a weighted mean abundance or biomass in a
##' given period.
##'
##' `spawningStock` returns the mean abundance or biomass of mature
##' fish in a given time period.
##'
##' `exploitableStock` returns the mean abundance or biomass that
##' could be exploited by fishing in a given time period.
##'
##' `vulnerableStock` returns the effort weighted mean exploitable
##' abundance or biomass in a given time period.
##'
##' @example inst/examples/stockSummary.R
##' @rdname stockSummary
##' @export
## ---- initial
initial <- function(P) P[1L,]
## ----

##' @rdname stockSummary
##' @export
## ---- final
final <- function(P) P[nrow(P),]
## ----

##' @rdname stockSummary
##' @export
## ---- meanStock
meanStock <- function(P,ws=1,period) {
  if(identical(ws,1))
    sum(trapzMeans(P[period,,drop=FALSE]))
  else
    sum(trapzMeans(ws[period,,drop=FALSE]*P[period,,drop=FALSE]))
}
## ----

##' @rdname stockSummary
##' @export
## ---- spawningStock
spawningStock <- function(P,gs,period) {
  sum(trapzMeans((gs*P)[period,,drop=FALSE]))
}
## ----

##' @rdname stockSummary
##' @export
## ---- exploitableStock
exploitableStock <- function(P,ss,period) {
  sum(trapzMeans((ss*P)[period,,drop=FALSE]))
}
## ----

##' @rdname stockSummary
##' @export
## ---- vulnerableStock
vulnerableStock <- function(P,ss,fwy,period) {
  m <- mean(fwy[period])
  if(m==0) 0 else sum(trapzMeans((fwy*ss*P)[period,,drop=FALSE]))/m
}
## ----



##' Deterministically estimate virgin spawning stock biomass.
##'
##' Estimate virgin spawning stock biomass by deterministically
##' computing an initial age structure and projecting forward in time
##' to compute spawning stock biomass.
##'
##' @title Deterministic Spawning Biomass
##' @param gs a matrix of maturity at age
##' @param ws matrix of weight at age
##' @param MMs vector of *unscaled* integrated natural mortality
##' @param spawn the increments in which spawning occurs.
##' @param R the mean numbers of recruits
##' @param plus the number of additional age classes in the plus class
##'   (`plus=0` implies no plus class).
##' @return The estimated spawning stock biomass
##' @example inst/examples/spawningB0D.R
##' @export
## ---- spawningB0D
spawningB0D <- function(gs,ws,MMs,spawn,R=1,plus=0) {
  MMsf <- MMs[nrow(MMs),]

  ## Deterministically compute the initial age structure , and project
  ## forward to determine the spawning biomass
  N0 <- ageStructureD(MMsf,R=R,plus=plus)
  pr <- project(ws,MMs,Nref=N0)
  spawningStock(pr$B,gs,spawn)
}
## ----


##' Stochastically estimate virgin spawning stock biomass.
##'
##' Estimate virgin spawning stock biomass by repeatedly simulating an
##' initial age structure and projecting forward in time to compute
##' spawning stock biomass.
##'
##' The matrix of random numbers of recruits `R` is used to simulate
##' initial age structures.  One spawning biomass is simulated for
##' every row of `R`, and `R` must have at least as many columns as
##' there are age classes with non neglible numbers (including age
##' classes aggregated into the plus group).
##'
##' The annual scaling of natural mortality may vary inter-annually,
##' in which case `M` should be a matrix of the same size as `R`.
##'
##' @title Stochastic Median Spawning Biomass
##' @param R a matrix of yearly recruitment numbers
##' @param gs a matrix of maturity at age
##' @param ws matrix of weight at age
##' @param Ms vector of *unscaled* integrated natural mortality
##' @param M matrix of annual natural mortalities
##' @param spawn the increments in which spawning occurs.
##' @param plus whether the final class is a "plus" class
##' @return The median, mean and standard deviation of the spawning
##'     stock biomass estimated through simulation, and the number of
##'     simulations used.
##' @example inst/examples/spawningB0S.R
##' @importFrom stats median sd
##' @export
## ---- spawningB0S
spawningB0S <- function(R,gs,ws,Ms,M,spawn,plus=FALSE) {

  Msf <- Ms[nrow(Ms),]
  SB0 <- double(nrow(R))
  gs <- gs[spawn,,drop=FALSE]

  if(length(M)==1L) {
    ## Constant M - project once and compute spawning biomass by
    ## rescaling by random age structures.  This could be optimized further.
    pr <- project(ws,M*Ms)
    bs <- trapzMeans(gs*pr$B[spawn,,drop=FALSE])
    S <- exp(-M*Msf)
    S <- cumprod(c(1,S[-length(Msf)],rep(if(plus) S[length(Msf)] else 0,ncol(R)-length(Msf))))
    for(k in 1:nrow(R)) {
      ## Simulate random age structure
      N0 <- rev(R[k,])*S
      if(plus) N0[length(Msf)] <- sum(N0[length(Msf):length(N0)])
      SB0[k] <- sum(N0[1:length(Msf)]*bs)
    }

  } else {
    ## Varying M - simulate random age structure and project to compute
    ## spawning biomass
    for(k in 1:nrow(R)) {
      N0 <- ageStructureS(R[k,],Msf,M[k,],plus=plus)
      pr <- project(ws,M[k,ncol(R)]*Ms,Nref=N0)
      SB0[k] <- sum(trapzMeans(gs*pr$B[spawn,,drop=FALSE]))
    }
  }

  list(median=median(SB0),mean=mean(SB0),sd=sd(SB0),n=nrow(R))
}
## ----



##' Generate random recruit numbers according to the proportional
##' recruitment model as implemented in the GYM.
##'
##' These functions reproduce the proportional recruits model as
##' implemented in the GYM and described in de la Mare (1994).  The
##' implementation in the GYM differs slightly from that described in de
##' la Mare (1994) in that a different bias correction is used.
##'
##' The proportional recruit model defines R, the number of
##' individuals in age class r as a proportion of the population of
##' the same age or older. Given estimates of the mean and variance of
##' R from independent surveys, the model determines the scaling of
##' natural mortality `M` and random recruitment numbers that will
##' reproduce the observed mean and variance of R.  The model
##' generates random recruitment numbers as a stochastic adds of the
##' expected size of the population above recruitment age and the GYM
##' assumes the corresponding proportion is Beta distributed.
##'
##' The `prRecruitParGYM` function generates a list of parameters
##' required by `prRecruitsGYM`.
##'
##' @title Proportional recruitment (GYM)
##' @param n number of recruitment events
##' @param q vector of quantiles
##' @param ps a list of parameters returned by `prRecruitParGYM`
##' @param mnA mean number of recruits
##' @param Msf matrix of unscaled final integrated natural mortalities
##' @param mnR the mean of the proportion R
##' @param vrR the variance of the proportion R
##' @param p (ignored)
##' @param r the reference age class
##' @param plus the number of additional age classes in the plus class
##'   (`plus=0` implies no plus class).
##' @param Mbrak the values of natural mortality to assess when
##'   bracketing a root for the nonlinear solver
##' @param tol the tolerance on mean R.
##' @return The function `prRecruitParGYM` returns a list with
##'   elements
##'
##' * `a`, `b` parameters of the Beta distribution
##' * `s0` summed survivorship (T in de la Mare (94))
##' * `B` bias correction
##' * `M` the estimated scaling of natural mortality
##'
##' The function `prRecruitsQuantileGYM` returns a vector of
##' quantiles from the distribution of recruit numbers.
##'
##' The function `prRecruitsGYM` returns a vector of
##' random numbers of recruits.
##'
##' @references
##' de la Mare, WK. (1994) "Modelling Krill Recruitment." *CCAMLR Science* 1: 49-54.
##' @example inst/examples/prRecruitsGYM.R
##' @importFrom stats rbeta
##' @export
## ---- prRecruitsGYM
prRecruitsGYM <- function(n,ps,mnA=1) {
  R <- rbeta(n,ps$a,ps$b)
  Q <- (R/(1-R)-ps$B)
  ## Reject negative or Inf deviates
  fail <- Q<0 | !is.finite(Q)
  while(any(fail)) {
    R <- rbeta(sum(fail),ps$a,ps$b)
    Q[fail] <- (R/(1-R)-ps$B)
    fail <- Q<0 | !is.finite(Q)
  }
  (ps$s0*mnA)*Q
}
## ----


##' @rdname prRecruitsGYM
##' @importFrom stats qbeta
##' @export
## ---- prRecruitsQuantileGYM
prRecruitsQuantileGYM <- function(q,ps,mnA=1) {
  r50 <- qbeta(q,ps$a,ps$b)
  (ps$s0*mnA)*(r50/(1-r50)-ps$B)
}
## ----


##' @rdname prRecruitsGYM
##' @importFrom stats uniroot
##' @export
## ---- prRecruitsParsGYM
prRecruitParsGYM <- function(Msf,mnR,vrR,r=1,p=0,plus=0,Mbrak=seq(0,2.5,0.1),tol=1.0E-8) {

  err <- function(M) {
    ## One year survivals for each age class
    S <- exp(-M*Msf)
    ## Cumulative survivals from r adjusted for plus class
    S <- cumprod(c(1,S[r:(length(S)-1L)],rep(S[length(S)],plus)))

    mnR+vrR-1/sum(S)
  }

  ## Solve err(M) == 0
  ys <- double(length(Mbrak))
  ys[1L] <- err(Mbrak[1L])
  for(k in 2L:length(Mbrak)) {
    ys[k] <- err(Mbrak[k])
    if(ys[k-1L]*ys[k] <= 0) break
  }
  sol <- uniroot(err,Mbrak[(k-1L):k],tol=tol)
  ## Raise error if uniroot fails?

  ## Set M and calculate sums
  M <- sol$root
  ## One year survivals for each age class
  S <- exp(-M*Msf)
  ## Cumulative survivals from 0 adjusted for plus class
  S0 <- cumprod(c(S[-length(S)],rep(S[length(S)],plus)))
  s0 <- sum(S0)
  ## Cumulative survivals from r adjusted for plus class
  Sr <- cumprod(c(S[r:(length(S)-1L)],rep(S[length(S)],plus)))
  s1 <- sum(Sr)
  s2 <- sum(Sr^2)

  deflate <- s1^2/(s1^2+s2)
  vrR <- deflate*vrR

  ## Beta shape parameters for deflated variance
  a <- mnR*(mnR*(1-mnR)-vrR)/vrR
  b <- (1-mnR)*(mnR*(1-mnR)-vrR)/vrR

  ## Bias correction from dlM(94)
  ## B <- vrR/(1-mnR)^3

  ## Bias correction (as implemented in GYM)
  B <- mnR*vrR*(1/(1-mnR)^2+mnR/(1-mnR)^3)

  list(a=a,b=b,s0=s0,B=B,M=M)
}
## ----



##' Generate random recruit numbers according to the proportional
##' recruitment model in which the distribution of the fraction of
##' recruits is inverse Beta (IB), Gamma (G), log Normal (LN) or a
##' zero inflated variant of these (DIB,DG,DLN).
##'
##' The proportional recruit model defines R, the number of
##' individuals in age class r as a proportion of the population of
##' the same age or older. Given estimates of the mean and variance of
##' R from independent surveys, the model determines the scaling of
##' natural mortality `M` and random recruitment numbers that will
##' reproduce the observed mean and variance of R.  The model
##' generates random recruitment numbers as a stochastic fraction of
##' the expected size of the population above recruitment age; these
##' functions generate random numbers of recruits when that fraction
##' is inverse Beta (IB), Gamma (G) or log Normal (LN), or a delta
##' variant (DIB, DG or DLN) where zero recruitment occurs with
##' probability `p`.
##'
##' `prRecruitPars` computes the scaling of natural mortality and mean
##' and variance of the stochastic fraction required to reproduce the
##' mean and variance observed in R.
##'
##' The functions `prRecruitParsIB`, `prRecruitsQuantileIB` and
##' `prRecruitsIB` compute the model parameters, quantiles of the
##' distribution of recruits, and random numbers of recruits when the
##' stochastic fraction is inverse Beta.
##'
##' The functions `prRecruitParsG`, `prRecruitsQuantileG` and
##' `prRecruitsG` compute the model parameters, quantiles of the
##' distribution of recruits, and random numbers of recruits when the
##' stochastic fraction is Gamma.
##'
##' The functions `prRecruitParsLN`, `prRecruitsQuantileLN` and
##' `prRecruitsLN` compute the model parameters, quantiles of the
##' distribution of recruits, and random numbers of recruits when the
##' stochastic fraction is log Normal.
##'
##' The functions `prRecruitParsDIB`, `prRecruitsQuantileDIB` and
##' `prRecruitsDIB` compute the model parameters, quantiles of the
##' distribution of recruits, and random numbers of recruits when the
##' stochastic fraction is delta inverse Beta.
##'
##' The functions `prRecruitParsDG`, `prRecruitsQuantileDG` and
##' `prRecruitsDG` compute the model parameters, quantiles of the
##' distribution of recruits, and random numbers of recruits when the
##' stochastic fraction is delta Gamma.
##'
##' The functions `prRecruitParsDLN`, `prRecruitsQuantileDLN` and
##' `prRecruitsDLN` compute the model parameters, quantiles of the
##' distribution of recruits, and random numbers of recruits when the
##' stochastic fraction is delta log Normal.
##'
##' @title Proportional recruitment
##' @param n number of recruitment events
##' @param q vector of quantiles
##' @param ps a list of model parameters
##' @param mnA mean number of recruits
##' @param Msf matrix of final unscaled integrated natural mortalities
##' @param mnR the mean of the proportion R
##' @param vrR the variance of the proportion R
##' @param p the probability of a zero recruitment (ignored by
##'   non-delta variants).
##' @param r the reference age class
##' @param plus the number of additional age classes in the plus class
##'   (`plus=0` implies no plus class).
##' @param Mbrak the values of natural mortality to assess when
##'   bracketing a root for the nonlinear solver
##' @param tol the tolerance on mean R.
##' @return the family of `prRecruitPar` functions each return a list
##'   with elements
##'
##' * `IB`, `G`, `LN`, `DIB`, `DG`, `DLN` parameters of the inverse beta,
##'   gamma or log normal or a delta variant.
##' * `mnQ`, `vrQ` mean and variance of the stochastic fraction
##' * `M` the estimated scaling of natural mortality
##' * `s0`, `s1`, `s2` summed survivorships
##'
##' The family of `prRecruitsQuantile` functions returns a vector of
##' quantiles from the distribution of recruit numbers.
##'
##' The family of `prRecruits` functions each return a vector of
##' random numbers of recruits.
##'
##' @note These routines must fail if `mnR` is smaller than the
##'   reciporical of the number of age classes of the same age or
##'   older than the reference age class.
##' @example inst/examples/prRecruits.R
##' @rdname prRecruits
##' @importFrom stats uniroot
##' @export
## ---- prRecruitPars0
prRecruitPars0 <- function(Msf,mnR,vrR,r=1,plus=0,Mbrak=seq(0,2.5,0.1),tol=1.0E-10) {

  err <- function(M) {
    ## One year survivals for each age class
    S <- exp(-M*Msf)
    ## Cumulative survivals from r adjusted for plus class
    Sr <- cumprod(c(S[r:(length(S)-1L)],rep(S[length(S)],plus)))
    s1 <- sum(Sr)
    s2 <- sum(Sr^2)
    1/(1+s1)+(s2-s1)*(1+s1)/(s1^2+s2)*vrR-mnR
  }

  ## Solve err(M) == 0
  ys <- double(length(Mbrak))
  ys[1L] <- err(Mbrak[1L])
  for(k in 2L:length(Mbrak)) {
    ys[k] <- err(Mbrak[k])
    if(ys[k-1L]*ys[k] <= 0) break
  }
  sol <- uniroot(err,Mbrak[(k-1L):k],tol=tol)
  ## Raise error if uniroot fails?

  ## Set M and calculate sums
  M <- sol$root
  ## One year survivals for each age class
  S <- exp(-M*Msf)
  ## Cumulative survivals from 0 adjusted for plus class
  S0 <- cumprod(c(S[-length(S)],rep(S[length(S)],plus)))
  s0 <- sum(S0)
  ## Cumulative survivals from r adjusted for plus class
  Sr <- cumprod(c(S[r:(length(S)-1L)],rep(S[length(S)],plus)))
  s1 <- sum(Sr)
  s2 <- sum(Sr^2)
  mnQ <- 1/s0
  vrQ <- (1+s1)^4/(s0^2*(s1^2+s2))*vrR

  list(mnQ=mnQ,vrQ=vrQ,M=M,s0=s0,s1=s1,s2=s2)
}
## ----


##' @rdname prRecruits
##' @importFrom stats rbeta
##' @export
## ---- prRecruitsIB
prRecruitsIB <- function(n,ps,mnA=1) {
  R <- rbeta(n,ps$IB[1L],ps$IB[2L])
  (ps$s0*mnA)*(R/(1-R))
}
## ----


##' @rdname prRecruits
##' @importFrom stats qbeta
##' @export
## ---- prRecruitsQuantileIB
prRecruitsQuantileIB <- function(q,ps,mnA=1) {
  r50 <- qbeta(q,ps$IB[1L],ps$IB[2L])
  (ps$s0*mnA)*(r50/(1-r50))
}
## ----


##' @rdname prRecruits
##' @export
## ---- prRecruitParsIB
prRecruitParsIB <- function(Msf,mnR,vrR,r=1,p=0,plus=0,Mbrak=seq(0,2.5,0.1),tol=1.0E-10) {
  ps <- prRecruitPars0(Msf,mnR,vrR,r,plus,Mbrak,tol)
  a <- ps$mnQ*(ps$mnQ*(1+ps$mnQ)+ps$vrQ)/ps$vrQ
  b <- ps$mnQ*(1+ps$mnQ)/ps$vrQ+2
  c(IB=list(c(a=a,b=b)),ps)
}
## ----


##' @rdname prRecruits
##' @importFrom stats rgamma
##' @export
## ---- prRecruitsG
prRecruitsG <- function(n,ps,mnA=1) {
  (ps$s0*mnA)*rgamma(n,ps$G[1L],ps$G[2L])
}
## ----


##' @rdname prRecruits
##' @importFrom stats qgamma
##' @export
## ---- prRecruitsQuantileG
prRecruitsQuantileG <- function(q,ps,mnA=1) {
  (ps$s0*mnA)*qgamma(q,ps$G[1L],ps$G[2L])
}
## ----


##' @rdname prRecruits
##' @export
## ---- prRecruitParsG
prRecruitParsG <- function(Msf,mnR,vrR,r=1,p=0,plus=0,Mbrak=seq(0,2.5,0.1),tol=1.0E-10) {
  ps <- prRecruitPars0(Msf,mnR,vrR,r,plus,Mbrak,tol)
  a <- ps$mnQ^2/ps$vrQ
  b <- ps$mnQ/ps$vrQ
  c(G=list(c(a=a,b=b)),ps)
}
## ----


##' @rdname prRecruits
##' @importFrom stats rlnorm
##' @export
## ---- prRecruitsLN
prRecruitsLN <- function(n,ps,mnA=1) {
  (ps$s0*mnA)*rlnorm(n,ps$LN[1L],ps$LN[2L])
}
## ----


##' @rdname prRecruits
##' @importFrom stats qlnorm
##' @export
## ---- prRecruitsQuantileLN
prRecruitsQuantileLN <- function(q,ps,mnA=1) {
  (ps$s0*mnA)*qlnorm(q,ps$LN[1L],ps$LN[2L])
}
## ----


##' @rdname prRecruits
##' @export
## ---- prRecruitParsLN
prRecruitParsLN <- function(Msf,mnR,vrR,r=1,p=0,plus=0,Mbrak=seq(0,2.5,0.1),tol=1.0E-10) {
  ps <- prRecruitPars0(Msf,mnR,vrR,r,plus,Mbrak,tol)
  m <- log(ps$mnQ^2/sqrt(ps$mnQ^2 + ps$vrQ))
  s <- sqrt(log(1 + ps$vrQ/ps$mnQ^2))
  c(LN=list(c(meanlog=m,sdlog=s)),ps)
}
## ----



##' @rdname prRecruits
##' @importFrom stats rbeta rbinom
##' @export
## ---- prRecruitsDIB
prRecruitsDIB <- function(n,ps,mnA=1) {
  R <- rbeta(n,ps$DIB[1L],ps$DIB[2L])
  (ps$s0*mnA)*rbinom(n,1,1-ps$DIB[3L])*(R/(1-R))
}
## ----


##' @rdname prRecruits
##' @importFrom stats qbeta
##' @export
## ---- prRecruitsQuantileDIB
prRecruitsQuantileDIB <- function(q,ps,mnA=1) {
  p <- ps$DIB[3L]
  rq <- qbeta((q-p)/(1-p),ps$DIB[1L],ps$DIB[2L])
  ifelse(q < p,0,(ps$s0*mnA)*(rq/(1-rq)))
}
## ----


##' @rdname prRecruits
##' @export
## ---- prRecruitParsDIB
prRecruitParsDIB <- function(Msf,mnR,vrR,r=1,p=0,plus=0,Mbrak=seq(0,2.5,0.1),tol=1.0E-10) {
  ps <- prRecruitPars0(Msf,mnR,vrR,r,plus,Mbrak,tol)
  mnQ <- ps$mnQ/(1-p)
  vrQ <- ps$vrQ/(1-p)-p*mnQ^2
  if(vrQ > 0) {
    a <- mnQ*(mnQ*(1+mnQ)+vrQ)/vrQ
    b <- mnQ*(1+mnQ)/vrQ+2
    c(DIB=list(c(a=a,b=b,p=p)),ps)
  }
}
## ----


##' @rdname prRecruits
##' @importFrom stats rgamma rbinom
##' @export
## ---- prRecruitsDG
prRecruitsDG <- function(n,ps,mnA=1) {
  (ps$s0*mnA)*rbinom(n,1,1-ps$DG[3L])*rgamma(n,ps$DG[1L],ps$DG[2L])
}
## ----


##' @rdname prRecruits
##' @importFrom stats qgamma
##' @export
## ---- prRecruitsQuantileDG
prRecruitsQuantileDG <- function(q,ps,mnA=1) {
  p <- ps$DG[3L]
  ifelse(q < p,0,(ps$s0*mnA)*qgamma((q-p)/(1-p),ps$DG[1L],ps$DG[2L]))
}
## ----


##' @rdname prRecruits
##' @export
## ---- prRecruitParsDG
prRecruitParsDG <- function(Msf,mnR,vrR,r=1,p=0,plus=0,Mbrak=seq(0,2.5,0.1),tol=1.0E-10) {
  ps <- prRecruitPars0(Msf,mnR,vrR,r,plus,Mbrak,tol)
  mnQ <- ps$mnQ/(1-p)
  vrQ <- ps$vrQ/(1-p)-p*mnQ^2
  if(vrQ > 0) {
    a <- mnQ^2/vrQ
    b <- mnQ/vrQ
    c(DG=list(c(a=a,b=b,p=p)),ps)
  }
}
## ----




##' @rdname prRecruits
##' @importFrom stats rlnorm rbinom
##' @export
## ---- prRecruitsDLN
prRecruitsDLN <- function(n,ps,mnA=1) {
  (ps$s0*mnA)*rbinom(n,1,1-ps$DLN[3L])*rlnorm(n,ps$DLN[1L],ps$DLN[2L])
}
## ----


##' @rdname prRecruits
##' @importFrom stats qlnorm
##' @export
## ---- prRecruitsQuantileDLN
prRecruitsQuantileDLN <- function(q,ps,mnA=1) {
  p <- ps$DLN[3L]
  ifelse(q < p,0,(ps$s0*mnA)*qlnorm((q-p)/(1-p),ps$DLN[1L],ps$DLN[2L]))
}
## ----


##' @rdname prRecruits
##' @export
## ---- prRecruitParsDLN
prRecruitParsDLN <- function(Msf,mnR,vrR,r=1,p=0,plus=0,Mbrak=seq(0,2.5,0.1),tol=1.0E-10) {
  ps <- prRecruitPars0(Msf,mnR,vrR,r,plus,Mbrak,tol)
  mnQ <- ps$mnQ/(1-p)
  vrQ <- ps$vrQ/(1-p)-p*mnQ^2
  if(vrQ > 0) {
    m <- log(mnQ^2/sqrt(mnQ^2 + vrQ))
    s <- sqrt(log(1 + vrQ/mnQ^2))
    c(DLN=list(c(meanlog=m,sdlog=s,p=p)),ps)
  }
}
## ----



##' Resample mean and variance of R
##'
##' This is the resampling procedure for the mean and variance of R
##' implemented in the GYM. This procedure is extremely inadvisable
##' and should never be used.
##'
##' @title Resample R
##' @param mnR the mean of the proportion R
##' @param vrR the variance of the proportion R
##' @param n number of independent surveys used in the estimation of
##'   the mean and variance of R
##' @return A list with elements
##' * `mnR` the resampled mean R
##' * `vrR` the resampled variance of R
##' @example inst/examples/resampleRGYM.R
##' @importFrom stats rchisq rnorm
##' @export
## ---- resampleRGYM
resampleRGYM <- function(mnR,vrR,n) {
  vrR <- rchisq(1,n-1L)/(n-1L)*vrR
  mnR <- rnorm(1,mnR,sqrt(vrR)/n)
  list(mnR=pmax(0,pmin(1,mnR)),vrR=vrR)
}
## ----


##' Simulate the mean and variance for the proportional recruits model
##' by parametric bootstrap when the distribution of the fraction of
##' recruits is inverse Beta (IB), Gamma (G) or log Normal (LN).
##'
##' The proportional recruit model defines R, the number of
##' individuals in age class r as a proportion of the population of
##' the same age or older, and assumes estimates of the mean and
##' variance of R are available from independent surveys. To allow for
##' the variability in these estimates, these functions use
##' parameteric bootstrap to simulate new values for the mean and
##' variance of R that are consistent with the observed values and the
##' chosen distributional model.
##'
##' @title Proportional Recruit Bootstrap
##' @param prRec a function for generating random recruits by a
##'   proportional recruitment model
##' @param ps a list of model parameters
##' @param n number of independent surveys used in the estimation of
##'   the mean and variance of R
##' @param Msf matrix of final unscaled integrated natural mortalities
##' @param mnA mean number of recruits
##' @param r the reference age class
##' @param plus the number of additional age classes in the plus class
##'   (`plus=0` implies no plus class).
##' @return Return a list with elements
##' * `mnR` simulated mean of R
##' * `vrR` simulated variance of R
##' * `p` simulated probability of zero recruitment
##' * `n` number os surveys used in simulation
##' * `r` reference age class
##' @example inst/examples/prBootstrap.R
##' @importFrom stats var
##' @export
## ---- prBootstrap
prBootstrap <- function(prRec,ps,n,Msf,mnA=1,r=1,plus=0) {
  ## Generate n boostrap samples
  Rs <- sapply(seq_len(n),function(k) {
    ## Generate a random age structure and compute R
    N <- ageStructureS(prRec(length(Msf)+plus,ps,mnA),Msf,ps$M,plus=plus)
    N[r]/sum(N[r:length(N)])})

  list(mnR=mean(Rs),vrR=var(Rs),p=mean(Rs==0),n=n,r=r)
}
## ----





##' Compute the scalings required to adjust surveyed age class
##' abundances to initial abundances at a reference age.
##'
##' Given the age class, the year and the time steps within the year
##' at which the age class was surveyed, this function computes the
##' total survival from the start of the year in which the cohort were
##' in the reference age class to the survey period.  If the surveyed
##' age class is younger than the reference class, the reciporical of
##' the total survival from the survey period to the start of the year
##' that the cohort will be in the reference class is computed.
##'
##' If there is inter-annual variability in natural or fishing
##' mortality, the survey years must be labelled so that `yr==1`
##' corresponds to the first element of the vector of `M` and/or `F`,
##' and it is no possible to compute the survival for cohorts that
##' recruit before year 1.
##'
##' If there is no inter-annual variablity in natural or fishing
##' mortality, the survey year is irrelevant.
##'
##' @title Survival to a survey period.
##' @param yr vector of survey year
##' @param cls vector of survey age class
##' @param s1 vector of the first time step in the survey
##' @param s2 vector of the final time step in the survey
##' @param Ms matrix of *unscaled* integrated natural mortality
##' @param M vector of annual natural mortalities
##' @param Fs matrix of *unscaled* integrated fishing mortality
##' @param F vector of annual fishing mortalities
##' @param rcls the reference age class to adjust to
##' @return Returns a vector of the mean survival from time of
##'   recruitment to the survey period.
##' @example inst/examples/surveySurvival.R
##' @export
## ---- surveySurvival
surveySurvival <- function(yr,cls,s1,s2,Ms,M,Fs=0,F=0,rcls=1) {


  if(length(M)==1L && length(F)==1L) {
    ## Precompute survivals
    N <- exp(-M*Ms-F*Fs)
    ## Survival from first class to start of survey year
    S <- cumprod(c(1,N[nrow(N),]))[cls]
    ## Adjust for survival within survey year
    for(k in seq_along(S)) {
      S[k] <- S[k]*mean(N[s1[k]:s2[k],cls[k]])
    }
    ## Adjust if target class is not 1
    if(rcls > 1) S <- S/cumprod(c(1,N[nrow(N),]))[rcls]

  } else {

    ryr <- yr-(cls-rcls)
    if(length(M)==1) M <- rep.int(M,max(yr,ryr))
    if(length(F)==1) F <- rep.int(F,max(yr,ryr))

    ## Cannot project outside the range of Ms, Fs
    S <- ifelse(pmin(yr,ryr) >= 1 & pmax(yr,ryr) <= min(length(F),length(M)),1,NA)

    ## Iterate over years in which surveyed individuals were alive
    for(y in seq.int(max(1,min(yr,ryr)),min(length(M),max(yr,ryr)))) {
      ## Compute survival
      N <- exp(-M[y]*Ms-F[y]*Fs)
      ## Age class of each individual in year y
      a <- cls-yr+y
      ## Back project survival to target class
      for(k in which(a>=rcls & a<cls))
        S[k] <- S[k]*N[nrow(N),a[k]]
      ## Forward project survival to target class
      for(k in which(a<rcls & a>=cls))
        S[k] <- S[k]/N[nrow(N),a[k]]
      ## Survival in survey year
      for(k in which(y==yr))
        S[k] <- S[k]*mean(N[s1[k]:s2[k],a[k]])
    }
  }
  S
}
## ----



##' Adjust surveyed age class abundances to initial abundances at a
##' reference age class as implemented in the GYM.
##'
##' This function repreoduces the procedure used in the GYM to adjust
##' surveyed abundances back to recruitment to a reference age class.
##' In particular, where abundance is estimated from multiple surveys,
##' abundance is estimated as a geometric mean weighted by the inverse
##' of the squared *coefficient of variation* of the surveyed density.
##'
##' The survey data must be provided as a dataframe with columns
##'
##' * `Year` - year of survey
##' * `AgeClass` - age class surveyed
##' * `Density` - estimated density
##' * `SE` - standard error of density
##' * `Scale` - factor to scale density to population abundance
##' * `Start` - increment at which survey starts
##' * `End` - increment at which survey ends
##'
##' Estimates are generated for a contiguous sequence of years from
##' the earliest year to the latest year for which the abundance of
##' the reference class can be estimated, with `NA` returned for
##' years for which no estimate can be made.
##'
##' @title Estimate recruitment from surveyed abundance
##' @param survey.df a dataframe of surved abundance (see details)
##' @param Ms matrix of *unscaled* integrated natural mortality
##' @param M vector of annual natural mortalities
##' @param rcls the reference age class to adjust to
##' @return A dataframe with columns
##' * `Year` the recruitment year
##' * `RAbund` the abundance of the reference class
##' @example inst/examples/surveyAdjustGYM.R
##' @export
## ---- surveyAdjustGYM
surveyAdjustGYM <- function(survey.df,Ms,M,rcls) {
  ## Scale density to abundance
  ab.mn <- survey.df$Scale*survey.df$Density
  ab.se <- survey.df$Scale*survey.df$SE

  ## Compute log survival adjustment
  S <- surveySurvival(survey.df$Year,survey.df$AgeClass,survey.df$Start,survey.df$End,Ms,M,rcls=rcls)

  ## Weight rescaled abundance by 1/cv^2
  ab.wt <- (ab.mn/ab.se)^2
  ## Compute the year of "recruitment" to the target age class
  rc.yr <- survey.df$Year-survey.df$AgeClass+rcls
  yr <- seq.int(min(rc.yr),max(rc.yr))
  rec.yf <- factor(rc.yr,yr)
  ## Compute the weighted geometric means
  data.frame(
    Year=yr,
    RAbund=exp(tapply(ab.wt*(log(ab.mn)-log(S)),rec.yf,sum)/tapply(ab.wt,rec.yf,sum)))
}
## ----
