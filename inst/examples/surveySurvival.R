## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(from=0, to=1, length=nsteps+1)
h <- 1/nsteps

## Constant intra-annual natural mortality
ms <- matrix(data=1, nrow=nsteps+1, ncol=length(x=Ages))
ms <- ms/mean(x=trapz(fs=ms, h=h))
Ms <- ctrapz(fs=ms, h=h)

## Survey year,  period and age classes
svy <- data.frame(yr=3:5, s1=c(190, 220, 150), s2=c(201, 231, 161))
svy <- cbind(svy[rep(x=1:3, each=7), ], cls=1:7)
head(svy)

## Constant mortality
M <- 0.2

## Survival to the survey period from age class 1
surveySurvival(yr=svy$yr, cls=svy$cls, s1=svy$s1, s2=svy$s2, Ms=Ms, M=M)

## Survival to the survey period from age class 3
surveySurvival(yr=svy$yr, cls=svy$cls, s1=svy$s1, s2=svy$s2, Ms=Ms, M=M, rcls=3)

## Variable mortality
M <- rgamma(n=10, shape=20, rate=100)
M

## Survival cannot be projected outside the period for which mortality
## is specified.
surveySurvival(yr=svy$yr, cls=svy$cls, s1=svy$s1, s2=svy$s2, Ms=Ms, M=M)
