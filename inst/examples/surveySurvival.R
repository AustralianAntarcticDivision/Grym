## Daily time steps and 7 age classes
nsteps <- 365
Ages <- 2:8
Days <- seq(0,1,length=nsteps+1)
h <- 1/nsteps

## Constant intra-annual natural mortality
ms <- matrix(1,nsteps+1,length(Ages))
ms <- ms/mean(trapz(ms,h))
Ms <- ctrapz(ms,h)

## Survey year, period and age classes
svy <- data.frame(yr=3:5,s1=c(190,220,150),s2=c(201,231,161))
svy <- cbind(svy[rep(1:3,each=7),],cls=1:7)
head(svy)

## Constant mortality
M <- 0.2

## Survival to the survey period from age class 1
surveySurvival(svy$yr,svy$cls,svy$s1,svy$s2,Ms,M)

## Survival to the survey period from age class 3
surveySurvival(svy$yr,svy$cls,svy$s1,svy$s2,Ms,M,rcls=3)

## Variable mortality
M <- rgamma(10,20,100)
M

## Survival cannot be projected outside the period for which mortality
## is specified.
surveySurvival(svy$yr,svy$cls,svy$s1,svy$s2,Ms,M)
