
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Grym

An implementation of the generalized yield model (Constable and de la
Mare, 1996) in R. This implementation is still considered a work in
progress up until it has been endorced by the Scientific Committee.

Once accepted by CCAMLR the accepted version of the Grym Packahosted as
a fork on the CCAMLR github account.  
ge will be \#\# Installing

The package is easily installed from GitHub, using the remotes package.

``` r
remotes::install_github("AustralianAntarcticDivision/Grym", build_vignettes=TRUE)
```

If you don’t have `remotes` installed already, install it first.

``` r
install.packages("remotes")
```

Grym does not otherwise require `remotes` for normal use.

## TODO

  - Parametric bootstrap from a vector of recruits

  - Check examples for consistency of implementation

  - Build GrymExamples

# Grym

<!-- badges: start -->

<!-- badges: end -->

The goal of Grym is to have provide an R based implementation of the
Generalised Yield Model by Constable and de la Mare (1996). This package
provides specific functions which users can use in combination to build
their own projection functions and packages.

## Example

Grym is accompanied by the GrymExamples package which provides vignettes
of examples of implementation for a number of species within the CCAMLR
region. The GrymExamples package can be found here:
github.com/AustralianAntarcticDivision/GrymExamples

Unless otherwise stated all examples are works in progress. A list of
examples can be accessed with the following code:

``` r
vignette(package="GrymExamples")
```

Specific examples can be accessed with the same function and specifying
the topic for example

``` r
vignette(topic = "Icefish_2019", package="GrymExamples")
```

## References

Constable, AJ, and WK de la Mare. 1996. “A Generalised Model for
Evaluating Yield and the Long-Term Status of Fish Stocks Under
Conditions of Uncertainty.” CCAMLR Science 3: 31–54.
