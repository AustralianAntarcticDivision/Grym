
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Grym

An implementation of the generalized yield model (Constable and de la
Mare, 1996) in R. This implementation is still considered a work in
progress up until the paper of comparisons has been submitted to
WG-SAM-2020.

Once accepted by CCAMLR the accepted version of the Grym Package will be
hosted as a fork on the CCAMLR github account.

## Installing

The package is easily installed from GitHub, using the remotes package.

``` r
remotes::install_github("AustralianAntarcticDivision/Grym")
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

Grym is accompanied by the GrymExamples package which provides examples
of implementation for a number of species within the CCAMLR region.

## References

Constable, AJ, and WK de la Mare. 1996. “A Generalised Model for
Evaluating Yield and the Long-Term Status of Fish Stocks Under
Conditions of Uncertainty.” CCAMLR Science 3: 31–54.
