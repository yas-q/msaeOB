
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msaeOB

<!-- badges: start -->

<!-- badges: end -->

Implements multivariate optimum benchmarking small area estimation. This
package provides optimum benchmarking estimation for univariate and
multivariate small area estimation and its MSE. In fact, MSE estimators
for optimum benchmark are not readily available, so resampling method
that called parametric bootstrap is applied. The optimum benchmark model
and parametric bootstrap in this package are based on the model proposed
in small area estimation (J.N.K Rao and Isabel Molina, 2015).

## Authors

Muhammad Yasqi Imanda, Zenda Oka Briantiko, Azka Ubaidillah

## Maintainer

Muhammad Yasqi Imanda <221810403@stis.ac.id>

## Installation

You can install the released version of msaeOB from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("msaeOB")
```

## Functions

  - est\_saeOB : Produces EBLUPs Optimum Benchmarking based on a
    Univariate Fay-Herriot (Model 1)
  - mse\_saeOB : Parametric Bootstrap Mean Squared Error Estimators of
    Optimum Benchmarking for Univariate Small Area Estimation
  - est\_msaeOB : Produces EBLUPs Optimum Benchmarking based on a
    Multivariate Fay-Herriot (Model 1)
  - mse\_msaeOB : Parametric Bootstrap Mean Squared Error Estimators of
    Optimum Benchmarking for Multivariate Small Area Estimation
  - est\_saeOBns : Produces EBLUPs Optimum Benchmarking for Non Sampled
    Area based on a Univariate Fay-Herriot (Model 1)
  - mse\_saeOBns : Parametric Bootstrap Mean Squared Error Estimators of
    Optimum Benchmarking for Univariate Non Sampled Area in Small Area
    Estimation
  - est\_msaeOBns : Produces EBLUPs Optimum Benchmarking for Non Sampled
    Area based on a Multivariate Fay Herriot (Model 1)
  - mse\_msaeOBns : Parametric Bootstrap Mean Squared Error Estimators
    of Optimum Benchmarking for Multivariate Non Sampled Area in Small
    Area Estimation

## References

  - Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New
    York: John Wiley and Sons, Inc.
  - Benavent, Roberto & Morales, Domingo. (2015). ???Multivariate
    Fay-Herriot models for small area estimation???. Computational
    Statistics and Data Analysis 94 2016 372-390. DOI:
    10.1016/j.csda.2015.07.013.
  - Ubaidillah, Azka et al.??(2019). Multivariate Fay-Herriot models for
    small area estimation with application to household consumption per
    capita expenditure in Indonesia. Journal of Applied Statistics.
    46:15. 2845-2861. DOI: 10.1080/02664763.2019.1615420.
  - Wang, J., Fuller, W.A., and Qu, Y. (2008). Small Area Estimation
    Under Restriction. Survey Methodology. 34. 29???36.
  - You, Y., Rao, J.N.K., and Hidiroglou, M.A.??(2013). On the
    Performance of Self-Benchmarked Small Area Estimators Under the
    Fay-Herriot Area Level Model. Survey Methodology, 39, 217???229.
  - Krzciuk, M. K. (2018). On the Simulation Study of Jackknife and
    Bootstrap MSE Estimators of a Domain Mean Predictor for Fay???Herriot
    Model. Acta Universitatis Lodziensis. Folia Oeconomica, 5(331),
    169-183.
