% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mse_saeOBns.R
\name{mse_saeOBns}
\alias{mse_saeOBns}
\title{Parametric Bootstrap Mean Squared Error Estimators of Optimum Benchmarking for Univariate Non Sampled Area in Small Area Estimation}
\usage{
mse_saeOBns(
  formula,
  vardir,
  weight,
  cluster,
  samevar = FALSE,
  B = 100,
  MAXITER = 100,
  PRECISION = 1e-04,
  data
)
}
\arguments{
\item{formula}{an object of class list of formula describe the fitted model}

\item{vardir}{vector containing sampling variances of direct estimators}

\item{weight}{vector containing proportion of units in small areas}

\item{cluster}{vector containing cluster of auxiliary variable}

\item{samevar}{logical. If \code{TRUE}, the varians is same. Default is \code{FALSE}}

\item{B}{number of bootstrap. Default is 1000}

\item{MAXITER}{maximum number of iterations for Fisher-scoring. Default is 100}

\item{PRECISION}{coverage tolerance limit for the Fisher Scoring algorithm. Default value is \code{1e-4}}

\item{data}{dataframe containing the variables named in formula, vardir, and weight}
}
\value{
\item{mse.eblup}{estimated mean squared errors of the EBLUPs for the small domains based on Prasad Rao}
\item{pbmse.eblupOB}{parametric bootstrap mean squared error estimates of the optimum benchmark}
\item{running.time}{time for running function}
}
\description{
Calculates the parametric bootstrap mean squared error estimates of optimum benchmarking for univariate non sampled area in small area estimation
}
\examples{
\donttest{
## load dataset
data(datamsaeOBns)

# Compute MSE EBLUP and Optimum Benchmark

## Using parameter 'data'
mse_sae = mse_saeOBns(Y1 ~ X1 + X2, v1, w1, c1, data = datamsaeOBns)

## Without parameter 'data'
mse_sae = mse_saeOBns(datamsaeOBns$Y1 ~ datamsaeOBns$X1 + datamsaeOBns$X2,
datamsaeOBns$v1, datamsaeOBns$w1, datamsaeOBns$c1)

## Return
mse_sae$pbmse.eblupOB # to see the MSE Optimum Benchmark estimators
}
}
