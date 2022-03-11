#' @title EBLUPs Optimum Benchmarking based on a Multivariate Fay Herriot (Model 1)
#'
#' @description This function gives EBLUPs optimum benchmarking based on multivariate Fay-Herriot (Model 1)
#'
#' @param formula an object of class list of formula describe the fitted models
#' @param vardir matrix containing sampling variances of direct estimators. The order is: \code{var1, cov12, ..., cov1r, var2, cov23, ..., cov2r, ..., cov(r-1)(r), var(r)}
#' @param weight matrix containing proportion of units in small areas. The order is: \code{w1, w2, ..., w(r)}
#' @param samevar logical. If \code{TRUE}, the varians is same. Default is \code{FALSE}
#' @param MAXITER maximum number of iterations for Fisher-scoring. Default is 100
#' @param PRECISION coverage tolerance limit for the Fisher Scoring algorithm. Default value is \code{1e-4}
#' @param data dataframe containing the variables named in formula, vardir, and weight
#'
#' @return This function returns a list with following objects:
#' \item{eblup}{a list containing a value of estimators}
#' \itemize{
#'   \item est.eblup : a dataframe containing EBLUP estimators
#'   \item est.eblupOB : a dataframe containing optimum benchmark estimators
#' }
#'
#' \item{fit}{a list containing following objects:}
#' \itemize{
#'   \item method : fitting method, named "REML"
#'   \item convergence : logical value of convergence of Fisher Scoring
#'   \item iterations : number of iterations of Fisher Scoring algorithm
#'   \item estcoef : a data frame containing estimated model coefficients (\code{beta, std. error, t value, p-value})
#'   \item refvar : estimated random effect variance
#' }
#' \item{random.effect}{a data frame containing values of random effect estimators}
#' \item{agregation}{a data frame containing agregation of direct, EBLUP, and optimum benchmark estimation}
#'
#' @export est_msaeOB
#'
#' @import abind
#' @importFrom magic adiag
#' @importFrom Matrix forceSymmetric
#' @importFrom stats model.frame na.omit model.matrix median pnorm rnorm
#' @importFrom MASS mvrnorm
#'
#' @examples
#' ## load dataset
#' data(datamsaeOB)
#'
#' # Compute EBLUP and Optimum Benchmark using auxiliary variables X1 and X2 for each dependent variable
#'
#' ## Using parameter 'data'
#' Fo = list(f1 = Y1 ~ X1 + X2,
#'           f2 = Y2 ~ X1 + X2,
#'           f3 = Y3 ~ X1 + X2)
#' vardir = c("v1", "v12", "v13", "v2", "v23", "v3")
#' weight = c("w1", "w2", "w3")
#'
#' est_msae = est_msaeOB(Fo, vardir, weight, data = datamsaeOB)
#'
#' ## Without parameter 'data'
#' Fo = list(f1 = datamsaeOB$Y1 ~ datamsaeOB$X1 + datamsaeOB$X2,
#'           f2 = datamsaeOB$Y2 ~ datamsaeOB$X1 + datamsaeOB$X2,
#'           f3 = datamsaeOB$Y3 ~ datamsaeOB$X1 + datamsaeOB$X2)
#' vardir = datamsaeOB[, c("v1", "v12", "v13", "v2", "v23", "v3")]
#' weight = datamsaeOB[, c("w1", "w2", "w3")]
#'
#' est_msae = est_msaeOB(Fo, vardir, weight)
#'
#' ## Return
#' est_msae$eblup$est.eblupOB # to see the Optimum Benchmark estimators
#'
est_msaeOB<-function (formula, vardir, weight, samevar = FALSE, MAXITER = 100,
                      PRECISION = 1e-04, data)
{
  r = length(formula)
  if (r <= 1)
    stop("You should use est_saeOB() for univariate")
  R_function = function(vardir, n, r) {
    if (r == 1) {
      R = diag(vardir)
    }
    else {
      R = matrix(rep(0, times = n * r * n * r), nrow = n *
                   r, ncol = n * r)
      k = 1
      for (i in 1:r) {
        for (j in 1:r) {
          if (i <= j) {
            mat0 = matrix(rep(0, times = r * r), nrow = r,
                          ncol = r)
            mat0[i, j] = 1
            matVr = diag(vardir[, k], length(vardir[,
                                                    k]))
            R_hasil = kronecker(mat0, matVr)
            R = R + R_hasil
            k = k + 1
          }
        }
      }
      R = forceSymmetric(R)
      R = R
    }
    return(as.matrix(R))
  }
  if (!missing(data)) {
    formuladata = lapply(formula, function(x) model.frame(x,
                                                          na.action = na.omit, data))
    y = unlist(lapply(formula, function(x) model.frame(x,
                                                       na.action = na.omit, data)[[1]]))
    X = Reduce(adiag, lapply(formula, function(x) model.matrix(x,
                                                               data)))
    W = as.matrix(data[, weight])
    n = length(y)/r
    if (any(is.na(data[, vardir])))
      stop("Object vardir contains NA values.")
    if (!all(vardir %in% names(data)))
      stop("Object vardir is not appropriate with data.")
    if (length(vardir) != sum(1:r))
      stop("Length of vardir is not appropriate with data. The length must be ",
           sum(1:r))
    if (any(is.na(data[, weight])))
      stop("Object weight contains NA values.")
    if (!all(weight %in% names(data)))
      stop("Object weight is not appropriate with data.")
    if (length(weight) != r)
      stop("Length of weight is not appropriate with data. The length must be ",
           r)
    R = R_function(data[, vardir], n, r)
    vardir = data[, vardir]
  }
  else {
    formuladata = lapply(formula, function(x) model.frame(x,
                                                          na.action = na.omit))
    y = unlist(lapply(formula, function(x) model.frame(x,
                                                       na.action = na.omit)[[1]]))
    X = Reduce(adiag, lapply(formula, function(x) model.matrix(x)))
    W = as.matrix(weight)
    n = length(y)/r
    if (any(is.na(vardir)))
      stop("Object vardir contains NA values")
    if ((dim(vardir)[1] != n) || (dim(vardir)[2] != sum(1:r)))
      stop("Object vardir is not appropriate with data. It must be ",
           n, " x ", sum(1:r), " matrix.")
    if (any(is.na(weight)))
      stop("Object weight contains NA values.")
    if ((dim(weight)[1] != n) || (dim(weight)[2] != r))
      stop("Object weight is not appropriate with data. It must be ",
           n, " x ", r, " matrix.")
    R = R_function(vardir, n, r)
  }
  y_names = sapply(formula, "[[", 2)
  Ir = diag(r)
  In = diag(n)
  dV = list()
  dV1 = list()
  for (i in 1:r) {
    dV[[i]] = matrix(0, nrow = r, ncol = r)
    dV[[i]][i, i] = 1
    dV1[[i]] = kronecker(dV[[i]], In)
  }
  convergence = TRUE
  if (samevar) {
    Vu = median(diag(R))
    k = 0
    diff = rep(PRECISION + 1, r)
    while (any(diff > PRECISION) & (k < MAXITER)) {
      k = k + 1
      Vu1 = Vu
      Gr = kronecker(Vu1, Ir)
      Gn = kronecker(Gr, In)
      V = as.matrix(Gn + R)
      Vinv = solve(V)
      XtVinv = t(Vinv %*% X)
      Q = solve(XtVinv %*% X)
      P = Vinv - t(XtVinv) %*% Q %*% XtVinv
      Py = P %*% y
      s = (-0.5) %*% sum(diag(P)) + 0.5 %*% (t(Py) %*%
                                               Py)
      iF = 0.5 %*% sum(diag(P %*% P))
      Vu = Vu1 + solve(iF) %*% s
      diff = abs((Vu - Vu1)/Vu1)
    }
    Vu = as.vector((rep(max(Vu, 0), r)))
    names(Vu) = y_names
    if (k >= MAXITER && diff >= PRECISION) {
      convergence = FALSE
    }
    Gn = kronecker(diag(Vu), In)
    V = as.matrix(Gn + R)
    Vinv = solve(V)
    XtVinv = t(Vinv %*% X)
    Q = solve(XtVinv %*% X)
    P = Vinv - t(XtVinv) %*% Q %*% XtVinv
    Py = P %*% y
    beta = Q %*% XtVinv %*% y
    res = y - X %*% beta
    eblup = data.frame(matrix(X %*% beta + Gn %*% Vinv %*%
                                res, n, r))
    names(eblup) = y_names
    se.b = sqrt(diag(Q))
    t.value = beta/se.b
    p.value = 2 * pnorm(abs(as.numeric(t.value)), lower.tail = FALSE)
    coef = as.matrix(cbind(beta, se.b, t.value, p.value))
    colnames(coef) = c("beta", "std. error",
                       "t value", "p-value")
    rownames(coef) = colnames(X)
    Bi <- R %*% solve(V)
    m <- dim(X)[1]
    I <- diag(m)
    g1d <- diag(Bi %*% Gn)
    g2d <- diag(Bi %*% X %*% Q %*% t(X) %*%
                  t(Bi))
    dg <- Vinv - (I - Bi) %*% Vinv
    g3d <- diag(dg %*% V %*% t(dg))/iF
    mse <- g1d + g2d + 2 * g3d
    mse <- data.frame(matrix(mse, n, r))
    names(mse) = y_names
  }
  else {
    Vu = apply(matrix(diag(R), nrow = n, ncol = r), 2, median)
    k = 0
    diff = rep(PRECISION + 1, r)
    while (any(diff > rep(PRECISION, r)) & (k < MAXITER)) {
      k = k + 1
      Vu1 = Vu
      if (r == 1) {
        Gr = Vu1
      }
      else {
        Gr = diag(as.vector(Vu1))
      }
      Gn = kronecker(Gr, In)
      V = as.matrix(Gn + R)
      Vinv = solve(V)
      XtVinv = t(Vinv %*% X)
      Q = solve(XtVinv %*% X)
      P = Vinv - t(XtVinv) %*% Q %*% XtVinv
      Py = P %*% y
      s = sapply(dV1, function(x) (-0.5) * sum(diag(P %*%
                                                      x)) + 0.5 * (t(Py) %*% x %*% Py))
      iF = matrix(unlist(lapply(dV1, function(x) lapply(dV1,
                                                        function(y) 0.5 * sum(diag(P %*% x %*% P %*%
                                                                                     y))))), r)
      Vu = Vu1 + solve(iF) %*% s
      diff = abs((Vu - Vu1)/Vu1)
    }
    Vu = as.vector(sapply(Vu, max, 0))
    if (k >= MAXITER && diff >= PRECISION) {
      convergence = FALSE
    }
    if (r == 1) {
      Gr = Vu1
    }
    else {
      Gr = diag(as.vector(Vu1))
    }
    Gn = kronecker(Gr, In)
    V = as.matrix(Gn + R)
    Vinv = solve(V)
    XtVinv = t(Vinv %*% X)
    Q = solve(XtVinv %*% X)
    P = Vinv - t(XtVinv) %*% Q %*% XtVinv
    Py = P %*% y
    beta = Q %*% XtVinv %*% y
    res = y - X %*% beta
    eblup = data.frame(matrix(X %*% beta + Gn %*% Vinv %*%
                                res, n, r))
    names(eblup) = y_names
    se.b = sqrt(diag(Q))
    t.value = beta/se.b
    p.value = 2 * pnorm(abs(as.numeric(t.value)), lower.tail = FALSE)
    coef = as.matrix(cbind(beta, se.b, t.value, p.value))
    colnames(coef) = c("beta", "std. error",
                       "t value", "p-value")
    rownames(coef) = colnames(X)
    iF_inv <- solve(iF)
    Bi <- R %*% solve(V)
    m <- dim(X)[1]
    I <- diag(m)
    g1d <- diag(Bi %*% Gn)
    g2d <- diag(Bi %*% X %*% Q %*% t(X) %*%
                  t(Bi))
    dg <- lapply(dV1, function(x) x %*% Vinv - Gn %*%
                   Vinv %*% x %*% Vinv)
    g3d = list()
    for (i in 1:r) {
      for (j in 1:r) {
        g3d[[(i - 1) * r + j]] = iF_inv[i, j] * (dg[[i]] %*%
                                                   V %*% t(dg[[j]]))
      }
    }
    g3d <- diag(Reduce("+", g3d))
    mse <- g1d + g2d + 2 * g3d
    mse <- data.frame(matrix(mse, n, r))
    names(mse) = y_names
  }
  random.effect = data.frame(matrix(Gn %*% Vinv %*% res, n,
                                    r))
  names(random.effect) = y_names
  y.mat = matrix(y, nrow = n, ncol = r)
  eblup.mat = as.matrix(eblup)
  alfa = matrix(0, nrow = n, ncol = r)
  lambda = matrix(0, nrow = n, ncol = r)
  eblup.optimum = matrix(0, nrow = n, ncol = r)
  for (i in 1:r) {
    alfa[, i] = (sum(W[, i] * y.mat[, i]) - sum(W[, i] * eblup.mat[, i]))
  }
  for (i in 1:r) {
    for (j in 1:n) {
      lambda[j , i] = (mse[j , i] * W[j , i]) / sum(mse[, i] * W[, i]^2)
    }
  }
  for (i in 1:r) {
    for (j in 1:n) {
      eblup.optimum[j , i] = eblup.mat[j , i] + (lambda[j , i] * alfa[j , i])
    }
  }
  eblup.optimum = as.data.frame(eblup.optimum)
  names(eblup.optimum) = y_names
  agregation.direct = diag(t(W) %*% y.mat)
  agregation.eblup = diag(t(W) %*% eblup.mat)
  agregation.eblup.optimum = diag(t(W) %*% as.matrix(eblup.optimum))
  agregation = as.matrix(rbind(agregation.direct, agregation.eblup,
                               agregation.eblup.optimum))
  colnames(agregation) = y_names
  result = list(eblup = list(est.eblup = NA, est.eblupOB = NA),
                fit = list(method = NA, convergence = NA, iteration = NA,
                           estcoef = NA, refvar = NA), random.effect = NA, agregation = NA)
  result$eblup$est.eblup = eblup
  result$eblup$est.eblupOB = eblup.optimum
  result$fit$method = "REML"
  result$fit$convergence = convergence
  result$fit$iteration = k
  result$fit$estcoef = coef
  result$fit$refvar = t(Vu)
  result$random.effect = random.effect
  result$agregation = agregation
  return(result)
}
