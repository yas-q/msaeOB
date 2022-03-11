#' @title Parametric Bootstrap Mean Squared Error Estimators of Optimum Benchmarking for Multivariate Small Area Estimation
#'
#' @description Calculates the parametric bootstrap mean squared error estimates of optimum benchmarking for multivariate small area estimation
#'
#' @param formula an object of class list of formula describe the fitted models
#' @param vardir matrix containing sampling variances of direct estimators. The order is: \code{var1, cov12, ..., cov1r, var2, cov23, ..., cov2r, ..., cov(r-1)(r), var(r)}
#' @param weight matrix containing proportion of units in small areas. The order is: \code{w1, w2, ..., w(r)}
#' @param samevar logical. If \code{TRUE}, the varians is same. Default is \code{FALSE}
#' @param B number of bootstrap. Default is 1000
#' @param MAXITER maximum number of iterations for Fisher-scoring. Default is 100
#' @param PRECISION coverage tolerance limit for the Fisher Scoring algorithm. Default value is \code{1e-4}
#' @param data dataframe containing the variables named in formula, vardir, and weight
#'
#' @return
#' \item{mse.eblup}{estimated mean squared errors of the EBLUPs for the small domains based on Prasad Rao}
#' \item{pbmse.eblupOB}{parametric bootstrap mean squared error estimates of the optimum benchmark}
#' \item{running.time}{time for running function}
#'
#' @export mse_msaeOB
#'
#' @import abind
#' @importFrom magic adiag
#' @importFrom Matrix forceSymmetric
#' @importFrom stats model.frame na.omit model.matrix median pnorm rnorm cor
#' @importFrom MASS mvrnorm ginv
#'
#' @examples
#' \donttest{
#' ## load dataset
#' data(datamsaeOB)
#'
#' # Compute MSE EBLUP and Optimum Benchmark
#' # This is the long running example
#' ## Using parameter 'data'
#' Fo = list(f1 = Y1 ~ X1 + X2,
#'           f2 = Y2 ~ X1 + X2,
#'           f3 = Y3 ~ X1 + X2)
#' vardir = c("v1", "v12", "v13", "v2", "v23", "v3")
#' weight = c("w1", "w2", "w3")
#'
#' mse_msae = mse_msaeOB(Fo, vardir, weight, data = datamsaeOB)
#'
#' ## Without parameter 'data'
#' Fo = list(f1 = datamsaeOB$Y1 ~ datamsaeOB$X1 + datamsaeOB$X2,
#'           f2 = datamsaeOB$Y2 ~ datamsaeOB$X1 + datamsaeOB$X2,
#'           f3 = datamsaeOB$Y3 ~ datamsaeOB$X1 + datamsaeOB$X2)
#' vardir = datamsaeOB[, c("v1", "v12", "v13", "v2", "v23", "v3")]
#' weight = datamsaeOB[, c("w1", "w2", "w3")]
#'
#' mse_msae = mse_msaeOB(Fo, vardir, weight)
#'
#' ## Return
#' mse_msae$pbmse.eblupOB # to see the MSE of Optimum Benchmark
#' }
mse_msaeOB<-function (formula, vardir, weight, samevar = FALSE, B = 100,
                      MAXITER = 100, PRECISION = 1e-04, data)
{
  start_time <- Sys.time()
  r = length(formula)
  if (r <= 1)
    stop("You should use mse_saeOB() for univariate")
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
    }
    return(as.matrix(R))
  }
  temporary_vardir = function(vardir, n, r) {
    dt = list()
    for (h in 1:n) {
      var_mat = matrix(0, nrow = r, ncol = r)
      k = 1
      for (i in 1:r) {
        for (j in 1:r) {
          if (i <= j) {
            var_mat[i, j] = vardir[h, k]
            k = k + 1
          }
        }
      }
      var_mat = forceSymmetric(var_mat)
      dt[[h]] = as.matrix(var_mat)
    }
    return(dt)
  }
  eblup_inside = function(r, n, samevar, y, X, R, MAXITER = 100,
                          PRECISION = 1e-04) {
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
      g1 = diag(Gn %*% Vinv %*% R)
      g2 = diag(R %*% Vinv %*% X %*% Q %*% t(X) %*% t(R %*%
                                                        Vinv))
      dg = Vinv - Gn %*% Vinv %*% Vinv
      gg3 = (dg %*% V %*% t(dg))/iF
      g3 = diag(gg3)
      mse = g1 + g2 + 2 * g3
      mse.df = data.frame(matrix(data = mse, nrow = n,
                                 ncol = r))
      names(mse.df) = y_names
    }
    else {
      Vu = apply(matrix(diag(R), nrow = n, ncol = r), 2,
                 median)
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
      FI = solve(iF)
      g1 = diag(Gn %*% Vinv %*% R)
      g2 = diag(R %*% Vinv %*% X %*% Q %*% t(X) %*% t(R %*%
                                                        Vinv))
      dg = lapply(dV1, function(x) x %*% Vinv - Gn %*%
                    Vinv %*% x %*% Vinv)
      gg3 = list()
      for (i in 1:r) {
        for (j in 1:r) {
          gg3[[(i - 1) * r + j]] = FI[i, j] * (dg[[i]] %*%
                                                 V %*% t(dg[[j]]))
        }
      }
      g3 = diag(Reduce("+", gg3))
      mse = g1 + g2 + 2 * g3
      mse.df = data.frame(matrix(data = mse, nrow = n,
                                 ncol = r))
      names(mse.df) = y_names
    }
    result = list(eblup = NA, fit = list(estcoef = NA, refvar = NA),
                  mse = NA, mse_component = list(g1 = NA, g2 = NA,
                                                 g3 = NA))
    result$eblup = eblup
    result$fit$estcoef = coef
    result$fit$refvar = t(Vu)
    result$mse = mse.df
    result$mse_component$g1 = g1
    result$mse_component$g2 = g2
    result$mse_component$g3 = g3
    return(result)
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
    samevar = samevar
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
    samevar = samevar
  }
  y_names = sapply(formula, "[[", 2)
  eblup_first = eblup_inside(r = r, n = n, samevar = samevar,
                             y = y, X = X, R = R)
  beta = eblup_first$fit$estcoef[, 1]
  A = eblup_first$fit$refvar
  A_mat = kronecker(diag(as.vector(A)), diag(n))
  Vinv = solve(A_mat + R)
  XtVinv = t(Vinv %*% X)
  Q = solve(XtVinv %*% X)
  temp_vardir = temporary_vardir(vardir, n = n, r = r)
  mse_prasad = eblup_first$mse
  g1d = eblup_first$mse_component$g1
  g2d = eblup_first$mse_component$g2
  sumg1.pb = rep(0, n * r)
  sumg2.pb = rep(0, n * r)
  sumg3.pb = rep(0, n * r)
  boot <- 1
  while (boot <= B) {
    u.boot = mvrnorm(n = n, mu = rep(0, r), Sigma = (diag(as.vector(A),
                                                          nrow = r, ncol = r)))
    theta.boot = X %*% beta + as.vector(u.boot)
    e.boot = matrix(0, nrow = n, ncol = r)
    for (i in 1:n) {
      e.boot[i, ] = mvrnorm(1, mu = rep(0, r), Sigma = (temp_vardir[[i]]))
    }
    direct.boot = theta.boot + as.vector(e.boot)
    direct.boot.mat = matrix(direct.boot, nrow = n, ncol = r)
    resultEBLUP = eblup_inside(r = r, n = n, samevar = samevar,
                               y = direct.boot, X = X, R = R)
    sigma2.simula = resultEBLUP$fit$refvar
    beta.simula = resultEBLUP$fit$estcoef[, 1]
    mse.simula = resultEBLUP$mse
    Gn.simula = kronecker(diag(as.vector(sigma2.simula)),
                          diag(n))
    Vinv.simula = as.matrix(solve(Gn.simula + R))
    Xbeta.simula = X %*% beta.simula
    XtVi.simula = t(Vinv.simula %*% X)
    Q.simula = solve(XtVi.simula %*% X)
    thetaEBLUP.boot1 = Xbeta.simula + Gn.simula %*% Vinv.simula %*%
      (direct.boot - Xbeta.simula)
    thetaEBLUP.boot1.mat = matrix(thetaEBLUP.boot1, nrow = n,
                                  ncol = r)
    thetaOPTIMUM.boot1.mat = matrix(0, nrow = n, ncol = r)
    thetaALFA.boot1 = matrix(0, nrow = n, ncol = r)
    thetaLAMBDAboot1 = matrix(0, nrow = n, ncol = r)
    for (i in 1:r) {
      thetaALFA.boot1[, i] = (sum(W[, i] * direct.boot.mat[, i])
                              - sum(W[, i] * thetaEBLUP.boot1.mat[, i]))
    }
    for (i in 1:r) {
      for (j in 1:n) {
        thetaLAMBDAboot1[j , i] = (mse.simula[j , i] * W[j , i]) / sum(mse.simula[, i] * W[, i]^2)
      }
    }
    for (i in 1:r) {
      for (j in 1:n) {
        thetaOPTIMUM.boot1.mat[j , i] = thetaEBLUP.boot1.mat[j , i] +
          (thetaLAMBDAboot1[j , i] * thetaALFA.boot1[j , i])
      }
    }
    Bstim.eblup = solve(XtVinv %*% X) %*% XtVinv %*% direct.boot
    Xbeta.eblup = X %*% Bstim.eblup
    thetaEBLUP.boot2 = Xbeta.eblup + A_mat %*% Vinv %*% (direct.boot -
                                                           Xbeta.eblup)
    thetaEBLUP.boot2.mat = matrix(thetaEBLUP.boot2, nrow = n,
                                  ncol = r)
    thetaOPTIMUM.boot2.mat = matrix(0, nrow = n, ncol = r)
    thetaALFA.boot2 = matrix(0, nrow = n, ncol = r)
    thetaLAMBDAboot2 = matrix(0, nrow = n, ncol = r)
    for (i in 1:r) {
      thetaALFA.boot2[, i] = (sum(W[, i] * direct.boot.mat[, i])
                              - sum(W[, i] * thetaEBLUP.boot2.mat[, i]))
    }
    for (i in 1:r) {
      for (j in 1:n) {
        thetaLAMBDAboot2[j , i] = (mse.simula[j , i] * W[j , i]) / sum(mse.simula[, i] * W[, i]^2)
      }
    }
    for (i in 1:r) {
      for (j in 1:n) {
        thetaOPTIMUM.boot2.mat[j , i] = thetaEBLUP.boot2.mat[j , i] +
          (thetaLAMBDAboot2[j , i] * thetaALFA.boot2[j , i])
      }
    }
    g1boot = diag(Gn.simula %*% Vinv.simula %*% R)
    g2boot = diag(R %*% Vinv.simula %*% X %*% Q.simula %*%
                    t(X) %*% t(R %*% Vinv.simula))
    g3boot = (thetaOPTIMUM.boot1.mat - thetaOPTIMUM.boot2.mat)^2
    sumg1.pb = sumg1.pb + g1boot
    sumg2.pb = sumg2.pb + g2boot
    sumg3.pb = sumg3.pb + as.vector(g3boot)
    boot = boot + 1
  }
  g1.pb = sumg1.pb/B
  g2.pb = sumg2.pb/B
  g3.pb = sumg3.pb/B
  msebootoptimum = 2 * (g1d + g2d) - g1.pb - g2.pb + g3.pb
  msebootoptimum.df = data.frame(matrix(data = msebootoptimum,
                                        nrow = n, ncol = r))
  names(msebootoptimum.df) = y_names
  end_time <- Sys.time()
  running_time = end_time - start_time
  result1 = list(mse.eblup = NA, pbmse.eblupOB = NA, running.time = NA)
  result1$mse.eblup = mse_prasad
  result1$pbmse.eblupOB = msebootoptimum.df
  result1$running.time = running_time
  return(result1)
}
