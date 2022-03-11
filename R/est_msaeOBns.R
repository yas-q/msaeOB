#' @title EBLUPs Optimum Benchmarking for Non Sampled Area based on a Multivariate Fay Herriot (Model 1)
#'
#' @description This function gives EBLUPs optimum benchmarking for non sampled area based on multivariate Fay-Herriot (Model 1)
#'
#' @param formula an object of class list of formula describe the fitted models
#' @param vardir matrix containing sampling variances of direct estimators. The order is: \code{var1, cov12, ..., cov1r, var2, cov23, ..., cov2r, ..., cov(r-1)(r), var(r)}
#' @param weight matrix containing proportion of units in small areas. The order is: \code{w1, w2, ..., w(r)}
#' @param cluster matrix containing cluster of auxiliary variables. The order is: \code{c1, c2, ..., c(r)}
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
#' \item{fit}{a list contining following objects:}
#' \itemize{
#'   \item method : fitting method, named "REML"
#'   \item convergence : logical value of convergence of Fisher Scoring
#'   \item iterations : number of iterations of Fisher Scoring algorithm
#'   \item estcoef : a data frame containing estimated model coefficients (\code{beta, std. error, t value, p-value})
#'   \item refvar : estimated random effect variance
#' }
#' \item{random.effect}{a data frame containing values of random effect estimators}
#' \item{agregation}{a data frame containing agregation of direct, EBLUP, and optimum benchmark estimation}
#' @export est_msaeOBns
#'
#' @import abind
#' @importFrom magic adiag
#' @importFrom Matrix forceSymmetric
#' @importFrom stats model.frame na.omit model.matrix median pnorm rnorm
#' @importFrom MASS mvrnorm
#'
#' @examples
#' ## load dataset
#' data(datamsaeOBns)
#'
#' # Compute EBLUP and Optimum Benchmark using auxiliary variables X1 and X2 for each dependent variable
#'
#' ## Using parameter 'data'
#' Fo = list(f1 = Y1 ~ X1 + X2,
#'           f2 = Y2 ~ X1 + X2,
#'           f3 = Y3 ~ X1 + X2)
#' vardir = c("v1", "v12", "v13", "v2", "v23", "v3")
#' weight = c("w1", "w2", "w3")
#' cluster = c("c1", "c2", "c3")
#'
#' est_msae = est_msaeOBns(Fo, vardir, weight, cluster, data = datamsaeOBns)
#'
#' ## Without parameter 'data'
#' Fo = list(f1 = datamsaeOBns$Y1 ~ datamsaeOBns$X1 + datamsaeOBns$X2,
#'           f2 = datamsaeOBns$Y2 ~ datamsaeOBns$X1 + datamsaeOBns$X2,
#'           f3 = datamsaeOBns$Y3 ~ datamsaeOBns$X1 + datamsaeOBns$X2)
#' vardir = datamsaeOBns[, c("v1", "v12", "v13", "v2", "v23", "v3")]
#' weight = datamsaeOBns[, c("w1", "w2", "w3")]
#' cluster = datamsaeOBns[, c("c1", "c2", "c3")]
#'
#' est_msae = est_msaeOBns(Fo, vardir, weight, cluster)
#'
#' ## Return
#' est_msae$eblup$est.eblupOB # to see the Optimum Benchmark estimators
#'
est_msaeOBns<-function (formula, vardir, weight, cluster, samevar = FALSE,
                        MAXITER = 100, PRECISION = 1e-04, data)
{
  r = length(formula)
  if (r <= 1)
    stop("You should use est_saeOBns() for univariate")
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
    if (!all(cluster %in% names(data)))
      stop("Object cluster is not appropriate with data.")
    if (length(cluster) != r)
      stop("Length of cluster is not appropriate with data. The length must be ",
           r)
    vardir = data[, vardir]
    cluster = as.matrix(data[, cluster])
  }
  else {
    formuladata = lapply(formula, function(x) model.frame(x,
                                                          na.action = na.omit))
    y = unlist(lapply(formula, function(x) model.frame(x,
                                                       na.action = na.omit)[[1]]))
    X = Reduce(adiag, lapply(formula, function(x) model.matrix(x)))
    W = as.matrix(weight)
    n = length(y)/r
    if ((dim(vardir)[1] != n) || (dim(vardir)[2] != sum(1:r)))
      stop("Object vardir is not appropriate with data. It must be ",
           n, " x ", sum(1:r), " matrix.")
    if (any(is.na(weight)))
      stop("Object weight contains NA values.")
    if ((dim(weight)[1] != n) || (dim(weight)[2] != r))
      stop("Object weight is not appropriate with data. It must be ",
           n, " x ", r, " matrix.")
    if (any(is.na(cluster)))
      stop("Object cluster contains NA values.")
    if ((dim(cluster)[1] != n) || (dim(cluster)[2] != r))
      stop("Object cluster is not appropriate with data. It must be ",
           n, " x ", r, " matrix.")
    cluster = as.matrix(cluster)
  }
  y.matrix = matrix(as.vector(y), n, r)
  y.matrix[y.matrix == 0] = NA
  indexns = unique(which(is.na(y.matrix), arr.ind = TRUE)[,
                                                          1])
  indexs = c(1:n)[-indexns]
  n.s = length(indexs)
  n.ns = n - n.s
  y.s = y.matrix[-indexns, ]
  y.s.vec = as.vector(y.s)
  corr = cor(y.s)
  hasil_cor = c()
  for (i in 1:r) {
    for (j in 1:r) {
      if (i < j && i != j) {
        hasil_cor = c(hasil_cor, corr[i, j])
      }
    }
  }
  vardir[indexns, ] = NA
  index.vardir = sum(1:r)
  for (i in 2:r) {
    index.vardir[i] = index.vardir[1] - sum(2:i)
  }
  index.vardir = sort(index.vardir)
  index.covardir = c(1:sum(1:r))[-index.vardir]
  varians.direct = vardir[, index.vardir]
  covarians.direct = vardir[, -index.vardir]
  for (i in 1:r) {
    df.vardir = data.frame(cluster[, i], varians.direct[,
                                                        i])
    for (j in 1:n) {
      if (df.vardir[j, 2] %in% NA) {
        df.vardir[j, 2] = mean(df.vardir[df.vardir[,
                                                   1] == df.vardir[j, 1], 2], na.rm = TRUE)
      }
    }
    varians.direct[, i] = df.vardir[, 2]
  }
  index.cor = 0
  if (r > 2) {
    for (p in 1:r) {
      for (q in 1:r) {
        if (p < q) {
          index.cor = index.cor + 1
          for (col in 1:n) {
            if (covarians.direct[col, index.cor] %in%
                NA) {
              covarians.direct[col, index.cor] = hasil_cor[index.cor] *
                sqrt(varians.direct[col, p] * varians.direct[col,
                                                             q])
            }
          }
        }
      }
    }
  }
  else {
    index.cor = index.cor + 1
    for (col in 1:n) {
      if (covarians.direct[col] %in% NA) {
        covarians.direct[col] = hasil_cor[index.cor] *
          sqrt(varians.direct[col, 1] * varians.direct[col,
                                                       2])
      }
    }
  }
  varians.direct = rbind(index.vardir, varians.direct)
  if (r > 2) {
    covarians.direct = rbind(index.covardir, covarians.direct)
  }
  else {
    covarians.direct = c(index.covardir, covarians.direct)
  }
  vardir.full = cbind(varians.direct, covarians.direct)
  vardir.full = vardir.full[, order(vardir.full[1, ])]
  vardir.full = vardir.full[-1, ]
  R.ns = R_function(vardir.full[indexns, ], n.ns, r)
  vardir.s = vardir[-indexns, ]
  R = R_function(vardir.s, n.s, r)
  W.s = W[-indexns, ]
  W.s = prop.table(W.s, 2)
  Xindexns = c()
  for (i in 1:r) {
    Xindexns = c(Xindexns, indexns + rep(n, times = n.ns) *
                   (i - 1))
  }
  X.s = X[-Xindexns, ]
  X.ns = X[Xindexns, ]
  y_names = sapply(formula, "[[", 2)
  Ir = diag(r)
  In = diag(n.s)
  dV = list()
  dV1 = list()
  for (i in 1:r) {
    dV[[i]] = matrix(0, nrow = r, ncol = r)
    dV[[i]][i, i] = 1
    dV1[[i]] = kronecker(dV[[i]], In)
  }
  convergence = TRUE
  vardir.full = cbind(varians.direct, covarians.direct)
  vardir.full = vardir.full[, order(vardir.full[1, ])]
  vardir.full = vardir.full[-1, ]
  R.ns = R_function(vardir.full[indexns, ], n.ns, r)
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
      XstVinv = t(Vinv %*% X.s)
      Q = solve(XstVinv %*% X.s)
      P = Vinv - t(XstVinv) %*% Q %*% XstVinv
      Py = P %*% y.s.vec
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
    XstVinv = t(Vinv %*% X.s)
    Q = solve(XstVinv %*% X.s)
    P = Vinv - t(XstVinv) %*% Q %*% XstVinv
    Py = P %*% y.s.vec
    beta = Q %*% XstVinv %*% y.s.vec
    res = y.s.vec - X.s %*% beta
    random.effect = data.frame(matrix(Gn %*% Vinv %*% res,
                                      n.s, r))
    names(random.effect) = y_names
    se.b = sqrt(diag(Q))
    t.value = beta/se.b
    p.value = 2 * pnorm(abs(as.numeric(t.value)), lower.tail = FALSE)
    coef = as.matrix(cbind(beta, se.b, t.value, p.value))
    colnames(coef) = c("beta", "std. error",
                       "t value", "p-value")
    rownames(coef) = colnames(X)
    Gn.ns = kronecker(diag(Vu), diag(n.ns))
    Vinv.ns = solve(Gn.ns + R.ns)
    Q.ns = ginv(t(Vinv.ns %*% X.ns) %*% X.ns)
    g1.s = diag(Gn %*% Vinv %*% R)
    g2.s = diag(R %*% Vinv %*% X.s %*% Q %*% t(X.s) %*%
                  t(R %*% Vinv))
    g1.ns = diag(Gn.ns %*% Vinv.ns %*% R.ns)
    g2.ns = diag(R.ns %*% Vinv.ns %*% X.ns %*% Q.ns %*%
                   t(X.ns) %*% t(R.ns %*% Vinv.ns))
    g12.s = matrix(g1.s + g2.s, n.s, r)
    names(g12.s) = y_names
    g12.s = data.frame(index = indexs, g12.s)
    g12.ns = matrix(g1.ns + g2.ns, n.ns, r)
    names(g12.ns) = y_names
    g12.ns = data.frame(index = indexns, g12.ns)
    g12 = rbind(g12.s, g12.ns)
    g12 = g12[order(g12$index), ]
    rownames(g12) = g12$index
    g12 = g12[, -1]
    dg = Vinv - Gn %*% Vinv %*% Vinv
    gg3 = (dg %*% V %*% t(dg))/iF
    g3 = diag(gg3)
    g3.matrix = matrix(g3, n.s, r)
    g3.s = data.frame(index = indexs, g3.matrix)
    g3.ns = data.frame(index = indexns, matrix(NA, n.ns,
                                               r))
    names(g3.ns) = names(g3.s)
    g3.all = rbind(g3.s, g3.ns)
    g3.all = g3.all[order(g3.all$index), ]
    rownames(g3.all) = g3.all$index
    g3.all = g3.all[, -1]
    g3.full = matrix(NA, nrow = n, ncol = r)
    colnames(g3.full) = y_names
    for (i in 1:r) {
      df = data.frame(cluster[, i], g3.all[, i])
      for (j in 1:n) {
        if (df[j, 2] %in% NA) {
          df[j, 2] = mean(df[df[, 1] == df[j, 1], 2],
                          na.rm = TRUE)
        }
      }
      g3.full[, i] = df[, 2]
    }
    names(g3.full) = y_names
    mse = as.data.frame(g12 + 2 * g3.full)
    names(mse) = y_names
  }
  else {
    Vu = apply(matrix(diag(R), nrow = n.s, ncol = r), 2,
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
      XstVinv = t(Vinv %*% X.s)
      Q = solve(XstVinv %*% X.s)
      P = Vinv - t(XstVinv) %*% Q %*% XstVinv
      Py = P %*% y.s.vec
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
      Gr = Vu
    }
    else {
      Gr = diag(as.vector(Vu))
    }
    Gn = kronecker(Gr, In)
    V = as.matrix(Gn + R)
    Vinv = solve(V)
    XtVinv = t(Vinv %*% X.s)
    Q = as.matrix(solve(XtVinv %*% X.s))
    P = Vinv - t(XtVinv) %*% Q %*% XtVinv
    Py = P %*% y.s.vec
    beta = Q %*% XtVinv %*% y.s.vec
    res = y.s.vec - X.s %*% beta
    random.effect = data.frame(matrix(Gn %*% Vinv %*% res,
                                      n.s, r))
    names(random.effect) = y_names
    se.b = sqrt(diag(Q))
    t.value = beta/se.b
    p.value = 2 * pnorm(abs(as.numeric(t.value)), lower.tail = FALSE)
    coef = as.matrix(cbind(beta, se.b, t.value, p.value))
    colnames(coef) = c("beta", "std. error",
                       "t value", "p-value")
    rownames(coef) = colnames(X)
    FI = solve(iF)
    Gn.ns = kronecker(Gr, diag(n.ns))
    Vinv.ns = solve(Gn.ns + R.ns)
    XtVinv.ns = t(Vinv.ns %*% X.ns)
    Q.ns = ginv(XtVinv.ns %*% X.ns)
    g1.s = diag(Gn %*% Vinv %*% R)
    g2.s = diag(R %*% Vinv %*% X.s %*% Q %*% t(X.s) %*%
                  t(R %*% Vinv))
    g1.ns = diag(Gn.ns %*% Vinv.ns %*% R.ns)
    g2.ns = diag(R.ns %*% Vinv.ns %*% X.ns %*% Q.ns %*%
                   t(X.ns) %*% t(R.ns %*% Vinv.ns))
    g12.s = matrix(g1.s + g2.s, n.s, r)
    names(g12.s) = y_names
    g12.s = data.frame(index = indexs, g12.s)
    g12.ns = matrix(g1.ns + g2.ns, n.ns, r)
    names(g12.ns) = y_names
    g12.ns = data.frame(index = indexns, g12.ns)
    g12 = rbind(g12.s, g12.ns)
    g12 = g12[order(g12$index), ]
    rownames(g12) = g12$index
    g12 = g12[, -1]
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
    g3.matrix = matrix(g3, n.s, r)
    g3.s = data.frame(index = indexs, g3.matrix)
    g3.ns = data.frame(index = indexns, matrix(NA, n.ns,
                                               r))
    names(g3.ns) = names(g3.s)
    g3.all = rbind(g3.s, g3.ns)
    g3.all = g3.all[order(g3.all$index), ]
    rownames(g3.all) = g3.all$index
    g3.all = g3.all[, -1]
    g3.full = matrix(NA, nrow = n, ncol = r)
    colnames(g3.full) = y_names
    for (i in 1:r) {
      df = data.frame(cluster[, i], g3.all[, i])
      for (j in 1:n) {
        if (df[j, 2] %in% NA) {
          df[j, 2] = mean(df[df[, 1] == df[j, 1], 2],
                          na.rm = TRUE)
        }
      }
      g3.full[, i] = df[, 2]
    }
    names(g3.full) = y_names
    mse = as.data.frame(g12 + 2 * g3.full)
    names(mse) = y_names
  }
  random.effect = data.frame(index = indexs, random.effect)
  random.effect.ns = data.frame(index = indexns, matrix(NA,
                                                        n.ns, r))
  names(random.effect.ns) = names(random.effect)
  random.effect = rbind(random.effect, random.effect.ns)
  random.effect = random.effect[order(random.effect$index),
  ]
  rownames(random.effect) = random.effect$index
  random.effect = random.effect[, -1]
  random.effect.full = matrix(NA, nrow = n, ncol = r)
  colnames(random.effect.full) = y_names
  for (i in 1:r) {
    df = data.frame(cluster[, i], random.effect[, i])
    for (j in 1:n) {
      if (df[j, 2] %in% NA) {
        df[j, 2] = mean(df[df[, 1] == df[j, 1], 2], na.rm = TRUE)
      }
    }
    random.effect.full[, i] = df[, 2]
  }
  Xbeta = X %*% beta
  eblup = matrix(Xbeta, nrow = n, ncol = r) + random.effect.full
  eblup = as.data.frame(eblup)
  names(eblup) = y_names
  random.effect.full = as.data.frame(random.effect.full)
  y.s.mat = matrix(y.s, nrow = n.s, ncol = r)
  eblup.mat = as.matrix(eblup)
  alfa = matrix(0, nrow = n, ncol = r)
  lambda = matrix(0, nrow = n, ncol = r)
  eblup.optimum = matrix(0, nrow = n, ncol = r)
  for (i in 1:r) {
    alfa[, i] = (sum(W.s[, i] * y.s.mat[, i]) - sum(W[, i] * eblup.mat[, i]))
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
  agregation.direct = diag(t(W.s) %*% y.s)
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
  result$random.effect = random.effect.full
  result$agregation = agregation
}
