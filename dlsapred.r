dlsa.pred <- function (model, h = 1, orig = 0, Out.level = FALSE, output = TRUE,
                       Phi = model$theta,Ph0 = model$Ph0,cnst = TRUE) 
{
  x = model$data
  sig = model$sigma
  p = as.numeric(model$optlag)
  np = dim(Phi)[2]
  k = dim(x)[2]
  nT = dim(x)[1]
  k = dim(x)[2]
  if (orig <= 0) 
    orig = nT
  if (orig > nT) 
    orig = nT
  psi = VARpsi(Phi, h)$psi
  beta = t(Phi)
  if (length(Ph0) < 1) 
    Ph0 = rep(0, k)
  if (p > orig) {
    cat("Too few data points to produce forecasts", 
        "\n")
  }
  pred = NULL
  se = NULL
  MSE = NULL
  mse = NULL
  px = as.matrix(x[1:orig, ])
  Past = px[orig, ]
  if (p > 1) {
    for (j in 1:(p - 1)) {
      Past = c(Past, px[(orig - j), ])
    }
  }
  cat("orig ", orig, "\n")
  ne = orig - p
  xmtx = NULL
  P = NULL
  if (cnst) 
    xmtx = rep(1, ne)
  xmtx = cbind(xmtx, x[p:(orig - 1), ])
  ist = p + 1
  if (p > 1) {
    for (j in 2:p) {
      xmtx = cbind(xmtx, x[(ist - j):(orig - j), ])
    }
  }
  xmtx = as.matrix(xmtx)
  G = t(xmtx) %*% xmtx/ne
  Ginv = solve(G,tol = 2e-35)
  P = Phi
  vv = Ph0
  if (p > 1) {
    II = diag(rep(1, k * (p - 1)))
    II = cbind(II, matrix(0, (p - 1) * k, k))
    P = rbind(P, II)
    vv = c(vv, rep(0, (p - 1) * k))
  }
  if (cnst) {
    c1 = c(1, rep(0, np))
    P = cbind(vv, P)
    P = rbind(c1, P)
  }
  Sig = sig
  n1 = dim(P)[2]
  MSE = (n1/orig) * sig
  for (j in 1:h) {
    tmp = Ph0 + matrix(Past, 1, np) %*% beta
    px = rbind(px, tmp)
    if (np > k) {
      Past = c(tmp, Past[1:(np - k)])
    }
    else {
      Past = tmp
    }
    if (j > 1) {
      idx = (j - 1) * k
      wk = psi[, (idx + 1):(idx + k)]
      Sig = Sig + wk %*% sig %*% t(wk)
    }
    if (j > 1) {
      for (ii in 0:(j - 1)) {
        psii = diag(rep(1, k))
        if (ii > 0) {
          idx = ii * k
          psii = psi[, (idx + 1):(idx + k)]
        }
        P1 = P^(j - 1 - ii) %*% Ginv
        for (jj in 0:(j - 1)) {
          psij = diag(rep(1, k))
          if (jj > 0) {
            jdx = jj * k
            psij = psi[, (jdx + 1):(jdx + k)]
          }
          P2 = P^(j - 1 - jj) %*% G
          k1 = sum(diag(P1 %*% P2))
          MSE = (k1/orig) * psii %*% sig %*% t(psij)
        }
      }
    }
    se = rbind(se, sqrt(diag(Sig)))
    if (Out.level) {
      cat("Covariance matrix of forecast errors at horizon: ", 
          j, "\n")
      print(Sig)
      cat("Omega matrix at horizon: ", j, "\n")
      print(MSE)
    }
    MSE = MSE + Sig
    mse = rbind(mse, sqrt(diag(MSE)))
  }
  if (output) {
    cat("Forecasts at origin: ", orig, "\n")
    print(px[(orig + 1):(orig + h), ], digits = 4)
    cat("Standard Errors of predictions: ", "\n")
    print(se[1:h, ], digits = 4)
    pred = px[(orig + 1):(orig + h), ]
    cat("Root mean square errors of predictions: ", 
        "\n")
    print(mse[1:h, ], digits = 4)
  }
  if (orig < nT) {
    cat("Observations, predicted values,     errors, and MSE", 
        "\n")
    tmp = NULL
    jend = min(nT, (orig + h))
    for (t in (orig + 1):jend) {
      case = c(t, x[t, ], px[t, ], x[t, ] - px[t, ])
      tmp = rbind(tmp, case)
    }
    colnames(tmp) <- c("time", rep("obs", k), 
                       rep("fcst", k), rep("err", k))
    idx = c(1)
    for (j in 1:k) {
      idx = c(idx, c(0, 1, 2) * k + j + 1)
    }
    tmp = tmp[, idx]
    print(round(tmp, 4))
  }
  VARpred <- list(pred = pred, se.err = se, rmse = mse)
}


##VAR模型精确度改进
VAR.v <-function (x, p = 1, output = T, include.mean = T, fixed = NULL) 
{
  if (!is.matrix(x)) 
    x = as.matrix(x)
  Tn = dim(x)[1]
  k = dim(x)[2]
  if (p < 1) 
    p = 1
  idm = k * p
  ne = Tn - p
  ist = p + 1
  y = x[ist:Tn, ]
  if (include.mean) {
    idm = idm + 1
    xmtx = cbind(rep(1, ne), x[p:(Tn - 1), ])
  }
  else {
    xmtx = x[p:(Tn - 1), ]
  }
  if (p > 1) {
    for (i in 2:p) {
      xmtx = cbind(xmtx, x[(ist - i):(Tn - i), ])
    }
  }
  ndim = ncol(xmtx)
  if (length(fixed) == 0) {
    paridx = matrix(1, ndim, k)
  }
  else {
    paridx = fixed
  }
  res = NULL
  beta = matrix(0, ndim, k)
  sdbeta = matrix(0, ndim, k)
  npar = 0
  for (i in 1:k) {
    idx = c(1:ndim)[paridx[, i] == 1]
    resi = y[, i]
    if (length(idx) > 0) {
      xm = as.matrix(xmtx[, idx])
      npar = npar + dim(xm)[2]
      xpx = t(xm) %*% xm
      xpxinv = solve(xpx,tol = 2e-35)
      xpy = t(xm) %*% as.matrix(y[, i], ne, 1)
      betai = xpxinv %*% xpy
      beta[idx, i] = betai
      resi = y[, i] - xm %*% betai
      nee = dim(xm)[2]
      sse = sum(resi * resi)/(Tn - p - nee)
      dd = diag(xpxinv)
      sdbeta[idx, i] = sqrt(dd * sse)
    }
    res = cbind(res, resi)
  }
  sse = t(res) %*% res/(Tn - p)
  aic = 0
  bic = 0
  hq = 0
  Phi = NULL
  Ph0 = NULL
  jst = 0
  if (include.mean) {
    Ph0 = beta[1, ]
    se = sdbeta[1, ]
    if (output) {
      cat("Constant term:", "\n")
      cat("Estimates: ", Ph0, "\n")
      cat("Std.Error: ", se, "\n")
    }
    jst = 1
  }
  if (include.mean) {
    for (i in 1:k) {
      if (abs(Ph0[i]) > 1e-08) 
        npar = npar - 1
    }
  }
  if (output) 
    cat("AR coefficient matrix", "\n")
  for (i in 1:p) {
    phi = t(beta[(jst + 1):(jst + k), ])
    se = t(sdbeta[(jst + 1):(jst + k), ])
    if (output) {
      cat("AR(", i, ")-matrix", "\n")
      print(phi, digits = 3)
      cat("standard error", "\n")
      print(se, digits = 3)
    }
    jst = jst + k
    Phi = cbind(Phi, phi)
  }
  if (output) {
    cat(" ", "\n")
    cat("Residuals cov-mtx:", "\n")
    print(sse)
    cat(" ", "\n")
  }
  dd = det(sse)
  d1 = log(dd)
  aic = d1 + (2 * npar)/Tn
  bic = d1 + log(Tn) * npar/Tn
  hq = d1 + 2 * log(log(Tn)) * npar/Tn
  if (output) {
    cat("det(SSE) = ", dd, "\n")
    cat("AIC = ", aic, "\n")
    cat("BIC = ", bic, "\n")
    cat("HQ  = ", hq, "\n")
  }
  VAR <- list(data = x, cnst = include.mean, order = p, coef = beta, 
              aic = aic, bic = bic, hq = hq, residuals = res, secoef = sdbeta, 
              Sigma = sse, Phi = Phi, Ph0 = Ph0, fixed = fixed)
}

##VARPRED 精度修正
VAR.pred <-function (model, h = 1, orig = 0, Out.level = FALSE, output = TRUE) 
{
  x = model$data
  Phi = model$Phi
  sig = model$Sigma
  Ph0 = model$Ph0
  p = model$order
  cnst = model$cnst
  np = dim(Phi)[2]
  k = dim(x)[2]
  nT = dim(x)[1]
  k = dim(x)[2]
  if (orig <= 0) 
    orig = nT
  if (orig > nT) 
    orig = nT
  psi = VARpsi(Phi, h)$psi
  beta = t(Phi)
  if (length(Ph0) < 1) 
    Ph0 = rep(0, k)
  if (p > orig) {
    cat("Too few data points to produce forecasts", 
        "\n")
  }
  pred = NULL
  se = NULL
  MSE = NULL
  mse = NULL
  px = as.matrix(x[1:orig, ])
  Past = px[orig, ]
  if (p > 1) {
    for (j in 1:(p - 1)) {
      Past = c(Past, px[(orig - j), ])
    }
  }
  cat("orig ", orig, "\n")
  ne = orig - p
  xmtx = NULL
  P = NULL
  if (cnst) 
    xmtx = rep(1, ne)
  xmtx = cbind(xmtx, x[p:(orig - 1), ])
  ist = p + 1
  if (p > 1) {
    for (j in 2:p) {
      xmtx = cbind(xmtx, x[(ist - j):(orig - j), ])
    }
  }
  xmtx = as.matrix(xmtx)
  G = t(xmtx) %*% xmtx/ne
  Ginv = solve(G,tol = 2e-35)
  P = Phi
  vv = Ph0
  if (p > 1) {
    II = diag(rep(1, k * (p - 1)))
    II = cbind(II, matrix(0, (p - 1) * k, k))
    P = rbind(P, II)
    vv = c(vv, rep(0, (p - 1) * k))
  }
  if (cnst) {
    c1 = c(1, rep(0, np))
    P = cbind(vv, P)
    P = rbind(c1, P)
  }
  Sig = sig
  n1 = dim(P)[2]
  MSE = (n1/orig) * sig
  for (j in 1:h) {
    tmp = Ph0 + matrix(Past, 1, np) %*% beta
    px = rbind(px, tmp)
    if (np > k) {
      Past = c(tmp, Past[1:(np - k)])
    }
    else {
      Past = tmp
    }
    if (j > 1) {
      idx = (j - 1) * k
      wk = psi[, (idx + 1):(idx + k)]
      Sig = Sig + wk %*% sig %*% t(wk)
    }
    if (j > 1) {
      for (ii in 0:(j - 1)) {
        psii = diag(rep(1, k))
        if (ii > 0) {
          idx = ii * k
          psii = psi[, (idx + 1):(idx + k)]
        }
        P1 = P^(j - 1 - ii) %*% Ginv
        for (jj in 0:(j - 1)) {
          psij = diag(rep(1, k))
          if (jj > 0) {
            jdx = jj * k
            psij = psi[, (jdx + 1):(jdx + k)]
          }
          P2 = P^(j - 1 - jj) %*% G
          k1 = sum(diag(P1 %*% P2))
          MSE = (k1/orig) * psii %*% sig %*% t(psij)
        }
      }
    }
    se = rbind(se, sqrt(diag(Sig)))
    if (Out.level) {
      cat("Covariance matrix of forecast errors at horizon: ", 
          j, "\n")
      print(Sig)
      cat("Omega matrix at horizon: ", j, "\n")
      print(MSE)
    }
    MSE = MSE + Sig
    mse = rbind(mse, sqrt(diag(MSE)))
  }
  if (output) {
    cat("Forecasts at origin: ", orig, "\n")
    print(px[(orig + 1):(orig + h), ], digits = 4)
    cat("Standard Errors of predictions: ", "\n")
    print(se[1:h, ], digits = 4)
    pred = px[(orig + 1):(orig + h), ]
    cat("Root mean square errors of predictions: ", 
        "\n")
    print(mse[1:h, ], digits = 4)
  }
  if (orig < nT) {
    cat("Observations, predicted values,     errors, and MSE", 
        "\n")
    tmp = NULL
    jend = min(nT, (orig + h))
    for (t in (orig + 1):jend) {
      case = c(t, x[t, ], px[t, ], x[t, ] - px[t, ])
      tmp = rbind(tmp, case)
    }
    colnames(tmp) <- c("time", rep("obs", k), 
                       rep("fcst", k), rep("err", k))
    idx = c(1)
    for (j in 1:k) {
      idx = c(idx, c(0, 1, 2) * k + j + 1)
    }
    tmp = tmp[, idx]
    print(round(tmp, 4))
  }
  VARpred <- list(pred = pred, se.err = se, rmse = mse)
}

