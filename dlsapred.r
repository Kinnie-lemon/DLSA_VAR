dlsa.pred <- function (model, h = 1, orig = 0, Out.level = FALSE, output = TRUE,
          Phi = model$theta, Ph0 = model$Ph0,cnst = TRUE) 
{
  x = model$data
  Phi = model$theta
  sig = model$sigma
  Ph0 = model$Ph0
  p = as.numeric(model$optlag[1])
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
  Ginv = solve(G)
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