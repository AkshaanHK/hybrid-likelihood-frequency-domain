## =============================================================================
## 01_core_estimators.R
## Core frequency-domain building blocks:
##   - periodogram / DFT helpers
##   - spectral density (ARMA)
##   - Whittle, FDEL, Hybrid objectives and estimators
## =============================================================================

# ---- Periodogram & Fourier frequencies ----------------------------------------

#' Compute periodogram ordinates at positive Fourier frequencies
#' @param y numeric vector (time series)
#' @return numeric vector of length floor((n-1)/2)
per <- function(y) {
  n   <- length(y)
  per <- (Mod(fft(y))^2 / (2 * pi * n))[1:(n %/% 2 + 1)]
  per[2:((n + 1) %/% 2)]
}

#' Fourier frequencies matching the periodogram output of per()
#' @param y numeric vector (time series)
#' @return numeric vector omega_j = 2*pi*j/n, j=1,...,floor((n-1)/2)
freq <- function(y) {
  n     <- length(y)
  mhalfm <- (n - 1) %/% 2L
  2 * pi / n * (1:mhalfm)
}

# ---- Spectral density ---------------------------------------------------------

#' ARMA spectral density  f(lambda; eta) = (1/2pi) * |theta(e^{-i*lambda})|^2 /
#'                                                    |phi(e^{-i*lambda})|^2
#' @param eta  numeric vector: AR coefficients first, then MA coefficients
#' @param p    integer, AR order
#' @param q    integer, MA order
#' @param x    numeric vector of frequencies (in radians)
spectral <- function(eta = c(), p = 0, q = 0, x) {
  if (p > 0) {
    phi  <- eta[1:p]
    px   <- outer(x, 1:p)
    Rar  <- cos(px) %*% phi
    Iar  <- sin(px) %*% phi
    far  <- (1 - Rar)^2 + Iar^2
  } else {
    far <- 1
  }
  if (q > 0) {
    psi  <- eta[(p + 1):(p + q)]
    px   <- outer(x, 1:q)
    Rma  <- cos(px) %*% psi
    Ima  <- sin(px) %*% psi
    fma  <- (1 + Rma)^2 + Ima^2
  } else {
    fma <- 1
  }
  (1 / (2 * pi)) * (fma / far)
}

# ---- Score building blocks ----------------------------------------------------

#' d/d(beta) log f(omega; beta)  for AR(1)
dlnf_dbeta <- function(beta, omega) {
  denom <- 1 + beta^2 - 2 * beta * cos(omega)
  -2 * (beta - cos(omega)) / denom
}

#' d/d(phi) log f(omega; phi, theta)  for ARMA(1,1)
dlnf_dphi <- function(phi, omega) {
  denom <- 1 + phi^2 - 2 * phi * cos(omega)
  -2 * (phi - cos(omega)) / denom
}

#' d/d(theta) log f(omega; phi, theta)  for ARMA(1,1)
dlnf_dtheta <- function(theta, omega) {
  denom <- 1 + theta^2 + 2 * theta * cos(omega)
  2 * (theta + cos(omega)) / denom
}

#' Monti spectral estimating function for AR(1)
#'   psi_j(beta) = (I(lambda_j)/f(lambda_j; beta) - 1) * d/dbeta log f(lambda_j; beta)
psi_monti <- function(I, beta, omega) {
  g   <- spectral(eta = c(beta), p = 1, q = 0, x = omega)
  dln <- dlnf_dbeta(beta, omega)
  (I / g - 1) * dln
}

# ---- Lagrange multiplier solvers (FDEL inner problem) -------------------------

#' Solve  sum_j psi_j / (1 + xi * psi_j) = 0  for scalar xi  (AR(1) / ARFIMA)
solve_xi <- function(psi, tol = 1e-10) {
  if (all(psi >= 0) || all(psi <= 0)) return(NA)
  gfun <- function(xi) sum(psi / (1 + xi * psi))
  lo   <- -1 / max(psi) + 1e-12
  hi   <- -1 / min(psi) - 1e-12
  tryCatch(
    uniroot(gfun, lower = lo, upper = hi, tol = tol)$root,
    error = function(e) NA
  )
}

#' Solve  sum_j psi_j / (1 + xi' psi_j) = 0  for vector xi  (ARMA(1,1))
solve_xi_vec <- function(psi_mat, tol = 1e-8) {
  p    <- ncol(psi_mat)
  gfun <- function(xi) {
    xi    <- as.numeric(xi)
    denom <- 1 + drop(psi_mat %*% xi)
    if (any(denom <= 0)) return(rep(1e6, p))
    colSums(psi_mat / denom)
  }
  obj  <- function(xi) sum(gfun(xi)^2)
  res  <- optim(rep(0, p), obj, method = "BFGS",
                control = list(reltol = tol, maxit = 1000))
  xi   <- res$par
  denom <- 1 + drop(psi_mat %*% xi)
  if (any(denom <= 0)) return(rep(NA_real_, p))
  xi
}

# ---- AR(1) objectives ---------------------------------------------------------

#' Whittle negative log-likelihood for AR(1)
#' @param beta  scalar AR(1) coefficient
obj_Whittle <- function(beta, Iper, omega) {
  g <- spectral(eta = c(beta), p = 1, q = 0, x = omega)
  sum(log(g) + Iper / g)
}

#' FDEL objective for AR(1)   (returns Inf if infeasible)
obj_FDEL <- function(beta, Iper, omega) {
  psi   <- psi_monti(Iper, beta, omega)
  xi    <- solve_xi(psi)
  if (is.na(xi)) return(Inf)
  denom <- 1 + xi * psi
  if (any(denom <= 0)) return(Inf)
  2 * sum(log(denom))
}

#' Hybrid objective for AR(1):  alpha * Whittle + (1-alpha) * FDEL
#' @param alpha  tuning parameter in [0,1]
hybrid_obj <- function(beta, alpha, Iper, omega) {
  w <- obj_Whittle(beta, Iper, omega)
  e <- obj_FDEL(beta, Iper, omega)
  if (!is.finite(e)) return(Inf)
  alpha * w + (1 - alpha) * e
}

# ---- ARMA(1,1) objectives -----------------------------------------------------

#' Whittle negative log-likelihood for ARMA(1,1)
#' @param beta  c(phi, theta)
obj_Whittle_arma <- function(beta, Iper, omega) {
  g <- spectral(eta = beta, p = 1, q = 1, x = omega)
  sum(log(g) + Iper / g)
}

#' FDEL objective for ARMA(1,1)
obj_FDEL_arma <- function(beta, Iper, omega) {
  phi   <- beta[1]
  theta <- beta[2]
  g     <- spectral(eta = beta, p = 1, q = 1, x = omega)

  score_phi   <- dlnf_dphi(phi, omega)
  score_theta <- dlnf_dtheta(theta, omega)
  score_mat   <- cbind(score_phi, score_theta)   # m x 2

  w       <- as.numeric(Iper / g - 1)
  psi_mat <- sweep(score_mat, 1, w, FUN = "*")   # m x 2

  xi <- solve_xi_vec(psi_mat)
  if (any(!is.finite(xi))) return(Inf)

  denom <- 1 + drop(psi_mat %*% xi)
  if (any(denom <= 0)) return(Inf)
  2 * sum(log(denom))
}

#' Hybrid objective for ARMA(1,1)
hybrid_obj_arma <- function(beta, alpha, Iper, omega) {
  w <- obj_Whittle_arma(beta, Iper, omega)
  e <- obj_FDEL_arma(beta, Iper, omega)
  if (!is.finite(e)) return(Inf)
  alpha * w + (1 - alpha) * e
}

# ---- ARFIMA(1,d,0) objectives -------------------------------------------------

#' Spectral density of ARFIMA(1,d,0) as a function of phi (d known)
arfima_spectral_phi <- function(phi, d, omega) {
  ar_denom  <- 1 + phi^2 - 2 * phi * cos(omega)
  frac_part <- (2 * sin(omega / 2))^(-2 * d)
  (1 / (2 * pi)) * frac_part / ar_denom
}

#' Monti estimating function for ARFIMA(1,d,0) (phi only, d fixed)
psi_monti_arfima <- function(I, phi, omega, d) {
  g   <- arfima_spectral_phi(phi, d, omega)
  dln <- dlnf_dbeta(phi, omega)          # same AR part
  (I / g - 1) * dln
}

#' Whittle for ARFIMA(1,d,0)
obj_Whittle_arfima_phi <- function(phi, Iper, omega, d) {
  g <- arfima_spectral_phi(phi, d, omega)
  sum(log(g) + Iper / g)
}

#' FDEL for ARFIMA(1,d,0)
obj_FDEL_arfima_phi <- function(phi, Iper, omega, d) {
  psi   <- psi_monti_arfima(Iper, phi, omega, d)
  xi    <- solve_xi(psi)
  if (is.na(xi)) return(Inf)
  denom <- 1 + xi * psi
  if (any(denom <= 0)) return(Inf)
  2 * sum(log(denom))
}

#' Hybrid for ARFIMA(1,d,0)
hybrid_obj_arfima_phi <- function(phi, alpha, Iper, omega, d) {
  w <- obj_Whittle_arfima_phi(phi, Iper, omega, d)
  e <- obj_FDEL_arfima_phi(phi, Iper, omega, d)
  if (!is.finite(e)) return(Inf)
  alpha * w + (1 - alpha) * e
}

# ---- Single-dataset estimators ------------------------------------------------

#' Estimate Whittle and FDEL for AR(1) from one dataset
estimate_whittle_fdel <- function(Iper, omega, start = 0.1) {
  wh <- tryCatch(
    optim(par = start, fn = obj_Whittle, Iper = Iper, omega = omega,
          method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
    error = function(e) NA_real_
  )
  fd <- tryCatch(
    optim(par = start, fn = obj_FDEL,    Iper = Iper, omega = omega,
          method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
    error = function(e) NA_real_
  )
  c(whittle = wh, fdel = fd)
}

#' Estimate the Hybrid for AR(1) at a given alpha
estimate_hybrid_once <- function(alpha, Iper, omega, start = 0.1) {
  res <- tryCatch(
    optim(par = start,
          fn  = function(b) hybrid_obj(b, alpha = alpha, Iper = Iper, omega = omega),
          method = "L-BFGS-B", lower = -0.9, upper = 0.9),
    error = function(e) list(par = NA)
  )
  as.numeric(res$par)
}

#' Estimate Whittle, FDEL and all Hybrid(alpha) for ARMA(1,1)
estimate_whittle_arma <- function(Iper, omega, start = c(0.1, 0.1)) {
  res <- optim(start, obj_Whittle_arma, Iper = Iper, omega = omega,
               method = "L-BFGS-B", lower = c(-0.95, -0.95), upper = c(0.95, 0.95))
  res$par
}

estimate_fdel_arma <- function(Iper, omega, start = c(0.1, 0.1)) {
  res <- optim(start, obj_FDEL_arma, Iper = Iper, omega = omega,
               method = "L-BFGS-B", lower = c(-0.95, -0.95), upper = c(0.95, 0.95))
  res$par
}

estimate_hybrid_arma <- function(alpha, Iper, omega, start) {
  res <- optim(start,
               fn     = function(b) hybrid_obj_arma(b, alpha = alpha, Iper = Iper, omega = omega),
               method = "L-BFGS-B",
               lower  = c(-0.95, -0.95), upper = c(0.95, 0.95))
  res$par
}

#' Estimate Whittle, FDEL and all Hybrid(alpha) for ARFIMA(1,d,0)
estimate_arfima_all <- function(Iper, omega, d, alpha_grid,
                                start = 0.1, lower = -0.95, upper = 0.95) {
  wh <- tryCatch(
    optim(start, obj_Whittle_arfima_phi, Iper = Iper, omega = omega, d = d,
          method = "L-BFGS-B", lower = lower, upper = upper)$par,
    error = function(e) NA_real_
  )
  fd <- tryCatch(
    optim(start, obj_FDEL_arfima_phi,    Iper = Iper, omega = omega, d = d,
          method = "L-BFGS-B", lower = lower, upper = upper)$par,
    error = function(e) NA_real_
  )
  K   <- length(alpha_grid)
  hyb <- numeric(K)
  for (k in seq_along(alpha_grid)) {
    a <- alpha_grid[k]
    if (a == 0) {
      hyb[k] <- fd
    } else if (a == 1) {
      hyb[k] <- wh
    } else {
      hyb[k] <- tryCatch(
        optim(wh,
              fn     = hybrid_obj_arfima_phi,
              alpha  = a, Iper = Iper, omega = omega, d = d,
              method = "L-BFGS-B", lower = lower, upper = upper)$par,
        error = function(e) NA_real_
      )
    }
  }
  list(whittle = wh, fdel = fd, hybrid = hyb)
}

# ---- Spectral distance (misspecification section) -----------------------------

#' Spectral distance between true AR(2) spectrum and fitted AR(1) spectrum
#' Uses trapezoidal integration over [0, pi]
spectral_distance <- function(beta_hat, phi1_true = 0.6, phi2_true = -0.3,
                               grid = seq(0, pi, length.out = 2000)) {
  f_true <- spectral(eta = c(phi1_true, phi2_true), p = 2, q = 0, x = grid)
  f_fit  <- spectral(eta = c(beta_hat),              p = 1, q = 0, x = grid)
  pracma::trapz(grid, abs(f_true - f_fit))
}

# ---- RMSE helper (defined once) -----------------------------------------------

rmse <- function(est, true) sqrt(mean((est - true)^2, na.rm = TRUE))
