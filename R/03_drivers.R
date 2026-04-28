## =============================================================================
## 03_drivers.R
## Monte Carlo simulation drivers for each model class.
##
## All drivers share the same output format:
##   list(n, phi_true, alpha_grid, whittle, fdel, hybrid [matrix])
##
## Depends on: 01_core_estimators.R, 02_dgm.R
## =============================================================================

# ---- Internal helper: run one_rep in parallel or sequentially ----------------

.run_reps <- function(one_rep, nsim, ncores = 1) {
  if (ncores > 1) {
    parallel::mclapply(seq_len(nsim), one_rep, mc.cores = ncores)
  } else {
    lapply(seq_len(nsim), one_rep)
  }
}

# ---- Internal helper: collect scalar whittle/fdel + hybrid matrix ------------

.collect_results <- function(res_list, alpha_grid) {
  K       <- length(alpha_grid)
  whittle <- vapply(res_list, function(z) z$wh, numeric(1))
  fdel    <- vapply(res_list, function(z) z$fd, numeric(1))
  hybrid  <- t(vapply(res_list, function(z) z$hy, numeric(K)))
  colnames(hybrid) <- paste0("a", alpha_grid)
  list(whittle = whittle, fdel = fdel, hybrid = hybrid)
}

# ---- Internal helper: AR(1)-family hybrid estimation over alpha-grid ---------
# Used by all univariate-parameter drivers (AR1, ARCH1, AO, skewt, AR2)

.hybrid_over_grid <- function(wh, fd, alpha_grid, Iper, omega) {
  K  <- length(alpha_grid)
  hy <- numeric(K)
  for (k in seq_along(alpha_grid)) {
    a <- alpha_grid[k]
    if (a == 0) {
      hy[k] <- fd
    } else if (a == 1) {
      hy[k] <- wh
    } else {
      hy[k] <- estimate_hybrid_once(
        alpha = a, Iper = Iper, omega = omega, start = wh
      )
    }
  }
  hy
}

# ---- AR(1) driver ------------------------------------------------------------

#' Monte Carlo simulation for AR(1) with generic innovations
#'
#' @param n          sample size
#' @param nsim       number of replications
#' @param phi_true   true AR(1) coefficient
#' @param alpha_grid numeric vector of alpha values in [0,1]
#' @param innov      innovation distribution (see r_innov)
#' @param start      starting value for optimisers
#' @param ncores     number of cores (>1 uses mclapply)
#' @return list with elements: n, phi_true, alpha_grid, whittle, fdel, hybrid
run_hybrid_grid_ar1 <- function(n, nsim, phi_true,
                                alpha_grid,
                                innov  = c("gaussian", "t5", "chisq5", "exp1"),
                                start  = 0.1,
                                ncores = 1) {
  innov <- match.arg(innov)

  one_rep <- function(s) {
    Y     <- simulate_ar1_innov(n = n, phi = phi_true, innov = innov)
    Iper  <- per(Y)
    omega <- freq(Y)
    wh <- tryCatch(
      optim(start, obj_Whittle, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    fd <- tryCatch(
      optim(start, obj_FDEL, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    list(wh = wh, fd = fd, hy = .hybrid_over_grid(wh, fd, alpha_grid, Iper, omega))
  }

  res <- .collect_results(.run_reps(one_rep, nsim, ncores), alpha_grid)
  c(list(n = n, phi_true = phi_true, alpha_grid = alpha_grid), res)
}

# ---- ARMA(1,1) driver --------------------------------------------------------

#' Monte Carlo simulation for ARMA(1,1) with generic innovations
#'
#' @param n           sample size
#' @param nsim        number of replications
#' @param phi_true    true AR coefficient
#' @param theta_true  true MA coefficient
#' @param alpha_grid  numeric vector of alpha values
#' @param innov       innovation distribution
#' @param ncores      number of cores
#' @return list with elements: n, phi_true, theta_true, alpha_grid,
#'         phi_wh, theta_wh, phi_fd, theta_fd, phi_hyb, theta_hyb
run_hybrid_grid_arma11 <- function(n, nsim,
                                   phi_true, theta_true,
                                   alpha_grid,
                                   innov  = c("gaussian", "t5", "chisq5", "exp1"),
                                   ncores = 1) {
  innov <- match.arg(innov)
  K     <- length(alpha_grid)

  one_rep <- function(s) {
    Y     <- simulate_arma11_innov(n = n, phi = phi_true,
                                   theta = theta_true, innov = innov)
    Iper  <- per(Y)
    omega <- freq(Y)
    wh    <- estimate_whittle_arma(Iper, omega, start = c(phi_true, theta_true))
    fd    <- estimate_fdel_arma(Iper, omega, start = wh)

    phi_hy <- numeric(K)
    th_hy  <- numeric(K)
    for (k in seq_along(alpha_grid)) {
      a   <- alpha_grid[k]
      est <- if (a == 0) fd else if (a == 1) wh else
        estimate_hybrid_arma(a, Iper, omega, start = wh)
      phi_hy[k] <- est[1]
      th_hy[k]  <- est[2]
    }
    list(phi_wh = wh[1], theta_wh = wh[2],
         phi_fd = fd[1], theta_fd = fd[2],
         phi_hy = phi_hy, th_hy = th_hy)
  }

  res_list <- .run_reps(one_rep, nsim, ncores)

  phi_hyb   <- t(vapply(res_list, function(z) z$phi_hy, numeric(K)))
  theta_hyb <- t(vapply(res_list, function(z) z$th_hy,  numeric(K)))
  colnames(phi_hyb) <- colnames(theta_hyb) <- paste0("a", alpha_grid)

  list(
    n          = n,
    phi_true   = phi_true,
    theta_true = theta_true,
    alpha_grid = alpha_grid,
    phi_wh     = vapply(res_list, function(z) z$phi_wh,   numeric(1)),
    theta_wh   = vapply(res_list, function(z) z$theta_wh, numeric(1)),
    phi_fd     = vapply(res_list, function(z) z$phi_fd,   numeric(1)),
    theta_fd   = vapply(res_list, function(z) z$theta_fd, numeric(1)),
    phi_hyb    = phi_hyb,
    theta_hyb  = theta_hyb
  )
}

# ---- ARFIMA(1,d,0) driver ----------------------------------------------------

#' Monte Carlo simulation for ARFIMA(1,d,0) — phi only, d treated as known
#'
#' @param n          sample size
#' @param nsim       number of replications
#' @param phi_true   true AR coefficient
#' @param d_true     true fractional differencing parameter
#' @param alpha_grid numeric vector of alpha values
#' @param innov      innovation distribution
#' @param ncores     number of cores
#' @return list with elements: n, phi_true, d_true, alpha_grid,
#'         whittle, fdel, hybrid
run_hybrid_grid_arfima_phi <- function(n, nsim,
                                       phi_true, d_true,
                                       alpha_grid,
                                       innov  = c("gaussian", "t5", "chisq5", "exp1"),
                                       ncores = 1) {
  innov <- match.arg(innov)
  K     <- length(alpha_grid)

  one_rep <- function(s) {
    Y     <- simulate_arfima10_innov(n = n, phi = phi_true,
                                     d = d_true, innov = innov)
    Iper  <- per(Y)
    omega <- freq(Y)
    est   <- estimate_arfima_all(Iper = Iper, omega = omega, d = d_true,
                                 alpha_grid = alpha_grid, start = phi_true)
    list(wh = est$whittle, fd = est$fdel, hy = est$hybrid)
  }

  res <- .collect_results(.run_reps(one_rep, nsim, ncores), alpha_grid)
  c(list(n = n, phi_true = phi_true, d_true = d_true, alpha_grid = alpha_grid), res)
}

# ---- AR(1)–ARCH(1) driver ----------------------------------------------------

#' Monte Carlo simulation for AR(1) with ARCH(1) volatility
#'
#' @param n          sample size
#' @param nsim       number of replications
#' @param phi_true   true AR coefficient
#' @param alpha_grid numeric vector of alpha values
#' @param innov      innovation distribution
#' @param alpha0     ARCH intercept  (default 0.5)
#' @param alpha1     ARCH slope      (default 0.4)
#' @param start      starting value
#' @param ncores     number of cores
run_hybrid_grid_arch <- function(n, nsim,
                                 phi_true,
                                 alpha_grid,
                                 innov  = c("gaussian", "t5", "chisq5", "exp1"),
                                 alpha0 = 0.5, alpha1 = 0.4,
                                 start  = 0.1,
                                 ncores = 1) {
  innov <- match.arg(innov)

  one_rep <- function(s) {
    Y     <- simulate_ar1_arch1(n = n, phi = phi_true,
                                alpha0 = alpha0, alpha1 = alpha1, innov = innov)
    Iper  <- per(Y)
    omega <- freq(Y)
    wh <- tryCatch(
      optim(start, obj_Whittle, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    fd <- tryCatch(
      optim(start, obj_FDEL, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    list(wh = wh, fd = fd, hy = .hybrid_over_grid(wh, fd, alpha_grid, Iper, omega))
  }

  res <- .collect_results(.run_reps(one_rep, nsim, ncores), alpha_grid)
  c(list(n = n, phi_true = phi_true, alpha_grid = alpha_grid), res)
}

# ---- AR(1) + Additive Outliers driver ----------------------------------------

#' Monte Carlo simulation for AR(1) contaminated with additive outliers
#'
#' @param n          sample size
#' @param nsim       number of replications
#' @param phi_true   true AR coefficient
#' @param alpha_grid numeric vector of alpha values
#' @param innov      innovation distribution
#' @param p_out      proportion of outliers (default 0.05)
#' @param ao_sd      outlier standard deviation (default 5)
#' @param start      starting value
#' @param ncores     number of cores
run_hybrid_grid_ar1_ao <- function(n, nsim, phi_true, alpha_grid,
                                   innov  = c("gaussian", "t5", "chisq5", "exp1"),
                                   p_out  = 0.05, ao_sd = 5,
                                   start  = 0.1,
                                   ncores = 1) {
  innov <- match.arg(innov)

  one_rep <- function(s) {
    Y     <- simulate_ar1_ao_innov(n = n, phi = phi_true,
                                   innov = innov, p_out = p_out, ao_sd = ao_sd)
    Iper  <- per(Y)
    omega <- freq(Y)
    wh <- tryCatch(
      optim(start, obj_Whittle, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    fd <- tryCatch(
      optim(start, obj_FDEL, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    list(wh = wh, fd = fd, hy = .hybrid_over_grid(wh, fd, alpha_grid, Iper, omega))
  }

  res <- .collect_results(.run_reps(one_rep, nsim, ncores), alpha_grid)
  c(list(n = n, phi_true = phi_true, alpha_grid = alpha_grid), res)
}

# ---- AR(1) with skewed heavy-tailed (split-t) driver -------------------------

#' Monte Carlo simulation for AR(1) with split-t innovations
#'
#' @param n          sample size
#' @param nsim       number of replications
#' @param phi_true   true AR coefficient
#' @param alpha_grid numeric vector of alpha values
#' @param df         degrees of freedom (default 5)
#' @param a          positive-tail scale (default 1.6)
#' @param b          negative-tail scale (default 0.8)
#' @param start      starting value
#' @param ncores     number of cores
run_hybrid_grid_ar1_skewt <- function(n, nsim, phi_true, alpha_grid,
                                      df     = 5, a = 1.6, b = 0.8,
                                      start  = 0.1,
                                      ncores = 1) {
  one_rep <- function(s) {
    Y     <- simulate_ar1_skewt(n = n, phi = phi_true, df = df, a = a, b = b)
    Iper  <- per(Y)
    omega <- freq(Y)
    wh <- tryCatch(
      optim(start, obj_Whittle, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    fd <- tryCatch(
      optim(start, obj_FDEL, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    list(wh = wh, fd = fd, hy = .hybrid_over_grid(wh, fd, alpha_grid, Iper, omega))
  }

  res <- .collect_results(.run_reps(one_rep, nsim, ncores), alpha_grid)
  c(list(n = n, phi_true = phi_true, alpha_grid = alpha_grid), res)
}

# ---- AR(2) misspecification driver -------------------------------------------

#' Monte Carlo simulation: true AR(2), fitted AR(1) — evaluates spectral distance
#'
#' @param n          sample size
#' @param nsim       number of replications
#' @param phi1_true  true AR(2) first coefficient  (default 0.6)
#' @param phi2_true  true AR(2) second coefficient (default -0.3)
#' @param alpha_grid numeric vector of alpha values
#' @param innov      innovation distribution
#' @param start      starting value
#' @param ncores     number of cores
#' @return list with elements: n, phi1_true, phi2_true, alpha_grid,
#'         dist_whittle, dist_fdel, dist_hybrid (matrix nsim x K)
run_hybrid_grid_ar2_misspec <- function(n, nsim,
                                        phi1_true = 0.6, phi2_true = -0.3,
                                        alpha_grid,
                                        innov  = c("gaussian", "t5", "chisq5", "exp1"),
                                        start  = 0.1,
                                        ncores = 1) {
  innov <- match.arg(innov)
  K     <- length(alpha_grid)

  one_rep <- function(s) {
    Y     <- simulate_ar2_innov(n = n, phi1 = phi1_true, phi2 = phi2_true,
                                innov = innov)
    Iper  <- per(Y)
    omega <- freq(Y)
    wh <- tryCatch(
      optim(start, obj_Whittle, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    fd <- tryCatch(
      optim(start, obj_FDEL, Iper = Iper, omega = omega,
            method = "L-BFGS-B", lower = -0.9, upper = 0.9)$par,
      error = function(e) NA_real_
    )
    hy    <- .hybrid_over_grid(wh, fd, alpha_grid, Iper, omega)

    # spectral distances
    d_wh  <- spectral_distance(wh,  phi1_true = phi1_true, phi2_true = phi2_true)
    d_fd  <- spectral_distance(fd,  phi1_true = phi1_true, phi2_true = phi2_true)
    d_hy  <- vapply(hy, spectral_distance, numeric(1),
                    phi1_true = phi1_true, phi2_true = phi2_true)

    list(wh = wh, fd = fd, hy = hy,
         dist_wh = d_wh, dist_fd = d_fd, dist_hy = d_hy)
  }

  res_list     <- .run_reps(one_rep, nsim, ncores)
  dist_whittle <- vapply(res_list, function(z) z$dist_wh, numeric(1))
  dist_fdel    <- vapply(res_list, function(z) z$dist_fd, numeric(1))
  dist_hybrid  <- t(vapply(res_list, function(z) z$dist_hy, numeric(K)))
  colnames(dist_hybrid) <- paste0("a", alpha_grid)

  list(
    n            = n,
    phi1_true    = phi1_true,
    phi2_true    = phi2_true,
    alpha_grid   = alpha_grid,
    dist_whittle = dist_whittle,
    dist_fdel    = dist_fdel,
    dist_hybrid  = dist_hybrid
  )
}
