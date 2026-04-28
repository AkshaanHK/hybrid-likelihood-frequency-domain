## =============================================================================
## 02_dgm.R
## Data-Generating Mechanisms (DGM) for all simulation scenarios:
##   - Innovation sampler (shared helper)
##   - AR(1)  — Gaussian and generic innovations
##   - ARMA(1,1) — generic innovations
##   - ARFIMA(1,d,0) — generic innovations
##   - AR(1)–ARCH(1) — generic innovations
##   - AR(1) + Additive Outliers — generic innovations
##   - AR(1) with skewed heavy-tailed (split-t) innovations
##   - AR(2) — generic innovations  (used for misspecification study)
## =============================================================================

# ---- Shared innovation sampler -----------------------------------------------

#' Draw n standardised innovations (mean 0, var ≈ 1) from the chosen family
#' @param n     sample size
#' @param innov one of "gaussian", "t5", "chisq5", "exp1"
r_innov <- function(n, innov = c("gaussian", "t5", "chisq5", "exp1")) {
  innov <- match.arg(innov)
  switch(
    innov,
    "gaussian" = rnorm(n, 0, 1),
    "t5" = {
      z <- rt(n, df = 5)
      z * sqrt((5 - 2) / 5)          # scale to var = 1
    },
    "chisq5" = {
      z <- rchisq(n, df = 5)
      (z - 5) / sqrt(10)             # center & scale
    },
    "exp1" = {
      z <- rexp(n, rate = 1)
      z - 1                          # center to mean 0
    }
  )
}

# ---- AR(1) -------------------------------------------------------------------

#' Simulate AR(1):  X_t = phi * X_{t-1} + eps_t,   eps_t ~ r_innov()
#' @param n      series length
#' @param phi    AR coefficient
#' @param innov  innovation distribution (see r_innov)
simulate_ar1_innov <- function(n, phi,
                               innov = c("gaussian", "t5", "chisq5", "exp1")) {
  innov <- match.arg(innov)
  eps   <- r_innov(n, innov = innov)
  x     <- numeric(n)
  x[1]  <- eps[1]
  for (t in 2:n) x[t] <- phi * x[t - 1] + eps[t]
  as.numeric(x)
}

# ---- ARMA(1,1) ---------------------------------------------------------------

#' Simulate ARMA(1,1):  X_t = phi * X_{t-1} + eps_t + theta * eps_{t-1}
#' @param n      series length
#' @param phi    AR coefficient
#' @param theta  MA coefficient
#' @param innov  innovation distribution
simulate_arma11_innov <- function(n, phi, theta,
                                  innov = c("gaussian", "t5", "chisq5", "exp1")) {
  innov <- match.arg(innov)
  eps   <- r_innov(n, innov = innov)
  x     <- numeric(n)
  x[1]  <- eps[1]
  for (t in 2:n) x[t] <- phi * x[t - 1] + eps[t] + theta * eps[t - 1]
  as.numeric(x)
}

# ---- ARFIMA(1,d,0) -----------------------------------------------------------

#' Fractional-differencing weights for (1 - B)^(-d)  (used as MA filter)
fracdiff_weights <- function(d, M) {
  w    <- numeric(M + 1)
  w[1] <- 1
  for (k in 1:M) w[k + 1] <- w[k] * ((k - 1 + d) / k)
  w
}

#' Simulate ARFIMA(1,d,0):  (1 - phi*B)(1 - B)^d  X_t = eps_t
#' @param n      series length
#' @param phi    AR coefficient
#' @param d      fractional differencing parameter
#' @param innov  innovation distribution
#' @param M      burn-in length for the fractional filter
simulate_arfima10_innov <- function(n, phi, d,
                                    innov = c("gaussian", "t5", "chisq5", "exp1"),
                                    M = 1000) {
  innov <- match.arg(innov)
  e     <- r_innov(n + M, innov = innov)

  # fractional integration step: (1-B)^{-d} e_t
  w <- fracdiff_weights(d, M)
  y <- stats::filter(e, w, method = "convolution", sides = 1)
  y <- as.numeric(y)[(M + 1):(M + n)]    # drop burn-in

  # AR(1) step:  X_t = phi * X_{t-1} + y_t
  x    <- numeric(n)
  x[1] <- y[1]
  for (t in 2:n) x[t] <- phi * x[t - 1] + y[t]
  as.numeric(x)
}

# ---- AR(1)–ARCH(1) -----------------------------------------------------------

#' Simulate AR(1)–ARCH(1):
#'   X_t = phi * X_{t-1} + sigma_t * z_t
#'   sigma_t^2 = alpha0 + alpha1 * X_{t-1}^2
#' @param n      series length
#' @param phi    AR coefficient
#' @param alpha0 ARCH intercept
#' @param alpha1 ARCH slope
#' @param innov  innovation distribution for z_t
simulate_ar1_arch1 <- function(n, phi,
                                alpha0 = 0.5, alpha1 = 0.4,
                                innov = c("gaussian", "t5", "chisq5", "exp1")) {
  innov   <- match.arg(innov)
  eps     <- r_innov(n, innov = innov)   # standardised z_t
  x       <- numeric(n)
  sig2    <- numeric(n)
  sig2[1] <- alpha0 / (1 - alpha1)      # unconditional variance
  x[1]    <- sqrt(sig2[1]) * eps[1]
  for (t in 2:n) {
    sig2[t] <- alpha0 + alpha1 * x[t - 1]^2
    x[t]    <- phi * x[t - 1] + sqrt(sig2[t]) * eps[t]
  }
  as.numeric(x)
}

# ---- AR(1) + Additive Outliers (AO) -----------------------------------------

#' Simulate AR(1) contaminated with additive outliers:
#'   Y_t = X_t + AO_t,  AO_t = N(0, ao_sd^2) at fraction p_out of time points
#' @param n      series length
#' @param phi    AR coefficient
#' @param innov  innovation distribution for the latent AR(1)
#' @param p_out  proportion of observations contaminated
#' @param ao_sd  standard deviation of the additive shock
simulate_ar1_ao_innov <- function(n, phi,
                                  innov = c("gaussian", "t5", "chisq5", "exp1"),
                                  p_out = 0.05,
                                  ao_sd = 5) {
  innov <- match.arg(innov)
  eps   <- r_innov(n, innov = innov)
  X     <- numeric(n)
  X[1]  <- eps[1]
  for (t in 2:n) X[t] <- phi * X[t - 1] + eps[t]

  # additive outliers
  AO  <- numeric(n)
  k   <- max(1, round(p_out * n))
  idx <- sort(sample.int(n, k))
  AO[idx] <- rnorm(k, mean = 0, sd = ao_sd)

  as.numeric(X + AO)
}

# ---- AR(1) with split-t (skewed heavy-tailed) innovations --------------------

#' Draw n split-t innovations (standardised: mean 0, var 1)
#' Positive tail scaled by a, negative tail by b
#' @param n   sample size
#' @param df  degrees of freedom
#' @param a   scale for positive part
#' @param b   scale for negative part
rsplit_t_std <- function(n, df = 5, a = 1.6, b = 0.8) {
  z <- rt(n, df = df)
  z <- ifelse(z >= 0, a * z, b * z)
  z <- z - mean(z)                     # centre
  z / sd(z)                            # scale to var ≈ 1
}

#' Simulate AR(1) with split-t innovations
#' @param n   series length
#' @param phi AR coefficient
#' @param df  degrees of freedom
#' @param a   scale for positive part
#' @param b   scale for negative part
simulate_ar1_skewt <- function(n, phi, df = 5, a = 1.6, b = 0.8) {
  eps  <- rsplit_t_std(n, df = df, a = a, b = b)
  x    <- numeric(n)
  for (t in 2:n) x[t] <- phi * x[t - 1] + eps[t]
  as.numeric(x)
}

# ---- AR(2) (for misspecification study) ---------------------------------------

#' Simulate AR(2):  X_t = phi1 * X_{t-1} + phi2 * X_{t-2} + eps_t
#' @param n     series length
#' @param phi1  first AR coefficient
#' @param phi2  second AR coefficient
#' @param innov innovation distribution
simulate_ar2_innov <- function(n, phi1 = 0.6, phi2 = -0.3,
                               innov = c("gaussian", "t5", "chisq5", "exp1")) {
  innov <- match.arg(innov)
  e     <- r_innov(n, innov = innov)
  x     <- numeric(n)
  x[1]  <- e[1]
  x[2]  <- e[2]
  for (t in 3:n) x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + e[t]
  as.numeric(x)
}
