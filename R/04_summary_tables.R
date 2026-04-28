## =============================================================================
## 04_summary_tables.R
## Post-processing helpers: summary statistics and LaTeX table generation.
##
## Three layers per model:
##   1. summary_*()        — Bias / SD / RMSE + oracle alpha* for one (innov, n)
##   2. make_big_table_*() — Assemble the long-format data frame
##   3. to_latex_table_*() — Print LaTeX tabular code
##
## Depends on: 01_core_estimators.R  (for rmse())
## =============================================================================

# ---- Internal helper: pick oracle alpha from interior of grid ----------------

.best_alpha <- function(rmse_vec, alpha_grid) {
  interior   <- alpha_grid > 0 & alpha_grid < 1
  alpha_grid[interior][which.min(rmse_vec[interior])]
}

# ---- Internal helper: build 3-stat block (Bias/SD/RMSE) for one estimator ---

.bsr <- function(est, true) {
  c(bias = mean(est - true, na.rm = TRUE),
    sd   = sd(est, na.rm = TRUE),
    rmse = rmse(est, true))
}

# =============================================================================
# AR(1) / AR(1)–ARCH(1) / AR(1)+AO / AR(1) skewt
# (univariate phi, same summary structure)
# =============================================================================

#' Summarise one (innov, n) case for any univariate-phi model
#' @param res      output of run_hybrid_grid_ar1() or similar
#' @param phi_true true parameter value (optional; overrides res$phi_true)
summary_case <- function(res, phi_true = NULL) {
  if (is.null(phi_true)) phi_true <- res$phi_true
  whittle    <- res$whittle
  fdel       <- res$fdel
  hybrid_mat <- res$hybrid
  alpha_grid <- res$alpha_grid

  rmse_alpha <- apply(hybrid_mat, 2, rmse, true = phi_true)
  a_best     <- .best_alpha(rmse_alpha, alpha_grid)
  hyb_best   <- hybrid_mat[, alpha_grid == a_best]

  s_wh  <- .bsr(whittle,  phi_true)
  s_fd  <- .bsr(fdel,     phi_true)
  s_hy  <- .bsr(hyb_best, phi_true)

  c(Whittle_bias  = s_wh["bias"],  Whittle_sd    = s_wh["sd"],
    Whittle_RMSE  = s_wh["rmse"],
    FDEL_bias     = s_fd["bias"],  FDEL_sd       = s_fd["sd"],
    FDEL_RMSE     = s_fd["rmse"],
    Hybrid_bias   = s_hy["bias"],  Hybrid_sd     = s_hy["sd"],
    Hybrid_RMSE   = s_hy["rmse"],
    alpha_best    = a_best)
}

# summary_arch_case() is identical to summary_case() because the driver now
# stores phi_true inside the result list — kept as an alias for clarity.
summary_arch_case <- summary_case

#' Build the long-format table for AR(1) / ARCH(1) / AO / skewt
#' @param results_all nested list: results_all[[innov]][[as.character(n)]]
#' @param phi_true    true parameter (passed to summary_case)
make_big_table_ar1 <- function(results_all, phi_true = 0.4) {
  rows <- list()
  for (innov in names(results_all)) {
    for (n in as.numeric(names(results_all[[innov]]))) {
      res  <- results_all[[innov]][[as.character(n)]]
      summ <- summary_case(res, phi_true)
      rows <- c(rows, list(
        data.frame(Innov = innov, T = n, Stat = "Bias",
                   Whittle = summ["Whittle_bias"], FDEL = summ["FDEL_bias"],
                   Hybrid = summ["Hybrid_bias"],   alpha = NA_real_),
        data.frame(Innov = innov, T = n, Stat = "SD",
                   Whittle = summ["Whittle_sd"],   FDEL = summ["FDEL_sd"],
                   Hybrid = summ["Hybrid_sd"],     alpha = NA_real_),
        data.frame(Innov = innov, T = n, Stat = "RMSE",
                   Whittle = summ["Whittle_RMSE"], FDEL = summ["FDEL_RMSE"],
                   Hybrid = summ["Hybrid_RMSE"],   alpha = summ["alpha_best"])
      ))
    }
  }
  `rownames<-`(do.call(rbind, rows), NULL)
}

# Convenience aliases — same structure, separate names for caller clarity
make_big_table_arch1  <- make_big_table_ar1
make_big_table_ar1_ao <- make_big_table_ar1

#' Print LaTeX table for AR(1)-family models (Bias / SD / RMSE + alpha*)
to_latex_table_ar1 <- function(df, caption = "Table") {
  cat("\\begin{table}[!ht]\n\\centering\n")
  cat(sprintf("\\caption{%s}\n", caption))
  cat("\\begin{tabular}{llcccc}\n\\toprule\n")
  cat("Innovation & $T$ & Statistic & Whittle & FDEL & Hybrid \\\\\n\\midrule\n")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    if (r$Stat == "RMSE") {
      cat(sprintf("%s & %d & %s & %.5f & %.5f & %.5f \\\\[-0.2em]\n",
                  r$Innov, r$T, r$Stat, r$Whittle, r$FDEL, r$Hybrid))
      cat(sprintf("\\multicolumn{6}{r}{$\\alpha^* = %.2f$} \\\\ \n", r$alpha))
      cat("\\addlinespace[0.3em]\n")
    } else {
      cat(sprintf("%s & %d & %s & %.5f & %.5f & %.5f \\\\\n",
                  r$Innov, r$T, r$Stat, r$Whittle, r$FDEL, r$Hybrid))
    }
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
}

# Reuse the same printer for ARCH(1) and AO sections
to_latex_table_arch1  <- to_latex_table_ar1
to_latex_table_ar1_ao <- to_latex_table_ar1

# =============================================================================
# ARMA(1,1)
# =============================================================================

#' Summarise one (innov, n) case for ARMA(1,1)
#' @param res  output of run_hybrid_grid_arma11()
summ_arma_case <- function(res) {
  phi_true   <- res$phi_true
  theta_true <- res$theta_true
  alpha_grid <- res$alpha_grid
  phi_hyb    <- res$phi_hyb
  theta_hyb  <- res$theta_hyb

  rmse_sum  <- apply(phi_hyb, 2, rmse, true = phi_true) +
               apply(theta_hyb, 2, rmse, true = theta_true)
  a_best    <- .best_alpha(rmse_sum, alpha_grid)
  idx_best  <- which(alpha_grid == a_best)

  s_phi_wh  <- .bsr(res$phi_wh,         phi_true)
  s_phi_fd  <- .bsr(res$phi_fd,         phi_true)
  s_phi_hy  <- .bsr(phi_hyb[, idx_best], phi_true)
  s_th_wh   <- .bsr(res$theta_wh,       theta_true)
  s_th_fd   <- .bsr(res$theta_fd,       theta_true)
  s_th_hy   <- .bsr(theta_hyb[, idx_best], theta_true)

  c(phi_Whittle_bias  = s_phi_wh["bias"], phi_Whittle_sd   = s_phi_wh["sd"],
    phi_Whittle_RMSE  = s_phi_wh["rmse"],
    phi_FDEL_bias     = s_phi_fd["bias"], phi_FDEL_sd      = s_phi_fd["sd"],
    phi_FDEL_RMSE     = s_phi_fd["rmse"],
    phi_Hybrid_bias   = s_phi_hy["bias"], phi_Hybrid_sd    = s_phi_hy["sd"],
    phi_Hybrid_RMSE   = s_phi_hy["rmse"],
    theta_Whittle_bias = s_th_wh["bias"], theta_Whittle_sd  = s_th_wh["sd"],
    theta_Whittle_RMSE = s_th_wh["rmse"],
    theta_FDEL_bias    = s_th_fd["bias"], theta_FDEL_sd     = s_th_fd["sd"],
    theta_FDEL_RMSE    = s_th_fd["rmse"],
    theta_Hybrid_bias  = s_th_hy["bias"], theta_Hybrid_sd   = s_th_hy["sd"],
    theta_Hybrid_RMSE  = s_th_hy["rmse"],
    sum_Whittle_RMSE   = s_phi_wh["rmse"] + s_th_wh["rmse"],
    sum_FDEL_RMSE      = s_phi_fd["rmse"] + s_th_fd["rmse"],
    sum_Hybrid_RMSE    = s_phi_hy["rmse"] + s_th_hy["rmse"],
    alpha_best         = a_best)
}

#' Build the long-format table for ARMA(1,1)
make_big_table_arma11 <- function(results_arma11) {
  rows <- list()
  for (innov in names(results_arma11)) {
    for (n in as.numeric(names(results_arma11[[innov]]))) {
      s <- summ_arma_case(results_arma11[[innov]][[as.character(n)]])
      rows <- c(rows, list(
        data.frame(Innovation = innov, T = n, Param = "phi", Stat = "Bias",
                   Whittle = s["phi_Whittle_bias"], FDEL = s["phi_FDEL_bias"],
                   Hybrid  = s["phi_Hybrid_bias"],  alpha = NA_real_),
        data.frame(Innovation = innov, T = n, Param = "phi", Stat = "SD",
                   Whittle = s["phi_Whittle_sd"],   FDEL = s["phi_FDEL_sd"],
                   Hybrid  = s["phi_Hybrid_sd"],    alpha = NA_real_),
        data.frame(Innovation = innov, T = n, Param = "phi", Stat = "RMSE",
                   Whittle = s["phi_Whittle_RMSE"], FDEL = s["phi_FDEL_RMSE"],
                   Hybrid  = s["phi_Hybrid_RMSE"],  alpha = NA_real_),
        data.frame(Innovation = innov, T = n, Param = "theta", Stat = "Bias",
                   Whittle = s["theta_Whittle_bias"], FDEL = s["theta_FDEL_bias"],
                   Hybrid  = s["theta_Hybrid_bias"],  alpha = NA_real_),
        data.frame(Innovation = innov, T = n, Param = "theta", Stat = "SD",
                   Whittle = s["theta_Whittle_sd"],   FDEL = s["theta_FDEL_sd"],
                   Hybrid  = s["theta_Hybrid_sd"],    alpha = NA_real_),
        data.frame(Innovation = innov, T = n, Param = "theta", Stat = "RMSE",
                   Whittle = s["theta_Whittle_RMSE"], FDEL = s["theta_FDEL_RMSE"],
                   Hybrid  = s["theta_Hybrid_RMSE"],  alpha = NA_real_),
        data.frame(Innovation = innov, T = n, Param = "phi+theta", Stat = "RMSE_sum",
                   Whittle = s["sum_Whittle_RMSE"],   FDEL = s["sum_FDEL_RMSE"],
                   Hybrid  = s["sum_Hybrid_RMSE"],    alpha = s["alpha_best"])
      ))
    }
  }
  `rownames<-`(do.call(rbind, rows), NULL)
}

#' Print LaTeX table for ARMA(1,1)
to_latex_table_arma <- function(df, caption = "Table") {
  cat("\\begin{table}[!ht]\n\\centering\n")
  cat(sprintf("\\caption{%s}\n", caption))
  cat("\\begin{tabular}{llllccc}\n\\toprule\n")
  cat("Innovation & $T$ & Parameter & Statistic & Whittle & FDEL & Hybrid \\\\\n\\midrule\n")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    param_tex <- switch(as.character(r$Param),
      "phi"       = "$\\phi$",
      "theta"     = "$\\theta$",
      "phi+theta" = "$\\phi + \\theta$",
      as.character(r$Param))
    stat_tex <- if (r$Stat == "RMSE_sum") "$\\text{RMSE}_{\\phi+\\theta}$" else
                  as.character(r$Stat)
    if (r$Stat == "RMSE_sum") {
      cat(sprintf("%s & %d & %s & %s & %.5f & %.5f & %.5f \\\\[-0.2em]\n",
                  r$Innovation, r$T, param_tex, stat_tex,
                  r$Whittle, r$FDEL, r$Hybrid))
      cat(sprintf("\\multicolumn{7}{r}{$\\alpha^* = %.2f$} \\\\ \n", r$alpha))
      cat("\\addlinespace[0.3em]\n")
    } else {
      cat(sprintf("%s & %d & %s & %s & %.5f & %.5f & %.5f \\\\\n",
                  r$Innovation, r$T, param_tex, stat_tex,
                  r$Whittle, r$FDEL, r$Hybrid))
    }
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
}

# =============================================================================
# ARFIMA(1,d,0)
# =============================================================================

#' Summarise one (d, innov, n) case for ARFIMA(1,d,0)
#' @param res  output of run_hybrid_grid_arfima_phi()
summ_arfima_case <- function(res) {
  phi0       <- res$phi_true
  alpha_grid <- res$alpha_grid
  hybrid_mat <- res$hybrid

  rmse_alpha <- apply(hybrid_mat, 2, rmse, true = phi0)
  a_best     <- .best_alpha(rmse_alpha, alpha_grid)
  hyb_best   <- hybrid_mat[, alpha_grid == a_best]

  s_wh <- .bsr(res$whittle, phi0)
  s_fd <- .bsr(res$fdel,    phi0)
  s_hy <- .bsr(hyb_best,    phi0)

  c(Whittle_bias = s_wh["bias"], Whittle_sd   = s_wh["sd"],
    Whittle_RMSE = s_wh["rmse"],
    FDEL_bias    = s_fd["bias"], FDEL_sd      = s_fd["sd"],
    FDEL_RMSE    = s_fd["rmse"],
    Hybrid_bias  = s_hy["bias"], Hybrid_sd    = s_hy["sd"],
    Hybrid_RMSE  = s_hy["rmse"],
    alpha_best   = a_best)
}

#' Build the long-format table for ARFIMA(1,d,0)
#' @param results_arfima1d0 nested list: [[d_name]][[innov]][[as.character(n)]]
make_big_table_arfima1d0 <- function(results_arfima1d0) {
  rows <- list()
  for (d_name in names(results_arfima1d0)) {
    d_true <- as.numeric(sub("^d", "", d_name))
    for (innov in names(results_arfima1d0[[d_name]])) {
      for (n in as.numeric(names(results_arfima1d0[[d_name]][[innov]]))) {
        s <- summ_arfima_case(
              results_arfima1d0[[d_name]][[innov]][[as.character(n)]])
        rows <- c(rows, list(
          data.frame(d = d_true, Innovation = innov, T = n, Stat = "Bias",
                     Whittle = s["Whittle_bias"], FDEL = s["FDEL_bias"],
                     Hybrid  = s["Hybrid_bias"],  alpha = NA_real_),
          data.frame(d = d_true, Innovation = innov, T = n, Stat = "SD",
                     Whittle = s["Whittle_sd"],   FDEL = s["FDEL_sd"],
                     Hybrid  = s["Hybrid_sd"],    alpha = NA_real_),
          data.frame(d = d_true, Innovation = innov, T = n, Stat = "RMSE",
                     Whittle = s["Whittle_RMSE"], FDEL = s["FDEL_RMSE"],
                     Hybrid  = s["Hybrid_RMSE"],  alpha = s["alpha_best"])
        ))
      }
    }
  }
  `rownames<-`(do.call(rbind, rows), NULL)
}

#' Print LaTeX table for ARFIMA(1,d,0)
to_latex_table_arfima <- function(df, caption = "Table") {
  cat("\\begin{table}[!ht]\n\\centering\n")
  cat(sprintf("\\caption{%s}\n", caption))
  cat("\\begin{tabular}{lllcccc}\n\\toprule\n")
  cat("$d$ & Innovation & $T$ & Statistic & Whittle & FDEL & Hybrid \\\\\n\\midrule\n")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    if (r$Stat == "RMSE") {
      cat(sprintf("%.1f & %s & %d & %s & %.5f & %.5f & %.5f \\\\[-0.2em]\n",
                  r$d, r$Innovation, r$T, r$Stat, r$Whittle, r$FDEL, r$Hybrid))
      cat(sprintf("\\multicolumn{7}{r}{$\\alpha^* = %.2f$} \\\\ \n", r$alpha))
      cat("\\addlinespace[0.3em]\n")
    } else {
      cat(sprintf("%.1f & %s & %d & %s & %.5f & %.5f & %.5f \\\\\n",
                  r$d, r$Innovation, r$T, r$Stat, r$Whittle, r$FDEL, r$Hybrid))
    }
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
}

# =============================================================================
# AR(2) misspecification
# =============================================================================

#' Summarise one (innov, n) case for the AR(2) misspecification study
#' @param res  output of run_hybrid_grid_ar2_misspec()
summary_case_spectral <- function(res) {
  alpha_grid   <- res$alpha_grid
  dist_hybrid  <- res$dist_hybrid

  mean_dist_alpha <- colMeans(dist_hybrid, na.rm = TRUE)
  a_best          <- .best_alpha(mean_dist_alpha, alpha_grid)
  dist_best       <- dist_hybrid[, alpha_grid == a_best]

  c(Whittle   = mean(res$dist_whittle, na.rm = TRUE),
    FDEL      = mean(res$dist_fdel,    na.rm = TRUE),
    Hybrid    = mean(dist_best,        na.rm = TRUE),
    alpha_best = a_best)
}

#' Build the long-format table for the AR(2) misspecification study
make_big_table_ar2_misspec <- function(results_all_misspec) {
  rows <- list()
  for (innov in names(results_all_misspec)) {
    for (n in as.numeric(names(results_all_misspec[[innov]]))) {
      s <- summary_case_spectral(results_all_misspec[[innov]][[as.character(n)]])
      rows <- c(rows, list(
        data.frame(Innov = innov, T = n,
                   Whittle = s["Whittle"], FDEL = s["FDEL"],
                   Hybrid  = s["Hybrid"],  alpha = s["alpha_best"],
                   row.names = NULL)
      ))
    }
  }
  `rownames<-`(do.call(rbind, rows), NULL)
}

#' Print LaTeX table for the AR(2) misspecification study
to_latex_table_ar2_misspec <- function(df, caption = "Table") {
  cat("\\begin{table}[!ht]\n\\centering\n")
  cat(sprintf("\\caption{%s}\n", caption))
  cat("\\begin{tabular}{lccccc}\n\\toprule\n")
  cat("Innovation & $T$ & Whittle & FDEL & Hybrid & \\\\\n\\midrule\n")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    cat(sprintf("%s & %d & %.5f & %.5f & %.5f \\\\[-0.2em]\n",
                r$Innov, r$T, r$Whittle, r$FDEL, r$Hybrid))
    cat(sprintf("\\multicolumn{6}{r}{$\\alpha^* = %.2f$} \\\\ \n", r$alpha))
    cat("\\addlinespace[0.3em]\n")
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
}
