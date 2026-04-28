# Hybrid Likelihood Methods in the Frequency Domain for Time Series Estimation

**Master Thesis — Akshaan Murugesu (21-327-671)**  
University of Geneva — Supervisors: Prof. Davide La Vecchia and Manon Felix

---

## Overview

This repository contains the R code for the Monte Carlo simulation study in the
thesis. The study evaluates a **Hybrid Likelihood estimator** that interpolates
between the **Whittle parametric likelihood** and the **Frequency-Domain
Empirical Likelihood (FDEL)** through a tuning parameter α ∈ [0, 1]:

ℓ_HL(θ; α) = α · ℓ_WL(θ) + (1 − α) · ℓ_EL(θ)

The simulations cover seven data-generating mechanisms:
1. AR(1) — Gaussian and non-Gaussian innovations
2. ARMA(1,1) — multiple innovation distributions
3. ARFIMA(1, d, 0) — long-memory processes (d = 0.1 and d = 0.4)
4. AR(1)–ARCH(1) — conditional heteroskedasticity
5. AR(1) with additive outliers — 5% contamination
6. AR(1) with skewed heavy-tailed innovations — split-t distribution
7. AR(2) DGM, AR(1) fit — model misspecification

---

## Repository structure

```
.
├── R/
│   ├── 01_core_estimators.R   # Periodogram, spectral density, Whittle/FDEL/Hybrid
│   ├── 02_dgm.R               # All data-generating mechanisms
│   ├── 03_drivers.R           # Monte Carlo simulation drivers (one per model)
│   └── 04_summary_tables.R    # Summary statistics + LaTeX table generation
│
├── simulations/
│   ├── 01_sim_ar1.Rmd         # AR(1) study
│   ├── 02_sim_arma11.Rmd      # ARMA(1,1) study
│   ├── 03_sim_arfima.Rmd      # ARFIMA(1,d,0) study
│   ├── 04_sim_arch1.Rmd       # AR(1)–ARCH(1) study
│   ├── 05_sim_ar1_ao.Rmd      # AR(1) + additive outliers
│   ├── 06_sim_ar1_skewt.Rmd   # AR(1) with split-t innovations
│   └── 07_sim_ar2_misspec.Rmd # AR(2) DGM → AR(1) fit (misspecification)
│
└── outputs/                   # .rds result files (gitignored, generated on first run)
```

---

## Reproducibility

All simulations use `set.seed(2025)` and 1 000 Monte Carlo replications.
Results are cached as `.rds` files in `outputs/`. Each Rmd checks whether
the cache exists and skips the (slow) simulation if it does:

```r
if (file.exists(rds_file)) {
  results <- readRDS(rds_file)
} else {
  # run simulation and saveRDS(results, rds_file)
}
```

To reproduce from scratch, delete the relevant `.rds` file(s) in `outputs/`
and re-knit the corresponding Rmd.

> **Runtime note.** Full reproduction of all simulations can take several hours
> on a single core. All drivers support `ncores > 1` via `parallel::mclapply`
> (macOS / Linux). Set `ncores <- 1` on Windows.

---

## Dependencies

```r
install.packages(c(
  "parallel",   # parallel computation (base package on most systems)
  "pracma",     # trapezoidal integration (spectral_distance)
  "here",       # project-relative paths
  "knitr",      # knitting Rmd files
  "rmarkdown"   # knitting Rmd files
))
```

The `here` package requires that you work from the project root (where this
`README.md` is located). Open the `.Rproj` file if using RStudio.

---

## Running the simulations

From R (or RStudio), knit each file in order:

```r
library(rmarkdown)

render("simulations/01_sim_ar1.Rmd")
render("simulations/02_sim_arma11.Rmd")
render("simulations/03_sim_arfima.Rmd")
render("simulations/04_sim_arch1.Rmd")
render("simulations/05_sim_ar1_ao.Rmd")
render("simulations/06_sim_ar1_skewt.Rmd")
render("simulations/07_sim_ar2_misspec.Rmd")
```

Or source only the R scripts to re-run simulations without producing HTML:

```r
source("R/01_core_estimators.R")
source("R/02_dgm.R")
source("R/03_drivers.R")
source("R/04_summary_tables.R")
```

---

## Key design parameters (thesis values)

| Model           | phi_true | theta_true | d    | n_values               | nsim |
|-----------------|----------|------------|------|------------------------|------|
| AR(1)           | 0.4      | —          | —    | 50, 75, 100, 200, 300  | 1000 |
| ARMA(1,1)       | 0.2      | 0.4        | —    | 50, 75, 100, 200, 300  | 1000 |
| ARFIMA(1,d,0)   | 0.4      | —          | 0.1, 0.4 | 50, 75, 100, 200, 300 | 1000 |
| AR(1)–ARCH(1)   | 0.4      | —          | —    | 50, 75, 100, 200, 300  | 1000 |
| AR(1) + AO      | 0.4      | —          | —    | 50, 75, 100, 200        | 1000 |
| AR(1) split-t   | 0.4      | —          | —    | 50, 75, 100, 200, 300  | 1000 |
| AR(2) → AR(1)   | (0.6, −0.3) | —       | —    | 50, 75, 100, 200, 300  | 1000 |

ARCH(1) parameters: (α₀, α₁) = (0.5, 0.4).  
Additive outliers: 5% contamination, shock sd = 5.  
Split-t: ν = 5, a = 1.6, b = 0.8.

---

## Code organisation notes

The original single-file `Thesis_functions.R` had several redundancies that
have been resolved in this version:

- `rmse()` was defined four times — now defined once in `01_core_estimators.R`
- `estimate_hybrid_once()` was defined twice — now defined once
- `solve_xi()` and `solve_xi_scalar()` were identical — unified as `solve_xi()`
- `simulate_ar2()` and `simulate_ar1_ao()` (old versions) duplicated newer
  generic variants — removed, only the generic `*_innov()` variants remain
- The `one_rep` / `run_hybrid_grid_*` pattern was copy-pasted six times —
  factored into shared internal helpers `.run_reps()`, `.collect_results()`,
  and `.hybrid_over_grid()` in `03_drivers.R`

---

## Citation

Murugesu, A. (2025). *Hybrid Likelihood Methods in the Frequency Domain for
Time Series Estimation*. Master Thesis, University of Geneva.
