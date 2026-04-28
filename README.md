# Hybrid Likelihood for Frequency-Domain Time Series Estimation

## Overview

This project proposes and evaluates a **Hybrid Likelihood estimator** for stationary
time series in the frequency domain. The estimator interpolates between the
**Whittle parametric likelihood** and the **Frequency-Domain Empirical Likelihood
(FDEL)** through a tuning parameter α ∈ [0, 1]:

> ℓ\_HL(θ; α) = α · ℓ\_Whittle(θ) + (1 − α) · ℓ\_FDEL(θ)

- α = 1 → pure Whittle (parametric, efficient under correct specification)
- α = 0 → pure FDEL (nonparametric, robust but unstable in small samples)
- α ∈ (0, 1) → hybrid combining both

The finite-sample behaviour is assessed through an extensive Monte Carlo study
across seven data-generating processes, 1 000 replications, and four innovation
distributions.

📄 [Read the full thesis (PDF)](Murugesu_2025_HybridLikelihood.pdf)

---

## Research Question

Can blending the Whittle and FDEL likelihoods produce an estimator that is both
more robust than Whittle and more stable than FDEL in finite samples?

---

## Main Results

The answer is nuanced. Across all seven models studied, the Hybrid estimator
**collapses to Whittle** as soon as α > 0 — even for very small values like α = 0.05.

Key findings:

- **FDEL is the fragile component.** In small samples (T = 50) under heavy-tailed,
  skewed, or heteroskedastic innovations, FDEL exhibits inflated variance and large
  RMSE. Whittle remains stable throughout.

- **The Hybrid stabilises FDEL, not the reverse.** Even a tiny Whittle weight is
  enough to regularise the estimator. There is no smooth interpolation — the RMSE
  surface drops sharply at α = 0 and stays flat for all α > 0.

- **Hybrid marginally outperforms Whittle in the worst cases.** Small but consistent
  gains appear at T = 50 under ARCH(1) volatility, skewed innovations, or long memory.
  These gains vanish as T grows.

- **Additive outliers affect all three estimators equally.** The dominant effect is
  systematic bias, not variance — none of the three approaches is robust to gross
  contamination.

- **Under model misspecification**, the Hybrid tracks Whittle in spectral distance,
  while FDEL produces substantially larger spectral errors, especially in small samples.

---

## Results Summary (AR(1), φ₀ = 0.4, T = 50)

Representative RMSE across innovation distributions. Lower is better.

| Innovation | Whittle | FDEL   | Hybrid (oracle α) |
|------------|---------|--------|-------------------|
| N(0,1)     | 0.1401  | 0.2894 | 0.1401            |
| t₅         | 0.1418  | 0.3104 | 0.1417            |
| χ²(5)      | 0.1397  | 0.3111 | 0.1397            |
| Exp(1)     | 0.1409  | 0.3201 | 0.1394            |

The Hybrid matches or marginally beats Whittle in all cases. FDEL is 2–3× worse.

---

## Simulation Design

Seven data-generating processes:

- **AR(1)** — φ₀ = 0.4, four innovation distributions
- **ARMA(1,1)** — (φ₀, θ₀) = (0.2, 0.4), four innovation distributions
- **ARFIMA(1, d, 0)** — d ∈ {0.1, 0.4}, long-range dependence
- **AR(1)–ARCH(1)** — conditional heteroskedasticity, (α₀, α₁) = (0.5, 0.4)
- **AR(1) + Additive Outliers** — 5% contamination, spike sd = 5
- **AR(1) with split-t innovations** — ν = 5, a = 1.6, b = 0.8 (skewed heavy-tailed)
- **AR(2) DGM → AR(1) fit** — model misspecification, evaluated via spectral distance

Sample sizes: T ∈ {50, 75, 100, 200, 300}. Replications: 1 000.  
Alpha grid: seq(0, 1, by = 0.05). Oracle tuning: α* minimises Monte Carlo RMSE.

---

## Repository Structure

```
.
├── Murugesu_2025_HybridLikelihood.pdf   Full thesis
│
├── R/
│   ├── 01_core_estimators.R   Periodogram, spectral density, Whittle / FDEL / Hybrid objectives
│   ├── 02_dgm.R               All data-generating mechanisms
│   ├── 03_drivers.R           Monte Carlo simulation drivers (one per model class)
│   └── 04_summary_tables.R    Summary statistics and LaTeX table generation
│
├── simulations/
│   ├── 01_sim_ar1.Rmd
│   ├── 02_sim_arma11.Rmd
│   ├── 03_sim_arfima.Rmd
│   ├── 04_sim_arch1.Rmd
│   ├── 05_sim_ar1_ao.Rmd
│   ├── 06_sim_ar1_skewt.Rmd
│   └── 07_sim_ar2_misspec.Rmd
│
└── outputs/                   Generated .rds files — gitignored, created on first run
```

---

## Reproducing the Results

**Dependencies**

```r
install.packages(c("pracma", "here", "knitr", "rmarkdown"))
```

**Open the project**

Double-click `thesis_repo.Rproj` in RStudio. This sets the working directory correctly.

**Run a simulation** (example: AR(1))

```r
source("R/01_core_estimators.R")
source("R/02_dgm.R")
source("R/03_drivers.R")

set.seed(2025)
res <- run_hybrid_grid_ar1(
  n = 50, nsim = 1000, phi_true = 0.4,
  alpha_grid = seq(0, 1, by = 0.05),
  innov = "exp1"
)
```

**Or knit the full report**

```r
rmarkdown::render("simulations/01_sim_ar1.Rmd")
```

Results are cached as `.rds` files in `outputs/`. Subsequent runs load the cache
instantly. To rerun from scratch, delete the relevant `.rds` file.

> Runtime: approximately 10–30 minutes per Rmd on a single core.  
> Multi-core supported on Mac/Linux via `ncores = detectCores() - 1`.

---

## Methods

The estimators are all built on the same **spectral estimating function**:

ψⱼ(θ) = (I(λⱼ) / f(λⱼ; θ) − 1) · ∇θ log f(λⱼ; θ)

- **Whittle** imposes these as an M-estimation criterion
- **FDEL** enforces them as empirical moment constraints (Monti, 1997)
- **Hybrid** combines both objectives via convex weighting

The shared structure is what makes the hybrid construction natural —
both methods live in the same frequency-domain estimating-equation framework.

---

## Key References

- Whittle, P. (1953). The analysis of multiple stationary time series. *JRSS-B*, 15, 125–139.
- Monti, A.C. (1997). Empirical likelihood confidence regions in time series. *Biometrika*, 84(2), 395–405.
- Hjort, N.L., McKeague, I.W., and Van Keilegom, I. (2018). Hybrid combinations of parametric and empirical likelihoods. *JRSS-B*, 80(2), 317–352.
- Dæhlen, I. and Hjort, N.L. (2024). Model-robust hybrid likelihood. *arXiv preprint*.

---

## Author

Akshaan Murugesu  
MSc Statistics — University of Geneva  
Supervisors: Prof. Davide La Vecchia and Manon Felix  
November 2025
