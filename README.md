# Shrinkage and Persistence in Dynamic Panel Forecasting

Replication code for

> Cornejo, M. and W. Sosa-Escudero. *Shrinkage and Persistence in Dynamic Panel Forecasting*.

The paper is in [`paper/`](paper/Cornejo_SosaEscudero_Shrinkage_Persistence_DPD.pdf).

The repository contains the code for the Monte Carlo experiments (Section 4 and
the appendix), the empirical application to firm leverage (Section 5) and the two
figures. All the code is in R. The Compustat data used in the application are proprietary and are **not**
included; Section [Data](#data) explains how to obtain them.

## Structure

```
paper/                  the paper (PDF)
code/
  montecarlo.R          Monte Carlo experiments: Tables 1, 4 and 5
  application.R         Empirical application: Tables 2 and 3
  figures.R             Figures 1 and 2
data/                   (not distributed) place df0_WRDS.csv here
results/                created by the scripts
```

All scripts are written in R and are run from the root of the repository.

## Requirements

**R** (version 4.2 or later) with the packages `dplyr`, `tidyr`, `Matrix`, `glmnet`,
`plm` and `FEShR`. `FEShR` implements the empirical Bayes and unbiased risk
shrinkage estimators of Kwon (2026) and is installed from GitHub:

```r
install.packages("remotes")
remotes::install_github("soonwookwon/FEShR")
```

If compilation fails with a C++11 error, download the package, replace
`CXX_STD = CXX11` by `CXX_STD = CXX14` in `src/Makevars` and `src/Makevars.win`,
and install it with `remotes::install_local()`.

## Monte Carlo experiments

```r
source("code/montecarlo.R")
```

The script runs the 48 designs (N = 20, 100; T = 5, 10, 20, 30; gamma = 0.2, 0.5,
0.8; phi = 0.25, 1, 4) with 1,000 replications each and stores one `.rds` file per
design in `results/montecarlo/`. Designs already computed are skipped, so the run
can be interrupted and resumed. The summary file `results/montecarlo/mc_results.csv`
contains, for each design and forecast, the mean squared forecast error (`mse`), the
ratio to fixed effects (`rel`) and its Monte Carlo standard error (`se`), together
with the theoretical limits `R_K`, `R_STAR`, `R_FE` and `R_OLS`.

The full grid takes several hours, mostly because of the GMM estimators. Setting
`skip_slow = TRUE` in `run_grid()` computes only OLS, fixed effects, RIDGE and LASSO.
Seeds are fixed, so the simulated samples are the same in both cases.

## Empirical application

```r
source("code/application.R")
```

The script reads `data/df0_WRDS.csv`, builds rolling windows of length T = 5, 10
and 20 with forecast origins 2010–2023, and computes the out-of-sample forecasts of
all methods for the AR(1) specification and for the specification with lagged firm
characteristics. Each window is cached in `results/application/`. The outputs are

- `summary.csv`: RMSE, MSE relative to fixed effects and Diebold–Mariano p-values
  (against fixed effects and against pooled OLS), Table 2;
- `table_app_ar1.tex`, `table_app_xlag.tex`: LaTeX versions of Table 2;
- `theory_map.csv`: calibrated parameters and limiting risks, Table 3;
- `by_origin.csv`: cross-sectional MSE of each method at each forecast origin.

## Figures

```r
source("code/figures.R")
```

Figure 1 is computed from the closed-form limits of Section 3. Figure 2 combines
these limits with a simulation (1,000 replications, a few minutes). Both figures
are saved in `results/figures/`. The script uses only base R.

## Data

The application uses Compustat North America Fundamentals Annual (`comp.funda`),
available through Wharton Research Data Services (WRDS) to subscribing
institutions. The file `data/df0_WRDS.csv` must contain one row per firm and fiscal
year, for fiscal years 1990–2024, with the following variables:

| Variable   | Description                                   |
|------------|-----------------------------------------------|
| `gvkey`    | Firm identifier                               |
| `datadate` | Date of the fiscal year end                   |
| `fyear`    | Fiscal year                                   |
| `sich`     | Historical SIC code                           |
| `at`       | Total assets                                  |
| `dltt`     | Long-term debt                                |
| `dlc`      | Debt in current liabilities                   |
| `oibdp`    | Operating income before depreciation          |
| `ppent`    | Property, plant and equipment (net)           |

The standard WRDS screens for annual industrial data are `indfmt = 'INDL'`,
`datafmt = 'STD'`, `popsrc = 'D'` and `consol = 'C'`. The script excludes
financial firms (SIC 6000–6999) and regulated utilities (SIC 4900–4999), keeps one
record per firm-year and constructs leverage, size, profitability and tangibility
as described in Section 5 of the paper.

## Contact

Magdalena Cornejo, Universidad Torcuato Di Tella (mcornejo@utdt.edu).
