# =============================================================================
# Empirical application: forecasting firm leverage (Section 5)
# "Shrinkage and Persistence in Dynamic Panel Forecasting"
# Magdalena Cornejo and Walter Sosa-Escudero
#
# Data: Compustat North America Fundamentals Annual (WRDS). The data are
# proprietary and are not distributed; see README.md for the extraction.
# Input:  data/df0_WRDS.csv
# Output: results/application/ (summary.csv, by_origin.csv, theory_map.csv,
#         table_app_ar1.tex, table_app_xlag.tex)
#
# Design:
#  - Rolling windows of fixed length T (T_SET). For each forecast origin tau,
#    parameters are estimated with years tau-T+1..tau (y in tau-T is the initial
#    lag) and the target is leverage in tau+1. Only information available at tau
#    is used: regressors enter with one lag (x_{t-1}), so forecasts never use
#    contemporaneous size, profitability or tangibility of the target year.
#  - In each window the panel is balanced within the window: firms with complete
#    data for years tau-T..tau+1. Firms therefore enter and exit across windows,
#    which avoids conditioning on survival over the whole 1990-2024 period.
#  - Forecast origins tau = 2010..2023 (targets 2011..2024) for every T, so all
#    window lengths are evaluated on the same years.
#  - Regressors (size, prof, tang) winsorized by year at the 1st/99th percentiles.
#  - Two specifications: "ar1" (the model of the theory) and "xlag"
#    (lev_it = mu + alpha_i + gamma*lev_{i,t-1} + beta'x_{i,t-1} + e_it).
#
# Estimators (same implementations as the Monte Carlo):
#   OLS, FE (LSDV), EBMLE, URE (Kwon, 2026; post-within), AH, GMM1, GMM2
#   (forecasts in levels), RIDGE (closed form), LASSO (glmnet, unpenalized
#   common intercept), and RW (random walk, y_{i,tau}) as a check.
#
# Evaluation:
#   MSE_tau = cross-sectional mean of squared errors at origin tau; reported
#   RMSE = sqrt(mean_tau MSE_tau).
#   Diebold-Mariano tests on the loss differentials d_tau = MSE_tau(A) - MSE_tau(B),
#   with the Harvey, Leybourne and Newbold (1997) small-sample correction. Using
#   one cross-sectional average per origin accounts for cross-sectional
#   dependence of forecast errors within a year (common shocks).
# =============================================================================

rm(list = ls()); options(scipen = 999)

library(dplyr)
library(tidyr)
library(Matrix)
library(glmnet)
library(plm)
library(FEShR)

# -----------------------------------------------------------------------------
# 0) Settings
# -----------------------------------------------------------------------------
T_SET        <- c(5, 10, 20)        # window lengths
ORIGINS      <- 2010:2023           # forecast origins (targets 2011..2024)
SPECS        <- c("ar1", "xlag")
XVARS        <- c("size", "prof", "tang")
GMM_COLLAPSE <- FALSE               # TRUE if the instrument count is a problem
WINS         <- c(0.01, 0.99)         # winsorization of the regressors
DATA_FILE    <- "data/df0_WRDS.csv"
OUT_DIR      <- "results/application"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MODELS <- c("OLS", "FE", "EBMLE", "URE", "AH", "GMM1", "GMM2",
            "RIDGE", "LASSO", "RW")

# -----------------------------------------------------------------------------
# 1) Data
# -----------------------------------------------------------------------------
df0 <- read.csv(DATA_FILE, colClasses = c(gvkey = "character"))

panel <- df0 %>%
  mutate(datadate = as.Date(datadate), fyear = as.integer(fyear),
         sic = as.integer(sich)) %>%
  filter(!is.na(gvkey), !is.na(fyear), !is.na(sic)) %>%
  filter(!(sic >= 6000 & sic <= 6999), !(sic >= 4900 & sic <= 4999)) %>%
  filter(at > 0) %>%
  mutate(lev  = pmin((dltt + dlc) / at, 1),     # leverage capped at one
         size = log(at),
         prof = oibdp / at,                      # operating income before depreciation / assets
         tang = ppent / at) %>%
  filter(is.finite(lev), is.finite(size), is.finite(prof), is.finite(tang)) %>%
  arrange(gvkey, fyear, datadate) %>%
  group_by(gvkey, fyear) %>% slice_tail(n = 1) %>% ungroup() %>%   # one record per firm-year
  dplyr::select(gvkey, fyear, lev, all_of(XVARS)) %>%
  # Winsorize the regressors at the 1st and 99th percentiles of each year's cross-section
  # (standard in corporate finance; uses only information of the same year)
  group_by(fyear) %>%
  mutate(across(all_of(XVARS), ~ pmin(pmax(.x, quantile(.x, WINS[1], na.rm = TRUE)),
                                      quantile(.x, WINS[2], na.rm = TRUE)))) %>%
  ungroup()
rm(df0)

# Wide matrices: rows = firms, columns = years
to_wide <- function(v) {
  w <- panel %>% dplyr::select(gvkey, fyear, all_of(v)) %>%
    pivot_wider(names_from = fyear, values_from = all_of(v)) %>% arrange(gvkey)
  m <- as.matrix(w[, -1]); rownames(m) <- w$gvkey
  m[, order(as.integer(colnames(m))), drop = FALSE]
}
W_lev <- to_wide("lev")
W_x   <- lapply(setNames(XVARS, XVARS), function(v) to_wide(v)[rownames(W_lev), colnames(W_lev)])

# Window for origin tau and length T.
# Y: N x (T+2), columns = years tau-T, ..., tau+1 (column T+2 is the target).
# X: list of N x (T+2) matrices aligned with Y: X[[k]][, j] is the regressor that
#    enters the equation of Y[, j], i.e. x_k in the previous year (x_{t-1}).
make_window <- function(tau, T, spec) {
  yrs  <- as.character((tau - T):(tau + 1))
  ok   <- rowSums(is.na(W_lev[, yrs, drop = FALSE])) == 0
  X    <- NULL
  if (spec == "xlag") {
    xyrs <- as.character((tau - T):tau)          # x needed for t-1 = tau-T..tau
    for (v in XVARS) ok <- ok & rowSums(is.na(W_x[[v]][, xyrs, drop = FALSE])) == 0
    X <- lapply(W_x, function(M) {
      out <- matrix(NA_real_, sum(ok), T + 2)
      out[, 2:(T + 2)] <- M[ok, xyrs, drop = FALSE]  # column j uses year of column j-1
      out
    })
  }
  list(Y = W_lev[ok, yrs, drop = FALSE], X = X, ids = rownames(W_lev)[ok])
}

# -----------------------------------------------------------------------------
# 2) Estimators (generalizations of the Monte Carlo functions to K regressors)
# -----------------------------------------------------------------------------
reg_list <- function(Y, X, Tw) {      # regressors w_it for t = 1..T: lag of y, then x's
  c(list(Y[, 1:Tw, drop = FALSE]),
    if (!is.null(X)) lapply(X, function(M) M[, 2:(Tw + 1), drop = FALSE]))
}

# Closed-form ridge path (Lemma 1; Appendix C)
ridge_prep <- function(Y, X, Tw) {
  yy <- Y[, 2:(Tw + 1), drop = FALSE]; W <- reg_list(Y, X, Tw); K <- length(W)
  ybar_i <- rowMeans(yy); wbar_i <- sapply(W, rowMeans)
  if (is.null(dim(wbar_i))) wbar_i <- matrix(wbar_i, ncol = K)
  ybar <- mean(ybar_i); wbar <- colMeans(wbar_i)
  Wd <- lapply(W, function(M) M - rowMeans(M)); yd <- yy - ybar_i
  SWww <- matrix(0, K, K); SWwy <- numeric(K)
  for (a in 1:K) { SWwy[a] <- sum(Wd[[a]] * yd); for (b in 1:K) SWww[a, b] <- sum(Wd[[a]] * Wd[[b]]) }
  Bw <- sweep(wbar_i, 2, wbar); By <- ybar_i - ybar
  list(SWww = SWww, SWwy = SWwy, SBww = crossprod(Bw), SBwy = as.numeric(crossprod(Bw, By)),
       Bw = Bw, By = By, ybar = ybar, wbar = wbar, wbar_i = wbar_i, ybar_i = ybar_i,
       Wd = Wd, yd = yd, Tw = Tw, K = K)
}
ridge_coef <- function(prep, kappa) {   # kappa = 1: FE; kappa = 0: pooled OLS
  om <- prep$Tw * (1 - kappa)
  theta <- as.numeric(solve(prep$SWww + om * prep$SBww, prep$SWwy + om * prep$SBwy))
  mu <- prep$ybar - sum(prep$wbar * theta)
  list(a = mu + kappa * (prep$By - as.numeric(prep$Bw %*% theta)), theta = theta)
}

# Forecast of the target (column T+2) given unit intercepts and slopes
forecast_target <- function(Y, X, a, theta, Tw) {
  f <- a + theta[1] * Y[, Tw + 1]
  if (!is.null(X)) for (k in seq_along(X)) f <- f + theta[k + 1] * X[[k]][, Tw + 2]
  f
}

# Unit intercepts in levels given slopes
levels_intercepts <- function(Y, X, Tw, theta) {
  r <- Y[, 2:(Tw + 1), drop = FALSE] - theta[1] * Y[, 1:Tw, drop = FALSE]
  if (!is.null(X)) for (k in seq_along(X)) r <- r - theta[k + 1] * X[[k]][, 2:(Tw + 1), drop = FALSE]
  rowMeans(r)
}

# Anderson-Hsiao: first differences, t = 2..T, instrument y_{t-2} (x's instrument themselves)
ah_fit <- function(Y, X, Tw) {
  ts <- 2:Tw
  dy <- as.vector(Y[, ts + 1] - Y[, ts])
  D  <- cbind(as.vector(Y[, ts] - Y[, ts - 1])); Z <- cbind(as.vector(Y[, ts - 1]))
  if (!is.null(X)) for (M in X) {
    dx <- as.vector(M[, ts + 1] - M[, ts]); D <- cbind(D, dx); Z <- cbind(Z, dx)
  }
  as.numeric(solve(crossprod(Z, D), crossprod(Z, dy)))
}

# Arellano-Bond via pgmm, forecasts in levels
gmm_fit <- function(Y, X, Tw, step) {
  N <- nrow(Y)
  df <- data.frame(id = rep(seq_len(N), times = Tw + 1), time = rep(0:Tw, each = N),
                   lev = as.vector(Y[, 1:(Tw + 1)]))
  xn <- character(0)
  if (!is.null(X)) for (k in seq_along(X)) {
    M <- X[[k]][, 1:(Tw + 1)]; M[, 1] <- M[, 2]   # x of the initial year is never used
    nm <- paste0("x", k); df[[nm]] <- as.vector(M); xn <- c(xn, nm)
  }
  rhs <- paste(c("lag(lev, 1)", xn), collapse = " + ")
  f <- if (length(xn)) as.formula(paste("lev ~", rhs, "| lag(lev, 2:99) |", paste(xn, collapse = " + ")))
       else as.formula(paste("lev ~", rhs, "| lag(lev, 2:99)"))
  g <- pgmm(f, data = pdata.frame(df, index = c("id", "time")), effect = "individual",
            model = step, transformation = "d", collapse = GMM_COLLAPSE)
  b <- coef(g)
  as.numeric(c(b[grep("^lag\\(lev", names(b))][1], b[xn]))
}

# Post-within shrinkage (FEShR). Intercept variances computed analytically:
# Var(a_i) = s2/T + wbar_i' V wbar_i, with V = s2 * SWww^{-1} (same as lm with dummies).
shrink_fe <- function(prep, type) {
  set.seed(1)                     # FEShR may use random starting values
  cf <- ridge_coef(prep, 1)
  res <- prep$yd; for (k in 1:prep$K) res <- res - cf$theta[k] * prep$Wd[[k]]
  N <- length(prep$ybar_i)
  s2 <- sum(res^2) / (N * prep$Tw - N - prep$K)
  V  <- s2 * solve(prep$SWww)
  v  <- s2 / prep$Tw + rowSums((prep$wbar_i %*% V) * prep$wbar_i)
  a  <- prep$ybar_i - as.numeric(prep$wbar_i %*% cf$theta)
  fit <- FEShR::fe_shrink(y = matrix(a, nrow = 1), M = lapply(v, function(x) matrix(x, 1, 1)),
                          centering = "gen", type = type, n_init_vals = 1,
                          all_init_vals = FALSE, optim_control = list(maxit = 500))
  list(a = as.numeric(fit$thetahat), theta = cf$theta)
}

# Corrected time-series cross-validation. Folds: train on t = 1..s, validate on
# t = s + 1, for s = max(2, ceiling(T/2))..T-1 (one-step ahead, expanding window);
# minimum rule.
# RIDGE: closed-form path (Lemma 1), grid over kappa, the same kappa in every fold.
KGRID <- sort(unique(c(0, 1, 1 - 10^seq(-3, 0, length.out = 40))))
ridge_cv <- function(Y, X, Tw, kgrid = KGRID) {
  L <- sapply(max(2, ceiling(Tw / 2)):(Tw - 1), function(s) {
    pr <- ridge_prep(Y, X, s)
    sapply(kgrid, function(k) {
      cf <- ridge_coef(pr, k)
      mean((Y[, s + 2] - forecast_target(Y, X, cf$a, cf$theta, s))^2)
    })
  })
  L <- matrix(L, nrow = length(kgrid))
  k <- kgrid[which.min(rowMeans(L))]
  c(ridge_coef(ridge_prep(Y, X, Tw), k), list(kappa = k))
}
# LASSO: glmnet, unit dummies penalized; common intercept, lag and x's unpenalized;
# standardize = FALSE so that a given lambda implies the same soft-threshold in every fold.
lasso_cv <- function(Y, X, Tw) {
  N <- nrow(Y); W <- reg_list(Y, X, Tw)
  y <- as.vector(Y[, 2:(Tw + 1)]); times <- rep(1:Tw, each = N)
  D <- sparseMatrix(i = seq_len(N * Tw), j = rep(seq_len(N), times = Tw), x = 1, dims = c(N * Tw, N))
  Xm <- cbind(Matrix(sapply(W, as.vector), sparse = TRUE), D)
  nth <- length(W); pf <- c(rep(0, nth), rep(1, N))
  lam <- glmnet(Xm, y, alpha = 1, penalty.factor = pf, intercept = TRUE, standardize = FALSE)$lambda
  L <- sapply(max(2, ceiling(Tw / 2)):(Tw - 1), function(s) {
    tr <- times <= s; va <- times == s + 1
    fit <- glmnet(Xm[tr, , drop = FALSE], y[tr], alpha = 1, lambda = lam, penalty.factor = pf,
                  intercept = TRUE, standardize = FALSE)
    colMeans((y[va] - as.matrix(predict(fit, newx = Xm[va, , drop = FALSE])))^2)
  })
  L <- matrix(L, nrow = length(lam))
  b <- as.numeric(coef(glmnet(Xm, y, alpha = 1, lambda = lam[which.min(rowMeans(L))],
                              penalty.factor = pf, intercept = TRUE, standardize = FALSE)))
  list(a = b[1] + b[(2 + nth):length(b)], theta = b[2:(1 + nth)])
}

# -----------------------------------------------------------------------------
# 3) Rolling-window forecasts
# -----------------------------------------------------------------------------
run_window <- function(tau, Tw, spec) {
  w <- make_window(tau, Tw, spec); Y <- w$Y; X <- w$X
  fits <- list()
  try_fit <- function(name, expr) {
    r <- try(expr, silent = TRUE)
    if (inherits(r, "try-error")) message("  ", name, " failed: ", conditionMessage(attr(r, "condition")))
    else fits[[name]] <<- r
  }
  prep <- ridge_prep(Y, X, Tw)
  try_fit("OLS", ridge_coef(prep, 0))
  try_fit("FE",  ridge_coef(prep, 1))
  try_fit("EBMLE", shrink_fe(prep, "EBMLE"))
  try_fit("URE",   shrink_fe(prep, "URE"))
  try_fit("AH", { th <- ah_fit(Y, X, Tw); list(a = levels_intercepts(Y, X, Tw, th), theta = th) })
  for (st in c("onestep", "twosteps"))
    try_fit(if (st == "onestep") "GMM1" else "GMM2",
            { th <- gmm_fit(Y, X, Tw, st); list(a = levels_intercepts(Y, X, Tw, th), theta = th) })
  try_fit("RIDGE", ridge_cv(Y, X, Tw))
  try_fit("LASSO", lasso_cv(Y, X, Tw))

  target <- Y[, Tw + 2]
  err <- sapply(MODELS, function(m) {
    if (m == "RW") return(target - Y[, Tw + 1])
    if (is.null(fits[[m]])) return(rep(NA_real_, length(target)))
    target - forecast_target(Y, X, fits[[m]]$a, fits[[m]]$theta, Tw)
  })
  gam <- sapply(MODELS, function(m) if (is.null(fits[[m]])) NA_real_ else fits[[m]]$theta[1])
  gam["RW"] <- 1

  list(err = err, gamma = gam, N = nrow(Y),
       kappa_ridge = if (!is.null(fits$RIDGE)) fits$RIDGE$kappa else NA_real_)
}

results <- list()
for (spec in SPECS) for (Tw in T_SET) for (tau in ORIGINS) {
  key <- paste(spec, Tw, tau, sep = "_")
  f <- file.path(OUT_DIR, paste0("win_", key, ".rds"))
  if (file.exists(f)) { results[[key]] <- readRDS(f); next }
  message("spec = ", spec, " | T = ", Tw, " | origin = ", tau)
  r <- run_window(tau, Tw, spec); r$spec <- spec; r$T <- Tw; r$tau <- tau
  saveRDS(r, f); results[[key]] <- r
}

# -----------------------------------------------------------------------------
# 4) Evaluation
# -----------------------------------------------------------------------------
by_origin <- do.call(rbind, lapply(results, function(r) {
  data.frame(spec = r$spec, T = r$T, tau = r$tau, N = r$N, model = MODELS,
             mse = colMeans(r$err[, MODELS]^2), mae = colMeans(abs(r$err[, MODELS])),
             bias = colMeans(r$err[, MODELS]), gamma = as.numeric(r$gamma[MODELS]))
}))
write.csv(by_origin, file.path(OUT_DIR, "by_origin.csv"), row.names = FALSE)

summary_tab <- by_origin %>%
  group_by(spec, T, model) %>%
  summarise(RMSE = sqrt(mean(mse)), MAE = mean(mae), Bias = mean(bias),
            gamma = mean(gamma), N = mean(N), .groups = "drop") %>%
  group_by(spec, T) %>%
  mutate(rel_MSE = RMSE^2 / RMSE[model == "FE"]^2,
         rel_RMSE = RMSE / RMSE[model == "FE"]) %>% ungroup()

# Diebold-Mariano on per-origin loss differentials, HLN correction (h = 1)
dm_test <- function(a, b) {
  d <- a - b; n <- length(d)
  stat <- mean(d) / sqrt(var(d) / n) * sqrt((n - 1) / n)
  c(stat = stat, pval = 2 * pt(-abs(stat), df = n - 1))
}
dm_tab <- function(benchmark) {
  by_origin %>% group_by(spec, T) %>% group_modify(function(d, k) {
    bm <- d %>% filter(model == benchmark) %>% arrange(tau)
    do.call(rbind, lapply(setdiff(MODELS, benchmark), function(m) {
      x <- d %>% filter(model == m) %>% arrange(tau)
      ok <- is.finite(x$mse) & is.finite(bm$mse)
      t <- if (sum(ok) > 2) dm_test(x$mse[ok], bm$mse[ok]) else c(stat = NA, pval = NA)
      data.frame(model = m, benchmark = benchmark, dm_stat = t[["stat"]], dm_pval = t[["pval"]])
    }))
  }) %>% ungroup()
}
dm_fe  <- dm_tab("FE")
dm_ols <- dm_tab("OLS")

summary_tab <- summary_tab %>%
  left_join(dm_fe %>% dplyr::select(spec, T, model, dm_pval_FE = dm_pval), by = c("spec", "T", "model")) %>%
  left_join(dm_ols %>% dplyr::select(spec, T, model, dm_pval_OLS = dm_pval), by = c("spec", "T", "model"))
write.csv(summary_tab, file.path(OUT_DIR, "summary.csv"), row.names = FALSE)

cat("\n=== Relative MSE (vs FE) ===\n")
print(summary_tab %>% dplyr::select(spec, T, model, rel_MSE) %>%
        pivot_wider(names_from = model, values_from = rel_MSE) %>%
        mutate(across(-c(spec, T), ~ round(.x, 3))), width = 200)
cat("\n=== RMSE ===\n")
print(summary_tab %>% dplyr::select(spec, T, model, RMSE) %>%
        pivot_wider(names_from = model, values_from = RMSE) %>%
        mutate(across(-c(spec, T), ~ round(.x, 4))), width = 200)
kap <- do.call(rbind, lapply(results, function(r) data.frame(spec = r$spec, T = r$T, tau = r$tau,
                                                           kappa = r$kappa_ridge)))
cat("\n=== RIDGE: kappa selected by cross-validation (average over origins) ===\n")
print(aggregate(kappa ~ spec + T, data = kap, FUN = mean))
cat("\n=== Average estimated gamma ===\n")
print(summary_tab %>% dplyr::select(spec, T, model, gamma) %>%
        pivot_wider(names_from = model, values_from = gamma) %>%
        mutate(across(-c(spec, T), ~ round(.x, 3))), width = 200)

# -----------------------------------------------------------------------------
# 5) Where the application sits in the theory (AR(1) specification)
#    Average calibrated parameters and implied limiting risks, relative to FE.
# -----------------------------------------------------------------------------
consts_T <- function(g, T) {
  k <- seq_len(T - 1)
  list(v = 1 + (2 / T) * sum((T - k) * g^k), h = sum((1 - g^(T - k)) / (1 - g)),
       c = g * (1 - g^T) / (T * (1 - g)), p = (1 - g^T) / (T * (1 - g)))
}
bias_T <- function(om, g, T, s2, sa2) {
  k <- consts_T(g, T); su2 <- s2 / (1 - g^2)
  (-s2 * k$h / T + om * (sa2 / (1 - g) + s2 * k$h / T^2)) /
    (su2 * (T - k$v) + om * (sa2 / (1 - g)^2 + su2 * k$v / T))
}
risk_T <- function(kappa, gf, g, T, s2, sa2) {
  k <- consts_T(g, T); su2 <- s2 / (1 - g^2); d <- g - gf
  s2 + kappa^2 * s2 / T + (1 - kappa)^2 * (1 - gf)^2 * sa2 / (1 - g)^2 +
    d^2 * su2 * (1 - 2 * kappa * k$c + kappa^2 * k$v / T) - 2 * d * kappa * s2 * (k$p - kappa * k$h / T^2)
}
# Calibration in each AR(1) window. gamma: bias-corrected within estimator, i.e. the
# solution in g of gamma_W = g + b_T(0; g) (Theorem 1 with omega = 0; this bias does
# not depend on sigma2_alpha). sigma2 from first differences; sigma2_alpha =
# Var(ybar_i - g*ybar_{i,-1}) - sigma2/T.
within_plim <- function(g, T) g + bias_T(0, g, T, 1 - g^2, 1)   # s2/su2 = 1 - g^2
gamma_bc <- function(gW, T) {
  f <- function(g) within_plim(g, T) - gW
  lo <- -0.99; hi <- 0.999
  if (f(hi) < 0) return(hi)
  uniroot(f, c(lo, hi), tol = 1e-10)$root
}
calibrate <- function(tau, Tw) {
  Y <- make_window(tau, Tw, "ar1")$Y
  gW <- ridge_coef(ridge_prep(Y, NULL, Tw), 1)$theta[1]
  g <- gamma_bc(gW, Tw); ts <- 2:Tw
  s2 <- mean(((Y[, ts + 1] - Y[, ts]) - g * (Y[, ts] - Y[, ts - 1]))^2) / 2
  sa2 <- max(var(rowMeans(Y[, 2:(Tw + 1)] - g * Y[, 1:Tw])) - s2 / Tw, 0)
  list(gamma = g, sigma2 = s2, sigma2_alpha = sa2, phi = sa2 / s2)
}
theory_map <- do.call(rbind, lapply(results, function(r) {
  if (r$spec != "ar1") return(NULL)
  p <- calibrate(r$tau, r$T); g <- p$gamma; s2 <- p$sigma2; sa2 <- p$sigma2_alpha; Tw <- r$T
  su2 <- s2 / (1 - g^2); om2 <- sa2 * su2 / (sa2 / (1 - g)^2 + su2)
  RFE <- risk_T(1, g + bias_T(0, g, Tw, s2, sa2), g, Tw, s2, sa2)
  data.frame(T = Tw, tau = r$tau, gamma = g, phi = p$phi,
             R_OLS = risk_T(0, g + bias_T(Tw, g, Tw, s2, sa2), g, Tw, s2, sa2) / RFE,
             R_K   = (s2 + sa2 * s2 / (Tw * sa2 + s2)) / RFE,
             R_STAR = (s2 + om2 * (s2 / Tw) / (om2 + s2 / Tw)) / RFE)
}))
if (!is.null(theory_map)) {
  tm <- theory_map %>% group_by(T) %>% summarise(across(-tau, mean), .groups = "drop")
  write.csv(tm, file.path(OUT_DIR, "theory_map.csv"), row.names = FALSE)
  cat("\n=== Calibrated theory (AR(1)): limits relative to FE ===\n"); print(round(as.data.frame(tm), 3))
}

# -----------------------------------------------------------------------------
# 6) LaTeX table: RMSE, relative MSE and DM p-values, for one spec and all T
# -----------------------------------------------------------------------------
latex_app <- function(spec_value, models = c("OLS", "FE", "EBMLE", "URE", "AH", "GMM1",
                                             "RIDGE", "LASSO", "RW")) {
  s <- summary_tab %>% filter(spec == spec_value, model %in% models)
  fmt <- function(x, d) formatC(x, format = "f", digits = d)
  star <- function(p) ifelse(is.na(p), "", ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$",
                       ifelse(p < 0.10, "$^{*}$", ""))))
  lines <- c("\\begin{table}[ht]", "\\centering", "\\small",
             paste0("\\caption{Out-of-sample forecast accuracy, firm leverage (", spec_value, ")}"),
             paste0("\\begin{tabular}{l", strrep("cc", length(T_SET)), "}"), "\\toprule",
             paste0(" & ", paste0("\\multicolumn{2}{c}{$T=", T_SET, "$}", collapse = " & "), " \\\\"),
             paste0("Estimator", strrep(" & RMSE & Rel. MSE", length(T_SET)), " \\\\"), "\\midrule")
  for (m in models) {
    cells <- sapply(T_SET, function(Tw) {
      r <- s %>% filter(T == Tw, model == m)
      if (nrow(r) == 0) return(" & ")
      paste0(fmt(r$RMSE, 4), " & ", fmt(r$rel_MSE, 3), if (m != "FE") star(r$dm_pval_FE) else "")
    })
    lines <- c(lines, paste0(gsub("_", "\\\\_", m), " & ", paste(cells, collapse = " & "), " \\\\"))
  }
  c(lines, "\\bottomrule", "\\end{tabular}",
    "\\begin{minipage}{0.95\\textwidth}\\footnotesize\\vspace{0.4em}",
    paste0("Notes: Rolling windows of length $T$; forecast origins ", min(ORIGINS), "--", max(ORIGINS),
           ". RMSE is the square root of the average over origins of the cross-sectional MSE. ",
           "Rel.\\ MSE is relative to fixed effects. Stars: Diebold--Mariano test of equal accuracy ",
           "with fixed effects (HLN correction), $^{*}$10\\%, $^{**}$5\\%, $^{***}$1\\%."),
    "\\end{minipage}", "\\end{table}")
}
for (sp in SPECS) writeLines(latex_app(sp), file.path(OUT_DIR, paste0("table_app_", sp, ".tex")))
