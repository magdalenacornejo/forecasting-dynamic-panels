# =============================================================================
# Monte Carlo experiments (Section 4 and appendix)
# "Shrinkage and Persistence in Dynamic Panel Forecasting"
# Magdalena Cornejo and Walter Sosa-Escudero
#
# Forecasts (one-step forecast of y_{i,T+1}; parameters estimated with t = 1..T):
#   OLS (pooled), LSDV (fixed effects)          endpoints of the ridge path
#   RIDGE, LASSO                                 penalized unit effects, unpenalized
#                                                common intercept; penalty chosen by
#                                                time-series cross-validation
#   EBMLE, URE (Kwon, 2026)                      shrinkage of the LSDV intercepts
#   AH, GMM1, GMM2 (Arellano-Bond)               consistent slope, forecasts in levels
#
# Theoretical limits reported with each design (no regressors):
#   R_K     known-persistence benchmark (Proposition 1)
#   R_STAR  best linear forecast (Proposition 4)
#   R_FE, R_OLS  limits of fixed effects and pooled OLS (Proposition 2, Theorem 1)
#
# DGP: y_it = mu + alpha_i + gamma*y_{i,t-1} [+ beta*x_it] + eps_it, Var(alpha_i) = phi*s2.
#   init = "stationary" (used in the paper): y_i0 drawn from its stationary
#   distribution given alpha_i. Other options: "uninformative", "zero" (burn-in).
#
# Cross-validation (RIDGE and LASSO): train on t = 1..s, validate the one-step
# forecast of t = s + 1, for s = max(2, ceiling(T/2))..T-1; the penalty minimizing
# the average validation loss is selected. RIDGE uses the closed-form path
# (Lemma 1) and a grid over kappa = T/(T + lambda) that is dense near one.
#
# Requires: plm, glmnet, FEShR (https://github.com/soonwookwon/FEShR).
# Output: one .rds file per design in results/montecarlo/ and the summary file
# results/montecarlo/mc_results.csv.
# =============================================================================

rm(list = ls()); options(scipen = 999)

library(plm)
library(glmnet)
library(FEShR)

# =============================================================================
# 1. Theory: constants, bias function and limiting risk (Section 3)
# =============================================================================
consts_T <- function(g, T) {
  k <- seq_len(T - 1)
  list(v = 1 + (2 / T) * sum((T - k) * g^k),
       h = sum((1 - g^(T - k)) / (1 - g)),
       c = g * (1 - g^T) / (T * (1 - g)),
       p = (1 - g^T) / (T * (1 - g)))
}

# Asymptotic bias of the ridge estimator of gamma, b_T(omega) (Theorem 1)
bias_T <- function(omega, g, T, s2, sa2) {
  k <- consts_T(g, T); su2 <- s2 / (1 - g^2)
  num <- -s2 * k$h / T + omega * (sa2 / (1 - g) + s2 * k$h / T^2)
  den <- su2 * (T - k$v) + omega * (sa2 / (1 - g)^2 + su2 * k$v / T)
  num / den
}

# Limiting risk R_T(kappa, g_forecast) (Theorem 2)
risk_T <- function(kappa, gf, g, T, s2, sa2) {
  k <- consts_T(g, T); su2 <- s2 / (1 - g^2); d <- g - gf
  s2 + kappa^2 * s2 / T + (1 - kappa)^2 * (1 - gf)^2 * sa2 / (1 - g)^2 +
    d^2 * su2 * (1 - 2 * kappa * k$c + kappa^2 * k$v / T) -
    2 * d * kappa * s2 * (k$p - kappa * k$h / T^2)
}

# Theoretical limits (N -> infinity, T fixed) for a design without regressors.
# R_STAR: best linear forecast; omega2 = Var(alpha_i | y_i0) (= sa2 if uninformative).
theory_limits <- function(gamma, T, s2, sa2, init) {
  su2 <- s2 / (1 - gamma^2)
  om2 <- if (init == "stationary") sa2 * su2 / (sa2 / (1 - gamma)^2 + su2) else sa2
  out <- c(R_K    = s2 + sa2 * s2 / (T * sa2 + s2),
           R_STAR = s2 + om2 * (s2 / T) / (om2 + s2 / T))
  if (init == "stationary") {
    out <- c(out,
             R_FE  = risk_T(1, gamma + bias_T(0, gamma, T, s2, sa2), gamma, T, s2, sa2),
             R_OLS = risk_T(0, gamma + bias_T(T, gamma, T, s2, sa2), gamma, T, s2, sa2))
  }
  out
}

# =============================================================================
# 2. Data generating process
#    Y, X: N x (T0 + 1) matrices, column j corresponds to t = j - 1
# =============================================================================
simulate_panel <- function(N, T0, gamma, mu, phi, s2, sparse, include_x,
                           beta, rho_x, sd_xi, init, burn) {
  alpha <- rnorm(N, 0, sqrt(phi * s2))
  n_zero <- round(sparse * N)
  if (n_zero > 0) alpha[sample(N, n_zero)] <- 0

  if (!include_x && init == "stationary") {
    # Exact draw from the stationary distribution (Assumption 1(iv))
    Y <- matrix(NA_real_, N, T0 + 1)
    Y[, 1] <- (mu + alpha) / (1 - gamma) + rnorm(N, 0, sqrt(s2 / (1 - gamma^2)))
    for (j in 2:(T0 + 1)) Y[, j] <- mu + alpha + gamma * Y[, j - 1] + rnorm(N, 0, sqrt(s2))
    return(list(Y = Y, X = NULL, alpha = alpha))
  }

  if (!include_x && init == "uninformative") {
    # y_i0 independent of alpha_i, with the same marginal variance as in the stationary case
    Y <- matrix(NA_real_, N, T0 + 1)
    Y[, 1] <- mu / (1 - gamma) + rnorm(N, 0, sqrt(s2 / (1 - gamma^2) + phi * s2 / (1 - gamma)^2))
    for (j in 2:(T0 + 1)) Y[, j] <- mu + alpha + gamma * Y[, j - 1] + rnorm(N, 0, sqrt(s2))
    return(list(Y = Y, X = NULL, alpha = alpha))
  }

  # Burn-in: "stationary" starts at the long-run mean (approximate stationarity
  # with regressors); "zero" starts at y_0 = x_0 = 0 and discards a burn-in.
  Ttot <- burn + T0 + 1
  Y <- matrix(NA_real_, N, Ttot)
  X <- if (include_x) matrix(NA_real_, N, Ttot) else NULL
  if (init == "stationary") {
    Y[, 1] <- (mu + alpha) / (1 - gamma)
    if (include_x) X[, 1] <- rnorm(N, 0, sd_xi / sqrt(1 - rho_x^2))
  } else {
    Y[, 1] <- 0
    if (include_x) X[, 1] <- 0
  }
  for (j in 2:Ttot) {
    xb <- 0
    if (include_x) {
      X[, j] <- rho_x * X[, j - 1] + rnorm(N, 0, sd_xi)
      xb <- beta * X[, j]
    }
    Y[, j] <- mu + alpha + gamma * Y[, j - 1] + xb + rnorm(N, 0, sqrt(s2))
  }
  keep <- (burn + 1):Ttot
  list(Y = Y[, keep], X = if (include_x) X[, keep] else NULL, alpha = alpha)
}

# =============================================================================
# 3. Estimators
# =============================================================================

# Closed-form ridge path (Lemma 1; Appendix C with regressors)
ridge_prep <- function(Y, X, Tstar) {
  yy <- Y[, 2:(Tstar + 1), drop = FALSE]
  W  <- list(Y[, 1:Tstar, drop = FALSE])
  if (!is.null(X)) W[[2]] <- X[, 2:(Tstar + 1), drop = FALSE]
  K <- length(W)
  ybar_i <- rowMeans(yy)
  wbar_i <- sapply(W, rowMeans); if (is.null(dim(wbar_i))) wbar_i <- matrix(wbar_i, ncol = 1)
  ybar <- mean(ybar_i); wbar <- colMeans(wbar_i)
  Wd <- lapply(W, function(M) M - rowMeans(M))
  yd <- yy - ybar_i
  SWww <- matrix(0, K, K); SWwy <- numeric(K)
  for (a in 1:K) {
    SWwy[a] <- sum(Wd[[a]] * yd)
    for (b in 1:K) SWww[a, b] <- sum(Wd[[a]] * Wd[[b]])
  }
  Bw <- sweep(wbar_i, 2, wbar); By <- ybar_i - ybar
  list(SWww = SWww, SWwy = SWwy, SBww = crossprod(Bw), SBwy = as.numeric(crossprod(Bw, By)),
       Bw = Bw, By = By, ybar = ybar, wbar = wbar, Tstar = Tstar)
}

# kappa = 1: fixed effects; kappa = 0: pooled OLS; kappa = T/(T + lambda) in between
ridge_coef <- function(prep, kappa) {
  om <- prep$Tstar * (1 - kappa)
  theta <- as.numeric(solve(prep$SWww + om * prep$SBww, prep$SWwy + om * prep$SBwy))
  mu <- prep$ybar - sum(prep$wbar * theta)
  alpha <- kappa * (prep$By - as.numeric(prep$Bw %*% theta))
  list(a = mu + alpha, theta = theta)          # a_i = unit-specific intercept
}

# Anderson-Hsiao, equations t = 2..Tstar, instrument y_{t-2}
ah_fit <- function(Y, X, Tstar) {
  ts <- 2:Tstar
  dy  <- as.vector(Y[, ts + 1] - Y[, ts])
  D <- cbind(as.vector(Y[, ts] - Y[, ts - 1]))
  Z <- cbind(as.vector(Y[, ts - 1]))
  if (!is.null(X)) {
    dx <- as.vector(X[, ts + 1] - X[, ts]); D <- cbind(D, dx); Z <- cbind(Z, dx)
  }
  theta <- as.numeric(solve(crossprod(Z, D), crossprod(Z, dy)))
  res <- dy - as.numeric(D %*% theta)
  list(theta = theta, s2 = mean(res^2) / 2)
}

# Unit intercepts in levels given slope estimates: a_i = mean_t (y - g*y_{-1} - b*x)
levels_intercepts <- function(Y, X, Tstar, theta) {
  r <- Y[, 2:(Tstar + 1), drop = FALSE] - theta[1] * Y[, 1:Tstar, drop = FALSE]
  if (!is.null(X)) r <- r - theta[2] * X[, 2:(Tstar + 1), drop = FALSE]
  rowMeans(r)
}

# One-step forecast errors for t = t_from..t_to, given a_i and theta
forecast_errors <- function(Y, X, a, theta, t_from, t_to) {
  cols <- (t_from:t_to) + 1
  yhat <- a + theta[1] * Y[, cols - 1, drop = FALSE]
  if (!is.null(X)) yhat <- yhat + theta[2] * X[, cols, drop = FALSE]
  Y[, cols, drop = FALSE] - yhat
}

# Post-within shrinkage with FEShR (EBMLE / URE of Kwon, 2026)
shrink_fe <- function(a_hat, v_hat, type) {
  M_list <- lapply(v_hat, function(v) matrix(v, 1, 1))
  fit <- try(FEShR::fe_shrink(y = matrix(a_hat, nrow = 1), M = M_list, centering = "gen",
                              type = type, n_init_vals = 1, all_init_vals = FALSE,
                              optim_control = list(maxit = 500)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  as.numeric(fit$thetahat)
}

lsdv_intercepts_var <- function(df, include_x) {
  f <- if (include_x) Y ~ lagY + X + factor(id) else Y ~ lagY + factor(id)
  fit <- lm(f, data = df)
  cf <- coef(fit); V <- vcov(fit)
  ids <- levels(factor(df$id)); nm <- paste0("factor(id)", ids)
  a <- numeric(length(ids)); v <- numeric(length(ids))
  a[1] <- cf["(Intercept)"]; v[1] <- V["(Intercept)", "(Intercept)"]
  for (j in 2:length(ids)) {
    a[j] <- cf["(Intercept)"] + cf[nm[j]]
    v[j] <- V["(Intercept)", "(Intercept)"] + V[nm[j], nm[j]] + 2 * V["(Intercept)", nm[j]]
  }
  list(a = a, v = v)
}

# Corrected time-series cross-validation (expanding window, one-step ahead).
# Folds: train on t = 1..s, validate on t = s + 1, for s = s_min..Tstar - 1,
# with s_min = max(2, ceiling(Tstar / 2)).
# RIDGE: closed-form path, grid over kappa (the same kappa in every fold and in
# the final fit). The minimum rule is used.
KGRID <- sort(unique(c(0, 1, 1 - 10^seq(-3, 0, length.out = 40))))
ridge_cv_kappa <- function(Y, X, Tstar, kgrid = KGRID, s_min = max(2, ceiling(Tstar / 2))) {
  folds <- s_min:(Tstar - 1)
  L <- sapply(folds, function(s) {
    pr <- ridge_prep(Y, X, s)
    sapply(kgrid, function(k) {
      cf <- ridge_coef(pr, k)
      mean(forecast_errors(Y, X, cf$a, cf$theta, s + 1, s + 1)^2)
    })
  })
  L <- matrix(L, nrow = length(kgrid))
  kgrid[which.min(rowMeans(L))]
}

# LASSO: glmnet with unpenalized intercept, lagY and X; penalized unit dummies;
# standardize = FALSE, so that a given lambda implies the same soft-threshold of
# the individual effects in every fold. Same folds and minimum rule as RIDGE.
lasso_cv <- function(Xm, y, times, pf, Tstar, s_min = max(2, ceiling(Tstar / 2))) {
  lam <- glmnet(Xm, y, alpha = 1, penalty.factor = pf, intercept = TRUE,
                standardize = FALSE)$lambda
  L <- sapply(s_min:(Tstar - 1), function(s) {
    tr <- times <= s; va <- times == s + 1
    fit <- glmnet(Xm[tr, , drop = FALSE], y[tr], alpha = 1, lambda = lam,
                  penalty.factor = pf, intercept = TRUE, standardize = FALSE)
    colMeans((y[va] - predict(fit, newx = Xm[va, , drop = FALSE]))^2)
  })
  L <- matrix(L, nrow = length(lam))
  lam[which.min(rowMeans(L))]
}

# =============================================================================
# 4. One design
# =============================================================================
run_design <- function(N = 100, T = 8, gamma = 0.8, phi = 1, s2 = 1, mu = 0,
                       sparse = 0, include_x = FALSE, beta = 1 - gamma,
                       rho_x = 0.5, sd_xi = 2, init = "stationary", burn = 50,
                       reps = 1000, seed = 123, skip_slow = FALSE) {
  Tstar <- T; T0 <- T + 1
  models <- c("OLS", "LSDV", "EBMLE", "URE", "AH", "GMM1", "GMM2", "RIDGE", "LASSO")
  mk <- function() matrix(NA_real_, reps, length(models), dimnames = list(NULL, models))
  g_hat <- mk(); o_bias <- mk(); o_var <- mk(); o_mse <- mk()
  kappa_cv <- rep(NA_real_, reps)

  set.seed(seed)
  for (it in 1:reps) {
    if (it %% 50 == 0) message("Iteration ", it, "/", reps)
    d <- simulate_panel(N, T0, gamma, mu, phi, s2, sparse, include_x, beta, rho_x, sd_xi, init, burn)
    Y <- d$Y; X <- d$X

    store <- function(name, a, theta) {
      if (is.null(a) || any(!is.finite(a)) || any(!is.finite(theta))) return(invisible())
      e <- forecast_errors(Y, X, a, theta, Tstar + 1, Tstar + 1)
      g_hat[it, name]  <<- theta[1]
      o_bias[it, name] <<- mean(e); o_var[it, name] <<- var(as.vector(e)); o_mse[it, name] <<- mean(e^2)
    }

    # Long data (training periods t = 1..Tstar) for lm, plm and glmnet
    tt <- 1:Tstar
    df <- data.frame(id = rep(1:N, times = Tstar), time = rep(tt, each = N),
                     Y = as.vector(Y[, tt + 1]), lagY = as.vector(Y[, tt]))
    if (include_x) df$X <- as.vector(X[, tt + 1])

    # --- Ridge path endpoints: OLS (kappa = 0), LSDV (kappa = 1)
    prep <- ridge_prep(Y, X, Tstar)
    cf <- ridge_coef(prep, 0); store("OLS", cf$a, cf$theta)
    cf <- ridge_coef(prep, 1); store("LSDV", cf$a, cf$theta)

    if (!skip_slow) {
    # --- Post-within shrinkage (EBMLE, URE)
    th_fe <- ridge_coef(prep, 1)$theta
    lv <- try(lsdv_intercepts_var(df, include_x), silent = TRUE)
    if (!inherits(lv, "try-error")) {
      store("EBMLE", shrink_fe(lv$a, lv$v, "EBMLE"), th_fe)
      store("URE",   shrink_fe(lv$a, lv$v, "URE"),   th_fe)
    }

    # --- AH in levels
    ah <- try(ah_fit(Y, X, Tstar), silent = TRUE)
    if (!inherits(ah, "try-error")) store("AH", levels_intercepts(Y, X, Tstar, ah$theta), ah$theta)

    # --- Arellano-Bond GMM, forecasts in levels
    pdf <- pdata.frame(df, index = c("id", "time"))
    f_gmm <- if (include_x) Y ~ lag(Y, 1) + X | lag(Y, 2:99) else Y ~ lag(Y, 1) | lag(Y, 2:99)
    for (step in c("onestep", "twosteps")) {
      gm <- try(pgmm(f_gmm, data = pdf, effect = "individual", model = step,
                     transformation = "d"), silent = TRUE)
      if (inherits(gm, "try-error")) next
      th <- c(coef(gm)["lag(Y, 1)"], if (include_x) coef(gm)["X"])
      store(if (step == "onestep") "GMM1" else "GMM2",
            levels_intercepts(Y, X, Tstar, th), as.numeric(th))
    }

    }  # end if (!skip_slow)

    # --- RIDGE: closed-form path, kappa chosen by corrected cross-validation
    kap <- try(ridge_cv_kappa(Y, X, Tstar), silent = TRUE)
    if (!inherits(kap, "try-error")) {
      cf <- ridge_coef(prep, kap); store("RIDGE", cf$a, cf$theta); kappa_cv[it] <- kap
    }

    # --- LASSO: glmnet, lambda chosen by corrected cross-validation
    D  <- model.matrix(~ factor(id) - 1, data = df)
    Xm <- cbind(lagY = df$lagY, if (include_x) cbind(X = df$X), D)
    pf <- c(0, if (include_x) 0, rep(1, N))
    lam <- try(lasso_cv(Xm, df$Y, df$time, pf, Tstar), silent = TRUE)
    if (!inherits(lam, "try-error")) {
      fit <- glmnet(Xm, df$Y, alpha = 1, lambda = lam, penalty.factor = pf,
                    intercept = TRUE, standardize = FALSE)
      b <- as.numeric(coef(fit))                     # (Intercept), lagY, [X], dummies
      nth <- if (include_x) 2 else 1
      store("LASSO", b[1] + b[(2 + nth):length(b)], b[2:(1 + nth)])
    }
  }

  summ <- function(M) apply(M, 2, mean, na.rm = TRUE)

  # Relative MSE with respect to LSDV and its Monte Carlo standard error.
  # Paired comparison: r = mean(A)/mean(B), se = sd(A - r*B) / (mean(B) * sqrt(R))
  rel_se <- sapply(models, function(m) {
    ok <- is.finite(o_mse[, m]) & is.finite(o_mse[, "LSDV"])
    if (sum(ok) < 2) return(NA_real_)
    a <- o_mse[ok, m]; b <- o_mse[ok, "LSDV"]; r <- mean(a) / mean(b)
    sd(a - r * b) / (mean(b) * sqrt(sum(ok)))
  })
  out <- rbind(Bias = summ(o_bias), Var = summ(o_var), MSE = summ(o_mse))
  gam <- rbind(Bias = summ(g_hat) - gamma,
               Var  = apply(g_hat, 2, var, na.rm = TRUE),
               MSE  = apply((g_hat - gamma)^2, 2, mean, na.rm = TRUE))
  theory <- if (!include_x && init %in% c("stationary", "uninformative"))
    theory_limits(gamma, Tstar, s2, (1 - sparse) * phi * s2, init) else NULL
  list(design = list(N = N, T = Tstar, gamma = gamma, phi = phi,
                     sparse = sparse, include_x = include_x, init = init),
       out_of_sample = out,
       relative_mse = out["MSE", ] / out["MSE", "LSDV"],
       relative_mse_se = rel_se,
       mse_by_rep = o_mse,
       gamma = gam,
       theory = theory,
       theory_relative = if (!is.null(theory)) theory / out["MSE", "LSDV"] else NULL,
       kappa_cv = mean(kappa_cv, na.rm = TRUE),
       failures = colSums(is.na(o_mse)))
}

# =============================================================================
# 5. Grid runner and table builder
# =============================================================================

# Runs every row of `grid` and saves one .rds file per design in `out_dir`.
# Designs already saved are skipped, so an interrupted run can be resumed.
run_grid <- function(grid, out_dir = "results/montecarlo", ...) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (r in seq_len(nrow(grid))) {
    g <- lapply(as.list(grid[r, , drop = FALSE]), function(v) if (is.factor(v)) as.character(v) else v)
    f <- file.path(out_dir, paste0(paste(names(g), sapply(g, as.character), sep = "",
                                         collapse = "_"), ".rds"))
    if (file.exists(f)) next
    message("Design ", r, "/", nrow(grid), ": ", basename(f))
    res <- do.call(run_design, c(g, list(...)))
    saveRDS(res, f)
  }
  invisible(out_dir)
}

# Collects the saved designs into a data frame of relative MSEs (vs simulated LSDV).
# Theoretical limits are added as rows (model = R_K, R_STAR, ...) with se = NA.
collect_results <- function(out_dir = "results/montecarlo") {
  files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    r <- readRDS(f); d <- r$design
    rel <- r$relative_mse; se <- r$relative_mse_se
    if (!is.null(r$theory_relative)) {
      rel <- c(rel, r$theory_relative); se <- c(se, rep(NA_real_, length(r$theory_relative)))
    }
    data.frame(N = d$N, T = d$T, gamma = d$gamma, phi = d$phi, init = d$init,
               model = names(rel), rel = as.numeric(rel), se = as.numeric(se),
               mse_lsdv = r$out_of_sample["MSE", "LSDV"])
  }))
}

# LaTeX table: rows (N, T, gamma), columns = models, entries rel. MSE (s.e.)
latex_table <- function(df, models, phi_value = 1, init_value = "stationary", digits = 3,
                        caption = "Out-of-sample MSE relative to fixed effects",
                        label = "tab:mc") {
  df <- df[df$phi == phi_value & df$init == init_value & df$model %in% models, ]
  keys <- unique(df[order(df$N, df$T, df$gamma), c("N", "T", "gamma")])
  fmt <- function(x) formatC(x, format = "f", digits = digits)
  lines <- c("\\begin{table}[ht]", "\\centering", "\\small",
             paste0("\\caption{", caption, "}\\label{", label, "}"),
             paste0("\\begin{tabular}{ccc", strrep("c", length(models)), "}"),
             "\\toprule",
             paste0("$N$ & $T$ & $\\gamma$ & ", paste(gsub("_", "\\\\_", models), collapse = " & "), " \\\\"),
             "\\midrule")
  for (k in seq_len(nrow(keys))) {
    sub <- df[df$N == keys$N[k] & df$T == keys$T[k] & df$gamma == keys$gamma[k], ]
    v  <- sub$rel[match(models, sub$model)]; se <- sub$se[match(models, sub$model)]
    lines <- c(lines,
      paste0(keys$N[k], " & ", keys$T[k], " & ", keys$gamma[k], " & ",
             paste(fmt(v), collapse = " & "), " \\\\"),
      paste0(" & & & ", paste(ifelse(is.na(se), "", paste0("(", fmt(se), ")")), collapse = " & "),
             " \\\\"))
  }
  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

# =============================================================================
# 5b. Risk along the degree of shrinkage (Figure 2 of the paper)
#     For a grid of kappa (1 - kappa = degree of shrinkage), computes the simulated
#     MSE of three forecasts and their theoretical limits:
#       ridge:       kappa with the ridge estimator of gamma at omega = T(1 - kappa)
#       post-within: kappa with the within estimator of gamma
#       known gamma: kappa with the true gamma
#     Values are mean squared forecast errors (not relative to fixed effects).
# =============================================================================
shrinkage_path <- function(N = 100, T = 5, gamma = 0.5, phi = 1, s2 = 1, reps = 1000,
                           kappa_grid = seq(0, 1, by = 0.05), seed = 2026) {
  set.seed(seed)
  K <- length(kappa_grid)
  sim <- matrix(0, K, 3, dimnames = list(NULL, c("ridge", "within", "known")))
  for (it in 1:reps) {
    Y <- simulate_panel(N, T + 1, gamma, 0, phi, s2, 0, FALSE, 0, 0, 0, "stationary", 0)$Y
    prep <- ridge_prep(Y, NULL, T)
    g_w  <- prep$SWwy / prep$SWww[1, 1]
    ybar_i <- rowMeans(Y[, 2:(T + 1)]); xbar_i <- rowMeans(Y[, 1:T])
    for (j in 1:K) {
      k <- kappa_grid[j]
      g_r <- ridge_coef(prep, k)$theta[1]
      q <- (Y[, T + 1] - mean(xbar_i)) - k * (xbar_i - mean(xbar_i))
      base <- mean(ybar_i) + k * (ybar_i - mean(ybar_i))
      tgt <- Y[, T + 2]
      sim[j, ] <- sim[j, ] + c(mean((tgt - base - g_r * q)^2),
                               mean((tgt - base - g_w * q)^2),
                               mean((tgt - base - gamma * q)^2))
    }
  }
  sim <- sim / reps
  sa2 <- phi * s2
  lim <- t(sapply(kappa_grid, function(k)
    c(ridge  = risk_T(k, gamma + bias_T(T * (1 - k), gamma, T, s2, sa2), gamma, T, s2, sa2),
      within = risk_T(k, gamma + bias_T(0, gamma, T, s2, sa2), gamma, T, s2, sa2),
      known  = risk_T(k, gamma, gamma, T, s2, sa2))))
  th <- theory_limits(gamma, T, s2, sa2, "stationary")
  data.frame(T = T, gamma = gamma, shrinkage = 1 - kappa_grid,
             sim_ridge = sim[, "ridge"], sim_within = sim[, "within"], sim_known = sim[, "known"],
             lim_ridge = lim[, "ridge"], lim_within = lim[, "within"], lim_known = lim[, "known"],
             R_STAR = th[["R_STAR"]], R_FE = th[["R_FE"]])
}

# Example (uncomment): the six panels of Figure 2
# path <- do.call(rbind, lapply(c(5, 10), function(Tw) do.call(rbind, lapply(c(0.2, 0.5, 0.8),
#           function(g) shrinkage_path(T = Tw, gamma = g)))))
# write.csv(path, "mc_shrinkage_path.csv", row.names = FALSE)
# library(ggplot2)
# cap <- function(v, top) ifelse(v > top, NA, v)          # truncate at 1.3 x limiting risk of FE
# for (v in grep("^(sim|lim)_", names(path), value = TRUE)) path[[v]] <- cap(path[[v]], 1.3 * path$R_FE)
# ggplot(path, aes(shrinkage)) +
#   geom_line(aes(y = lim_ridge)) + geom_point(aes(y = sim_ridge), size = 1) +
#   geom_line(aes(y = lim_within), linetype = "dashed", colour = "grey35") +
#   geom_point(aes(y = sim_within), shape = 0, size = 1, colour = "grey35") +
#   geom_line(aes(y = lim_known), linetype = "dotted") + geom_point(aes(y = sim_known), shape = 2, size = 1) +
#   geom_hline(aes(yintercept = R_STAR), linetype = "dotdash", colour = "grey50") +
#   facet_wrap(~ T + gamma, scales = "free_y", labeller = label_both) +
#   labs(x = expression(1 - kappa), y = "Mean squared forecast error") + theme_bw()

# =============================================================================
# 6. Run
# =============================================================================
# Quick check with one design:
# res <- run_design(N = 100, T = 8, gamma = 0.5, phi = 1, reps = 1000)
# print(round(res$relative_mse, 3)); print(round(res$theory_relative, 3))

# Designs of Table 1 and Tables 4 and 5 in the appendix.
grid1 <- expand.grid(N = c(20, 100), T = c(5, 10, 20, 30), gamma = c(0.2, 0.5, 0.8),
                     phi = 1, init = "stationary", stringsAsFactors = FALSE)
grid2 <- expand.grid(N = 100, T = c(5, 10, 20, 30), gamma = c(0.2, 0.5, 0.8),
                     phi = c(0.25, 4), init = "stationary", stringsAsFactors = FALSE)
run_grid(rbind(grid1, grid2), reps = 1000)
# skip_slow = TRUE computes only OLS, LSDV, RIDGE and LASSO (much faster).

# Summary. The tables in the paper report the mean squared forecast error,
# mse = rel * mse_lsdv (rel is the MSE relative to fixed effects).
df <- collect_results("results/montecarlo")
df$mse <- df$rel * df$mse_lsdv
write.csv(df, "results/montecarlo/mc_results.csv", row.names = FALSE)

# Figures 1 and 2 are produced by code/figures.R.
