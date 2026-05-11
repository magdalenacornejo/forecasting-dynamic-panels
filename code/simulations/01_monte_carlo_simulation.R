# ============================================================
# Monte Carlo simulations
# Paper: Forecasting Dynamic Panel Models with Shrinkage of Fixed Effects
# Authors: Magdalena Cornejo and Walter Sosa-Escudero
#
# This script runs one Monte Carlo design at a time. To replicate
# all designs in the paper, change the parameters in Section 1
# or call this script from a batch wrapper.
# ============================================================

rm(list = ls())
options(scipen = 999)

# -----------------------------
# 0) Packages
# -----------------------------
# install.packages("remotes")
# https://github.com/soonwookwon/FEShR #download zip file from here

# zipfile <- "C:/Users/mcornejo/Downloads/FEShR-master.zip"
# outdir  <- "C:/Users/mcornejo/Downloads/FEShR-master"
# unzip(zipfile, exdir = outdir)
# list.files(outdir)
# pkgdir <- file.path(outdir, "FEShR-master")
# 
# files <- c(
#   file.path(pkgdir, "src", "Makevars"),
#   file.path(pkgdir, "src", "Makevars.win")
# )
# 
# for (f in files) {
#   if (file.exists(f)) {
#     txt <- readLines(f)
#     txt <- gsub("CXX_STD[[:space:]]*=[[:space:]]*CXX11", "CXX_STD = CXX14", txt)
#     writeLines(txt, f)
#   }
# }
# 
# remotes::install_local(pkgdir, force = TRUE)

library(plm)
library(dplyr)
library(glmnet)
library(fastDummies)
library(FEShR)

# -----------------------------
# 1) Simulation parameters
# -----------------------------
sparse <- 0.1 # fraction of fixed effects set equal to zero
num_iterations <- 1000

N  <- 20 # 20 or 100
T0 <- 5  # 5, 10, 20 or 30
T  <- T0 + 10

rho   <- 0.5
gamma <- 0.2 # 0.2 or 0.8
beta  <- 1 - gamma

Ttr <- floor(T0 * 0.8)   # train cutoff in "time" after dropping first 10 periods

# -----------------------------
# 2) Helper functions
# -----------------------------
safe_colmean <- function(M) {
  # column means with NA handling (works for vectors/matrices)
  if (is.null(dim(M))) return(mean(M, na.rm = TRUE))
  apply(M, 2, function(x) mean(x, na.rm = TRUE))
}

make_year_blocks <- function(year, K = 5) {
  yrs <- sort(unique(year))
  split(yrs, cut(seq_along(yrs), breaks = K, labels = FALSE))
}

make_lambda_seq <- function(X, y, alpha, penalty_factor) {
  fit0 <- glmnet(X, y, alpha = alpha, penalty.factor = penalty_factor,
                 intercept = FALSE, standardize = TRUE)
  fit0$lambda
}

rolling_origin_select_lambda <- function(X, y, years, year_blocks, alpha,
                                         penalty_factor, lambda_seq, use_1se = TRUE) {
  K <- length(year_blocks)
  val_blocks <- 2:K
  L <- matrix(NA_real_, nrow = length(lambda_seq), ncol = length(val_blocks))
  
  for (jj in seq_along(val_blocks)) {
    k <- val_blocks[jj]
    tr <- years %in% unlist(year_blocks[1:(k - 1)])
    va <- years %in% unlist(year_blocks[k])
    
    fit <- glmnet(X[tr, , drop = FALSE], y[tr],
                  alpha = alpha, lambda = lambda_seq,
                  penalty.factor = penalty_factor,
                  intercept = FALSE, standardize = TRUE)
    
    pred <- predict(fit, newx = X[va, , drop = FALSE], s = lambda_seq)
    L[, jj] <- colMeans((y[va] - pred)^2)
  }
  
  mean_loss <- rowMeans(L, na.rm = TRUE)
  se_loss   <- apply(L, 1, sd, na.rm = TRUE) / sqrt(ncol(L))
  j_min     <- which.min(mean_loss)
  
  if (!use_1se) return(lambda_seq[j_min])
  
  thresh <- mean_loss[j_min] + se_loss[j_min]
  max(lambda_seq[mean_loss <= thresh])
}

eval_metrics <- function(y, yhat) {
  y    <- as.numeric(y)
  yhat <- as.numeric(yhat)
  ok <- is.finite(y) & is.finite(yhat)
  if (sum(ok) == 0) return(c(bias = NA_real_, var = NA_real_, mse = NA_real_))
  
  err <- y[ok] - yhat[ok]
  c(
    bias = mean(err),
    var  = var(yhat[ok]),   # variance of predictions (your definition)
    mse  = mean(err^2)
  )
}

compute_shrunk_lsdv <- function(in.sample.df, shrink_type = c("EBMLE", "URE"),
                                centering = c("gen", "0")) {
  shrink_type <- match.arg(shrink_type)
  centering   <- match.arg(centering)
  
  in.sample.df$id <- as.factor(in.sample.df$id)
  
  fit <- lm(Y ~ lagY + X + id, data = in.sample.df)
  
  ids   <- levels(in.sample.df$id)
  coefs <- coef(fit)
  V     <- vcov(fit)
  
  # FE estimados
  alpha_hat <- numeric(length(ids))
  names(alpha_hat) <- ids
  
  alpha_hat[ids[1]] <- coefs["(Intercept)"]
  
  if (length(ids) >= 2) {
    for (j in 2:length(ids)) {
      cname <- paste0("id", ids[j])
      alpha_hat[ids[j]] <- coefs["(Intercept)"] + coefs[cname]
    }
  }
  
  # Varianzas de los FE
  M_list <- vector("list", length(ids))
  names(M_list) <- ids
  
  M_list[[ids[1]]] <- matrix(V["(Intercept)", "(Intercept)"], 1, 1)
  
  if (length(ids) >= 2) {
    for (j in 2:length(ids)) {
      cname <- paste0("id", ids[j])
      var_alpha_j <- V["(Intercept)", "(Intercept)"] +
        V[cname, cname] +
        2 * V["(Intercept)", cname]
      M_list[[ids[j]]] <- matrix(var_alpha_j, 1, 1)
    }
  }
  
  y_mat <- matrix(alpha_hat, nrow = 1)
  
  shrink_fit <- FEShR::fe_shrink(
    y = y_mat,
    M = M_list,
    centering = centering,
    type = shrink_type,
    n_init_vals = 1,
    all_init_vals = FALSE,
    optim_control = list(maxit = 500)
  )
  
  alpha_shrunk <- as.numeric(shrink_fit$thetahat)
  names(alpha_shrunk) <- ids
  
  list(
    fit_lm       = fit,
    alpha_fe     = alpha_hat,
    alpha_shrunk = alpha_shrunk,
    shrink_fit   = shrink_fit
  )
}

# ------------------------------------------------------------
# Robust GMM in-sample reconstruction in levels for pgmm "d"
# ------------------------------------------------------------
predict_gmm_levels_1step <- function(gmm_obj, df) {
  g <- as.numeric(coef(gmm_obj)["lag(Y, 1)"])
  b <- as.numeric(coef(gmm_obj)["X"])
  
  # ΔYhat_t = g*ΔY_{t-1} + b*ΔX_t
  dYhat <- g * as.numeric(df$LDY) + b * as.numeric(df$DX)
  
  # Yhat_t = Y_{t-1} + ΔYhat_t
  yhat <- as.numeric(df$lagY) + dYhat
  
  yhat
}

# -----------------------------
# 3) Base panel: X and fixed effects are fixed across replications
# -----------------------------
panel <- data.frame(
  id   = rep(1:N, each = T),
  time = rep(1:T, times = N),
  X    = NA_real_
)

set.seed(123)

panel$X[panel$time == 1] <- 0
for (i in 1:N) {
  for (t in 2:T) {
    panel$X[(panel$id == i) & (panel$time == t)] <-
      rho * panel$X[(panel$id == i) & (panel$time == t - 1)] +
      rnorm(1, mean = 0, sd = 2)
  }
}

panel <- panel %>% arrange(id, time)
panel$FE <- rep(rnorm(N, mean = 0, sd = 1), each = T)

num_of_zeros <- round(sparse * N)
if (num_of_zeros > 0) {
  indices_to_replace <- sample(N, num_of_zeros)
  panel$FE[panel$id %in% indices_to_replace] <- 0
}

# -----------------------------
# 4) Storage: scalar metrics by iteration and estimator
# -----------------------------
models <- c("OLS", "LSDV", "EBMLE", "URE", "AH", "GMM1", "GMM2",
            "LASSO", "RIDGE", "ELASTIC_NET")

gamma_hat <- matrix(NA_real_, nrow = num_iterations, ncol = length(models),
                    dimnames = list(NULL, models))

in_bias  <- in_var  <- in_mse  <- matrix(NA_real_, nrow = num_iterations, ncol = length(models),
                                         dimnames = list(NULL, models))
out_bias <- out_var <- out_mse <- matrix(NA_real_, nrow = num_iterations, ncol = length(models),
                                         dimnames = list(NULL, models))

# -----------------------------
# 5) Diagnostics counters
# -----------------------------
gmm1_fail <- 0
gmm2_fail <- 0
gmm1_empty_align <- 0
gmm2_empty_align <- 0
gmm1_align_error <- 0
gmm2_align_error <- 0

# -----------------------------
# 6) Monte Carlo loop
# -----------------------------
for (it in 1:num_iterations) {
  set.seed(123 + it)
  if (it %% 50 == 0) message("Iteration: ", it, "/", num_iterations)
  
  # --- Simulate Y ---
  panel_new <- panel
  panel_new$Y <- NA_real_
  panel_new$Y[panel_new$time == 1] <- 0
  panel_new$e <- rnorm(N * T, mean = 0, sd = 1)
  
  for (j in 1:N) {
    for (t in 2:T) {
      idx <- (panel_new$id == j) & (panel_new$time == t)
      idx_lag <- (panel_new$id == j) & (panel_new$time == t - 1)
      
      panel_new$Y[idx] <-
        gamma * panel_new$Y[idx_lag] +
        beta  * panel_new$X[idx] +
        panel_new$FE[idx] +
        panel_new$e[idx]
    }
  }
  
  # --- Transformations ---
  panel_new <- panel_new %>%
    group_by(id) %>%
    mutate(
      lagY = dplyr::lag(Y, 1),
      DY   = Y - dplyr::lag(Y, 1),
      L2Y  = dplyr::lag(lagY, 1),
      LDY  = dplyr::lag(DY, 1),
      lagX = dplyr::lag(X, 1),
      DX   = X - lagX
    ) %>%
    ungroup()
  
  # Drop first 10 periods (burn-in)
  panel_new <- as.data.frame(panel_new) %>% filter(time > 10)
  panel_new$time <- panel_new$time - 10
  
  # FE dummies (used in regularized methods only)
  panel_new <- fastDummies::dummy_cols(panel_new, select_columns = "id",
                                       remove_first_dummy = FALSE,
                                       remove_selected_columns = FALSE)
  
  # Split train/test
  in.sample     <- panel_new %>% filter(time <= Ttr)
  out.of.sample <- panel_new %>% filter(time >  Ttr)
  
  in.sample.df  <- panel_new %>% filter(time <= Ttr)
  out.sample.df <- panel_new %>% filter(time >  Ttr)
  
  in.sample     <- pdata.frame(in.sample.df, index = c("id", "time"))
  out.of.sample <- pdata.frame(out.sample.df, index = c("id", "time"))
  
  y_in  <- as.numeric(in.sample$Y)
  y_out <- as.numeric(out.of.sample$Y)
  
  X_in_df  <- in.sample %>% dplyr::select(lagY, X, starts_with("id_"))
  X_out_df <- out.of.sample %>% dplyr::select(lagY, X, starts_with("id_"))
  
  X_in  <- as.matrix(X_in_df)
  X_out <- as.matrix(X_out_df)
  
  penalty_factors <- c(0, 0, rep(1, N))
  
  # -----------------------------
  # Estimation
  # -----------------------------
  OLS <- LSDV <- EBMLE <- URE <- AH <- GMM1 <- GMM2 <- NULL
  
  try({
    OLS   <- plm(Y ~ lagY + X, data = in.sample, model = "pooling")
    LSDV  <- plm(Y ~ lagY + X, data = in.sample, model = "within", effect = "individual")
    EBMLE <- compute_shrunk_lsdv(in.sample.df, shrink_type = "EBMLE", centering = "gen")
    URE   <- compute_shrunk_lsdv(in.sample.df, shrink_type = "URE",   centering = "gen")
    AH    <- plm(DY ~ LDY - 1 + DX | L2Y + DX, data = in.sample, model = "pooling")
  }, silent = TRUE)
  
  try({
    GMM1 <- pgmm(
      Y ~ lag(Y, 1) + X | lag(Y, 2:(Ttr - 2)),
      data = in.sample,
      effect = "individual",
      model = "onestep",
      transformation = "d",
      collapse = FALSE
    )
  }, silent = TRUE)
  
  try({
    GMM2 <- pgmm(
      Y ~ lag(Y, 1) + X | lag(Y, 2:(Ttr - 2)),
      data = in.sample,
      effect = "individual",
      model = "twosteps",
      transformation = "d",
      collapse = FALSE
    )
  }, silent = TRUE)
  
  if (is.null(GMM1)) gmm1_fail <- gmm1_fail + 1
  if (is.null(GMM2)) gmm2_fail <- gmm2_fail + 1
  
  # -----------------------------
  # Lambda selection (rolling origin CV)
  # -----------------------------
  years_vec <- as.numeric(attr(in.sample, "index")$time)
  blocks    <- make_year_blocks(years_vec, K = 5)
  p_factor <- penalty_factors
  
  l_seq_lasso       <- make_lambda_seq(X_in, y_in, alpha = 1, penalty_factor = p_factor)
  best_lambda_lasso <- rolling_origin_select_lambda(X_in, y_in, years_vec, blocks,
                                                    alpha = 1, penalty_factor = p_factor,
                                                    lambda_seq = l_seq_lasso, use_1se = TRUE)
  
  l_seq_ridge       <- make_lambda_seq(X_in, y_in, alpha = 0, penalty_factor = p_factor)
  best_lambda_ridge <- rolling_origin_select_lambda(X_in, y_in, years_vec, blocks,
                                                    alpha = 0, penalty_factor = p_factor,
                                                    lambda_seq = l_seq_ridge, use_1se = TRUE)
  
  alpha_en          <- 0.5
  l_seq_en          <- make_lambda_seq(X_in, y_in, alpha = alpha_en, penalty_factor = p_factor)
  best_lambda_en    <- rolling_origin_select_lambda(X_in, y_in, years_vec, blocks,
                                                    alpha = alpha_en, penalty_factor = p_factor,
                                                    lambda_seq = l_seq_en, use_1se = TRUE)
  
  LASSO <- glmnet(X_in, y_in, alpha = 1, lambda = best_lambda_lasso,
                  penalty.factor = penalty_factors, intercept = FALSE)
  RIDGE <- glmnet(X_in, y_in, alpha = 0, lambda = best_lambda_ridge,
                  penalty.factor = penalty_factors, intercept = FALSE)
  ELASTIC_NET <- glmnet(X_in, y_in, alpha = alpha_en, lambda = best_lambda_en,
                        penalty.factor = penalty_factors, intercept = FALSE)
  
  # -----------------------------
  # Store gamma estimates
  # -----------------------------
  if (!is.null(OLS))  gamma_hat[it, "OLS"]  <- coef(OLS)["lagY"]
  if (!is.null(LSDV)) gamma_hat[it, "LSDV"] <- coef(LSDV)["lagY"]
  if (!is.null(AH))   gamma_hat[it, "AH"]   <- coef(AH)["LDY"]
  if (!is.null(GMM1)) gamma_hat[it, "GMM1"] <- coef(GMM1)["lag(Y, 1)"]
  if (!is.null(GMM2)) gamma_hat[it, "GMM2"] <- coef(GMM2)["lag(Y, 1)"]
  if (!is.null(EBMLE)) gamma_hat[it, "EBMLE"] <- coef(EBMLE$fit_lm)["lagY"]
  if (!is.null(URE)) gamma_hat[it, "URE"] <- coef(URE$fit_lm)["lagY"]
  gamma_hat[it, "LASSO"]       <- as.numeric(coef(LASSO)["lagY", ])
  gamma_hat[it, "RIDGE"]       <- as.numeric(coef(RIDGE)["lagY", ])
  gamma_hat[it, "ELASTIC_NET"] <- as.numeric(coef(ELASTIC_NET)["lagY", ])
  
  # -----------------------------
  # In-sample predictions + scalar evaluation
  # -----------------------------
  if (!is.null(OLS)) {
    yhat <- as.numeric(coef(OLS)["(Intercept)"] + coef(OLS)["lagY"] * in.sample$lagY + coef(OLS)["X"] * in.sample$X)
    m <- eval_metrics(y_in, yhat)
    in_bias[it, "OLS"] <- m["bias"]; in_var[it, "OLS"] <- m["var"]; in_mse[it, "OLS"] <- m["mse"]
  }
  
  if (!is.null(LSDV)) {
    fe_i <- fixef(LSDV, effect = "individual")
    yhat <- as.numeric(fe_i[as.character(in.sample$id)] + coef(LSDV)["lagY"] * in.sample$lagY + coef(LSDV)["X"] * in.sample$X)
    m <- eval_metrics(y_in, yhat)
    in_bias[it, "LSDV"] <- m["bias"]; in_var[it, "LSDV"] <- m["var"]; in_mse[it, "LSDV"] <- m["mse"]
  }
  
  if (!is.null(EBMLE)) {
    beta_lagY <- coef(EBMLE$fit_lm)["lagY"]
    beta_X    <- coef(EBMLE$fit_lm)["X"]
    
    ids_in <- as.character(in.sample.df$id)
    
    yhat <- as.numeric(
      EBMLE$alpha_shrunk[ids_in] +
        beta_lagY * in.sample.df$lagY +
        beta_X    * in.sample.df$X
    )
    
    m <- eval_metrics(y_in, yhat)
    in_bias[it, "EBMLE"] <- m["bias"]
    in_var[it, "EBMLE"]  <- m["var"]
    in_mse[it, "EBMLE"]  <- m["mse"]
  }
  
  if (!is.null(URE)) {
    beta_lagY <- coef(URE$fit_lm)["lagY"]
    beta_X    <- coef(URE$fit_lm)["X"]
    
    ids_in <- as.character(in.sample.df$id)
    
    yhat <- as.numeric(
      URE$alpha_shrunk[ids_in] +
        beta_lagY * in.sample.df$lagY +
        beta_X    * in.sample.df$X
    )
    
    m <- eval_metrics(y_in, yhat)
    in_bias[it, "URE"] <- m["bias"]
    in_var[it, "URE"]  <- m["var"]
    in_mse[it, "URE"]  <- m["mse"]
  }
  
  if (!is.null(AH)) {
    yhat <- as.numeric(fitted(AH)) + as.numeric(in.sample$lagY)
    m <- eval_metrics(y_in, yhat)
    in_bias[it, "AH"] <- m["bias"]; in_var[it, "AH"] <- m["var"]; in_mse[it, "AH"] <- m["mse"]
  }
  
  if (!is.null(GMM1)) {
    yhat <- predict_gmm_levels_1step(GMM1, in.sample)
    m <- eval_metrics(y_in, yhat)
    in_bias[it, "GMM1"] <- m["bias"]; in_var[it, "GMM1"] <- m["var"]; in_mse[it, "GMM1"] <- m["mse"]
  }
  
  if (!is.null(GMM2)) {
    yhat <- predict_gmm_levels_1step(GMM2, in.sample)
    m <- eval_metrics(y_in, yhat)
    in_bias[it, "GMM2"] <- m["bias"]; in_var[it, "GMM2"] <- m["var"]; in_mse[it, "GMM2"] <- m["mse"]
  }
  
  yhat <- as.numeric(predict(LASSO, newx = X_in, s = best_lambda_lasso))
  m <- eval_metrics(y_in, yhat)
  in_bias[it, "LASSO"] <- m["bias"]; in_var[it, "LASSO"] <- m["var"]; in_mse[it, "LASSO"] <- m["mse"]
  
  yhat <- as.numeric(predict(RIDGE, newx = X_in, s = best_lambda_ridge))
  m <- eval_metrics(y_in, yhat)
  in_bias[it, "RIDGE"] <- m["bias"]; in_var[it, "RIDGE"] <- m["var"]; in_mse[it, "RIDGE"] <- m["mse"]
  
  yhat <- as.numeric(predict(ELASTIC_NET, newx = X_in, s = best_lambda_en))
  m <- eval_metrics(y_in, yhat)
  in_bias[it, "ELASTIC_NET"] <- m["bias"]; in_var[it, "ELASTIC_NET"] <- m["var"]; in_mse[it, "ELASTIC_NET"] <- m["mse"]
  
  # -----------------------------
  # Out-of-sample predictions + scalar evaluation
  # -----------------------------
  if (!is.null(OLS)) {
    yhat <- as.numeric(coef(OLS)["(Intercept)"] + coef(OLS)["lagY"] * out.of.sample$lagY + coef(OLS)["X"] * out.of.sample$X)
    m <- eval_metrics(y_out, yhat)
    out_bias[it, "OLS"] <- m["bias"]; out_var[it, "OLS"] <- m["var"]; out_mse[it, "OLS"] <- m["mse"]
  }
  
  if (!is.null(LSDV)) {
    fe_i <- fixef(LSDV, effect = "individual")
    yhat <- as.numeric(fe_i[as.character(out.of.sample$id)] + coef(LSDV)["lagY"] * out.of.sample$lagY + coef(LSDV)["X"] * out.of.sample$X)
    m <- eval_metrics(y_out, yhat)
    out_bias[it, "LSDV"] <- m["bias"]; out_var[it, "LSDV"] <- m["var"]; out_mse[it, "LSDV"] <- m["mse"]
  }
  
  if (!is.null(EBMLE)) {
    beta_lagY <- coef(EBMLE$fit_lm)["lagY"]
    beta_X    <- coef(EBMLE$fit_lm)["X"]
    
    ids_out <- as.character(out.sample.df$id)
    
    yhat <- as.numeric(
      EBMLE$alpha_shrunk[ids_out] +
        beta_lagY * out.sample.df$lagY +
        beta_X    * out.sample.df$X
    )
    
    m <- eval_metrics(y_out, yhat)
    out_bias[it, "EBMLE"] <- m["bias"]
    out_var[it, "EBMLE"]  <- m["var"]
    out_mse[it, "EBMLE"]  <- m["mse"]
  }
  
  if (!is.null(URE)) {
    beta_lagY <- coef(URE$fit_lm)["lagY"]
    beta_X    <- coef(URE$fit_lm)["X"]
    
    ids_out <- as.character(out.sample.df$id)
    
    yhat <- as.numeric(
      URE$alpha_shrunk[ids_out] +
        beta_lagY * out.sample.df$lagY +
        beta_X    * out.sample.df$X
    )
    
    m <- eval_metrics(y_out, yhat)
    out_bias[it, "URE"] <- m["bias"]
    out_var[it, "URE"]  <- m["var"]
    out_mse[it, "URE"]  <- m["mse"]
  }
  
  if (!is.null(AH)) {
    yhat <- as.numeric(coef(AH)["LDY"] * out.of.sample$LDY + coef(AH)["DX"] * out.of.sample$DX + out.of.sample$lagY)
    m <- eval_metrics(y_out, yhat)
    out_bias[it, "AH"] <- m["bias"]; out_var[it, "AH"] <- m["var"]; out_mse[it, "AH"] <- m["mse"]
  }
  
  if (!is.null(GMM1)) {
    yhat <- predict_gmm_levels_1step(GMM1, out.of.sample)
    m <- eval_metrics(y_out, yhat)
    out_bias[it, "GMM1"] <- m["bias"]; out_var[it, "GMM1"] <- m["var"]; out_mse[it, "GMM1"] <- m["mse"]
  }
  
  if (!is.null(GMM2)) {
    yhat <- predict_gmm_levels_1step(GMM2, out.of.sample)
    m <- eval_metrics(y_out, yhat)
    out_bias[it, "GMM2"] <- m["bias"]; out_var[it, "GMM2"] <- m["var"]; out_mse[it, "GMM2"] <- m["mse"]
  }
  
  yhat <- as.numeric(predict(LASSO, newx = X_out, s = best_lambda_lasso))
  m <- eval_metrics(y_out, yhat)
  out_bias[it, "LASSO"] <- m["bias"]; out_var[it, "LASSO"] <- m["var"]; out_mse[it, "LASSO"] <- m["mse"]
  
  yhat <- as.numeric(predict(RIDGE, newx = X_out, s = best_lambda_ridge))
  m <- eval_metrics(y_out, yhat)
  out_bias[it, "RIDGE"] <- m["bias"]; out_var[it, "RIDGE"] <- m["var"]; out_mse[it, "RIDGE"] <- m["mse"]
  
  yhat <- as.numeric(predict(ELASTIC_NET, newx = X_out, s = best_lambda_en))
  m <- eval_metrics(y_out, yhat)
  out_bias[it, "ELASTIC_NET"] <- m["bias"]; out_var[it, "ELASTIC_NET"] <- m["var"]; out_mse[it, "ELASTIC_NET"] <- m["mse"]
}

# -----------------------------
# 7) Diagnostics report
# -----------------------------
cat("\n--- GMM diagnostics ---\n")
cat("GMM1 NULL count:        ", gmm1_fail, "out of", num_iterations, "\n")
cat("GMM2 NULL count:        ", gmm2_fail, "out of", num_iterations, "\n")
cat("GMM1 empty alignment:   ", gmm1_empty_align, "out of", num_iterations, "\n")
cat("GMM2 empty alignment:   ", gmm2_empty_align, "out of", num_iterations, "\n")
cat("GMM1 alignment errors:  ", gmm1_align_error, "out of", num_iterations, "\n")
cat("GMM2 alignment errors:  ", gmm2_align_error, "out of", num_iterations, "\n\n")

# -----------------------------
# 8) Table 1: gamma estimation
# -----------------------------
gamma_bias <- gamma - safe_colmean(gamma_hat)
gamma_var  <- apply(gamma_hat, 2, function(x) mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
gamma_mse  <- apply(gamma_hat, 2, function(x) mean((x - gamma)^2, na.rm = TRUE))

table1 <- rbind(Bias = gamma_bias, Var = gamma_var, MSE = gamma_mse)
table1 <- as.data.frame(round(table1, 5))
table1 <- cbind(Measure = rownames(table1), table1)
rownames(table1) <- NULL
print(table1)

# -----------------------------
# 9) Table 2: in-sample prediction evaluation
# -----------------------------
table2 <- rbind(
  Bias = safe_colmean(in_bias),
  Var  = safe_colmean(in_var),
  MSE  = safe_colmean(in_mse)
)
table2 <- as.data.frame(round(table2, 5))
table2 <- cbind(Measure = rownames(table2), table2)
rownames(table2) <- NULL
print(table2)

# -----------------------------
# 10) Table 3: out-of-sample forecast evaluation
# -----------------------------
table3 <- rbind(
  Bias = safe_colmean(out_bias),
  Var  = safe_colmean(out_var),
  MSE  = safe_colmean(out_mse)
)
table3 <- as.data.frame(round(table3, 5))
table3 <- cbind(Measure = rownames(table3), table3)
rownames(table3) <- NULL
print(table3)


# Optional: save the three tables for the current design.
# Uncomment these lines if you want one CSV per scenario.
# output_dir <- file.path("output", "tables")
# dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# scenario_id <- paste0("N", N, "_T", T0, "_g", gamma, "_s", sparse)
# write.csv(table1, file.path(output_dir, paste0("table1_gamma_", scenario_id, ".csv")), row.names = FALSE)
# write.csv(table2, file.path(output_dir, paste0("table2_insample_", scenario_id, ".csv")), row.names = FALSE)
# write.csv(table3, file.path(output_dir, paste0("table3_oos_", scenario_id, ".csv")), row.names = FALSE)
