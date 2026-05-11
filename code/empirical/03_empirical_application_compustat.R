# ============================================================
# Empirical application: Compustat / WRDS
# Paper: Forecasting Dynamic Panel Models with Shrinkage of Fixed Effects
#
# Single split matched to the Monte Carlo exercise
# Split:
#   IS  = 1991–2010 (T=20)
#   OOS = 2011–2024 (H=14)
#
# Matching Monte Carlo:
#  - metrics: err=y-yhat; bias=mean(err); var=var(yhat); mse=mean(err^2)
#  - AH in levels: yhat = y_{t-1} + Δyhat
#  - GMM in levels (MC-style): yhat = y_{t-1} + (g*Δy_{t-1} + b'ΔX_t)
#  - Shrinkage: glmnet with intercept=FALSE; penalize only firm dummies;
#               lambda via rolling-origin blocks (1SE rule)
# ============================================================

rm(list = ls())
options(scipen = 999)

# -----------------------------
# Paths
# -----------------------------
# Raw Compustat data are proprietary and should not be committed to GitHub.
# Place the WRDS/Compustat extract locally at data/raw/df0_WRDS.csv.
input_file <- file.path("data", "raw", "df0_WRDS.csv")
table_dir  <- file.path("output", "tables", "compustat")
figure_dir <- file.path("output", "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(plm)
library(Matrix)
library(glmnet)
library(FEShR)
library(kableExtra)
library(tidyr)
library(ggplot2)

# -----------------------------
# 1) Settings: single split
# -----------------------------
IS_START  <- 1991
IS_END    <- 2010
OOS_START <- 2011
OOS_END   <- 2024

K_CV   <- 5
USE_1SE <- TRUE
GMM_COLLAPSE <- FALSE  # try TRUE if instruments explode

# -----------------------------
# 2) Helper functions: MC-style metrics and shrinkage utilities
# -----------------------------
eval_metrics_mc <- function(y, yhat) {
  y    <- as.numeric(y)
  yhat <- as.numeric(yhat)
  
  ok <- is.finite(y) & is.finite(yhat)
  
  if (sum(ok) == 0) {
    return(c(
      bias = NA_real_,
      var  = NA_real_,
      mse  = NA_real_,
      mae  = NA_real_
    ))
  }
  
  err <- y[ok] - yhat[ok]
  
  c(
    bias = mean(err),
    var  = var(yhat[ok]),
    mse  = mean(err^2),
    mae  = mean(abs(err))
  )
}

make_year_blocks <- function(year, K = 5) {
  yrs <- sort(unique(year))
  split(yrs, cut(seq_along(yrs), breaks = K, labels = FALSE))
}

make_lambda_seq <- function(X, y, alpha, penalty_factor) {
  fit0 <- glmnet(X, y,
                 alpha = alpha,
                 penalty.factor = penalty_factor,
                 intercept = FALSE,
                 standardize = TRUE)
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
    
    # guards
    if (sum(tr) < 50 || sum(va) < 20) next
    
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

compute_shrunk_lsdv_WRDS <- function(in.sample.df, shrink_type = c("EBMLE", "URE"),
                                     centering = c("gen", "0")) {
  
  shrink_type <- match.arg(shrink_type)
  centering   <- match.arg(centering)
  
  in.sample.df <- in.sample.df %>%
    mutate(gvkey = as.factor(gvkey))
  
  fit <- lm(lev_book ~ L1_lev + size + prof + tang + gvkey,
            data = in.sample.df)
  
  ids   <- levels(in.sample.df$gvkey)
  coefs <- coef(fit)
  V     <- vcov(fit)
  
  alpha_hat <- numeric(length(ids))
  names(alpha_hat) <- ids
  
  alpha_hat[ids[1]] <- coefs["(Intercept)"]
  
  if (length(ids) >= 2) {
    for (j in 2:length(ids)) {
      cname <- paste0("gvkey", ids[j])
      alpha_hat[ids[j]] <- coefs["(Intercept)"] + coefs[cname]
    }
  }
  
  M_list <- vector("list", length(ids))
  names(M_list) <- ids
  
  M_list[[ids[1]]] <- matrix(V["(Intercept)", "(Intercept)"], 1, 1)
  
  if (length(ids) >= 2) {
    for (j in 2:length(ids)) {
      cname <- paste0("gvkey", ids[j])
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

# -----------------------------
# 2) Sparse design for glmnet
# -----------------------------
build_Xy_glmnet <- function(df, firm_levels = NULL) {
  df <- df %>%
    mutate(gvkey = as.factor(gvkey)) %>%
    filter(is.finite(lev_book), is.finite(L1_lev), is.finite(size),
           is.finite(prof), is.finite(tang), !is.na(gvkey), !is.na(fyear))
  
  if (!is.null(firm_levels)) {
    df$gvkey <- factor(df$gvkey, levels = firm_levels)
  } else {
    firm_levels <- levels(df$gvkey)
  }
  
  X <- sparse.model.matrix(~ 0 + L1_lev + size + prof + tang + gvkey, data = df)
  y <- df$lev_book
  
  cn <- colnames(X)
  penalty_factor <- ifelse(grepl("^gvkey", cn), 1, 0)
  
  list(X = X, y = y, df = df,
       firm_levels = firm_levels,
       penalty_factor = penalty_factor)
}

# -----------------------------
# 3) Differences on FULL balanced panel (one dataset for AH/GMM)
# -----------------------------
build_diffs_full <- function(df_full) {
  df_full %>%
    arrange(gvkey, fyear) %>%
    group_by(gvkey) %>%
    mutate(
      lagY   = dplyr::lag(lev_book, 1),
      DY     = lev_book - dplyr::lag(lev_book, 1),
      LDY    = dplyr::lag(DY, 1),
      L2Y    = dplyr::lag(lagY, 1),
      D_size = size - dplyr::lag(size, 1),
      D_prof = prof - dplyr::lag(prof, 1),
      D_tang = tang - dplyr::lag(tang, 1)
    ) %>%
    ungroup() %>%
    filter(is.finite(lagY), is.finite(DY), is.finite(LDY), is.finite(L2Y),
           is.finite(D_size), is.finite(D_prof), is.finite(D_tang))
}

# -----------------------------
# 4) GMM MC-style level prediction on selected years
#     yhat = lagY + (g*LDY + b'DX_t)
# -----------------------------
predict_gmm_levels_mcstyle_WRDS <- function(gmm_obj, df_full, years_keep) {
  
  b <- coef(gmm_obj)
  nms <- names(b)
  
  # gamma / AR(1) coefficient name in pgmm
  g_idx <- which(grepl("lag\\(lev_book, 1\\)", nms) | grepl("^lag\\(lev_book, 1\\)$", nms))
  if (length(g_idx) == 0) {
    stop("Cannot find AR(1) coefficient in coef(gmm_obj). Inspect names(coef(gmm_obj)).")
  }
  g <- unname(b[g_idx[1]])
  
  needed <- c("size", "prof", "tang")
  if (!all(needed %in% nms)) {
    stop("Missing GMM coefficients for: ", paste(setdiff(needed, nms), collapse = ", "))
  }
  bb <- b[needed]
  
  dfD <- df_full %>%
    arrange(gvkey, fyear) %>%
    group_by(gvkey) %>%
    mutate(
      lagY   = dplyr::lag(lev_book, 1),
      DY     = lev_book - dplyr::lag(lev_book, 1),
      LDY    = dplyr::lag(DY, 1),
      D_size = size - dplyr::lag(size, 1),
      D_prof = prof - dplyr::lag(prof, 1),
      D_tang = tang - dplyr::lag(tang, 1)
    ) %>%
    ungroup() %>%
    filter(is.finite(lagY), is.finite(LDY),
           is.finite(D_size), is.finite(D_prof), is.finite(D_tang),
           fyear %in% years_keep)
  
  dYhat <- as.numeric(
    g * dfD$LDY +
      bb["size"] * dfD$D_size +
      bb["prof"] * dfD$D_prof +
      bb["tang"] * dfD$D_tang
  )
  
  yhat <- as.numeric(dfD$lagY + dYhat)
  
  list(y_true = dfD$lev_book, yhat = yhat, df_used = dfD)
}

# ============================================================
# 5) Load, clean, and build the balanced panel
# ============================================================
df0 <- read.csv(input_file)

df1 <- df0 %>%
  mutate(
    datadate = as.Date(datadate),
    fyear    = as.integer(fyear),
    sic      = as.integer(sich)
  ) %>%
  filter(!is.na(gvkey), !is.na(fyear)) %>%
  # Exclude financials and utilities
  filter(!(sic >= 6000 & sic <= 6999)) %>%
  filter(!(sic >= 4900 & sic <= 4999)) %>%
  arrange(gvkey, fyear, datadate) %>%
  group_by(gvkey) %>%
  mutate(has_consecutive = any(fyear - dplyr::lag(fyear) == 1, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(has_consecutive) %>%
  dplyr::select(-has_consecutive) %>%
  mutate(
    debt     = dltt + dlc,
    lev_book = debt / at,
    size     = log(at),
    tang     = ppent / at,
    prof     = ib / at
  ) %>%
  filter(is.finite(lev_book), is.finite(size), is.finite(tang), is.finite(prof))

summary(df1$lev_book) # outliers

df1 <- df1 %>%
  mutate(
    lev_book = pmin(debt / at, 1)   # leverage ratios above one capped at one
  )


summary(df1$lev_book) # outliers


df_panel <- df1 %>%
  arrange(gvkey, fyear, datadate) %>%
  group_by(gvkey) %>%
  mutate(
    L1_lev   = dplyr::lag(lev_book, 1),
    year_gap = fyear - dplyr::lag(fyear, 1)
  ) %>%
  ungroup() %>%
  filter(is.finite(L1_lev), year_gap == 1)

rm(df0, df1)

years_full <- sort(unique(df_panel$fyear))
T_full <- length(years_full)

cat("\nFULL years:", min(years_full), "-", max(years_full), " | T =", T_full, "\n")
cat("Pre-balance firms:", length(unique(df_panel$gvkey)), " | obs:", nrow(df_panel), "\n")

firms_balanced <- df_panel %>%
  distinct(gvkey, fyear) %>%
  count(gvkey) %>%
  filter(n == T_full) %>%
  pull(gvkey)

df_panel_bal <- df_panel %>% filter(gvkey %in% firms_balanced)
rm(df_panel)

cat("Balanced firms:", length(unique(df_panel_bal$gvkey)), " | obs:", nrow(df_panel_bal), "\n")

# Precompute one full-differences dataset (used for AH IS and OOS evaluation)
df_diffs_full <- build_diffs_full(df_panel_bal)

# ============================================================
# 6) Build IS/OOS samples
# ============================================================
in.sample  <- df_panel_bal %>% filter(fyear >= IS_START,  fyear <= IS_END)
out.sample <- df_panel_bal %>% filter(fyear >= OOS_START, fyear <= OOS_END)

cat("\n====================================================\n")
cat("SPLIT:\n")
cat("IS :", IS_START, "-", IS_END,  " | T_IS =", length(unique(in.sample$fyear)), "\n")
cat("OOS:", OOS_START, "-", OOS_END, " | H    =", length(unique(out.sample$fyear)), "\n")
cat("N firms:", length(unique(in.sample$gvkey)), "\n")
cat("====================================================\n\n")

# pdata for IS
pin <- pdata.frame(in.sample, index = c("gvkey","fyear"))

# ============================================================
# 7) Estimate models on IS
# ============================================================
fml <- lev_book ~ L1_lev + size + prof + tang

POLS <- plm(fml, data = pin, model = "pooling")
FE   <- plm(fml, data = pin, model = "within", effect = "individual")

# ---- AH on IS diffs (from full diffs)
df_AH_IS <- df_diffs_full %>% filter(fyear >= IS_START, fyear <= IS_END)
pin_AH   <- pdata.frame(df_AH_IS, index = c("gvkey","fyear"))
fml_AH   <- DY ~ LDY - 1 + D_size + D_prof + D_tang | L2Y + D_size + D_prof + D_tang
AH <- plm(fml_AH, data = pin_AH, model = "pooling")

# ---- GMM instruments: lag 2:(T_IS-2) similar spirit to MC
yrs_is <- sort(unique(in.sample$fyear))
Ttr <- length(yrs_is)
max_lag <- max(2, Ttr - 2)

gmm_formula <- paste0(
  "lev_book ~ lag(lev_book,1) + size + prof + tang | lag(lev_book,2:", max_lag, ")"
)
fml_gmm <- as.formula(gmm_formula)

GMM1 <- pgmm(fml_gmm, data = pin,
             effect = "individual",
             model = "onestep",
             transformation = "d",
             collapse = GMM_COLLAPSE)

GMM2 <- pgmm(fml_gmm, data = pin,
             effect = "individual",
             model = "twosteps",
             transformation = "d",
             collapse = GMM_COLLAPSE)

EBMLE <- compute_shrunk_lsdv_WRDS(in.sample, shrink_type = "EBMLE", centering = "gen")
URE   <- compute_shrunk_lsdv_WRDS(in.sample, shrink_type = "URE",   centering = "gen")

# ---- Shrinkage: rolling-origin CV inside IS
dat_in <- build_Xy_glmnet(in.sample)
X_in <- dat_in$X
y_in <- dat_in$y
years_vec <- dat_in$df$fyear
firm_levels <- dat_in$firm_levels
p_factor <- dat_in$penalty_factor

blocks <- make_year_blocks(years_vec, K = K_CV)

# LASSO
lseq_lasso <- make_lambda_seq(X_in, y_in, alpha = 1, penalty_factor = p_factor)
lam_lasso  <- rolling_origin_select_lambda(X_in, y_in, years_vec, blocks,
                                           alpha = 1, penalty_factor = p_factor,
                                           lambda_seq = lseq_lasso, use_1se = USE_1SE)
LASSO <- glmnet(X_in, y_in, alpha = 1, lambda = lam_lasso,
                penalty.factor = p_factor, intercept = FALSE, standardize = TRUE)

# RIDGE
lseq_ridge <- make_lambda_seq(X_in, y_in, alpha = 0, penalty_factor = p_factor)
lam_ridge  <- rolling_origin_select_lambda(X_in, y_in, years_vec, blocks,
                                           alpha = 0, penalty_factor = p_factor,
                                           lambda_seq = lseq_ridge, use_1se = USE_1SE)
RIDGE <- glmnet(X_in, y_in, alpha = 0, lambda = lam_ridge,
                penalty.factor = p_factor, intercept = FALSE, standardize = TRUE)

# EN
alpha_en <- 0.5
lseq_en <- make_lambda_seq(X_in, y_in, alpha = alpha_en, penalty_factor = p_factor)
lam_en  <- rolling_origin_select_lambda(X_in, y_in, years_vec, blocks,
                                        alpha = alpha_en, penalty_factor = p_factor,
                                        lambda_seq = lseq_en, use_1se = USE_1SE)
ELASTIC_NET <- glmnet(X_in, y_in, alpha = alpha_en, lambda = lam_en,
                      penalty.factor = p_factor, intercept = FALSE, standardize = TRUE)

cat("GMM max lag:", max_lag, "\n")
cat("Lambdas: lasso", round(lam_lasso, 6),
    "| ridge", round(lam_ridge, 6),
    "| en", round(lam_en, 6), "\n\n")

# ============================================================
# 8) Evaluation (MC definitions)
# ============================================================
models <- c("POLS","FE","AH","GMM1","GMM2","EBMLE","URE",
            "LASSO","RIDGE","ELASTIC_NET")

in_bias <- setNames(rep(NA_real_, length(models)), models)
in_var  <- setNames(rep(NA_real_, length(models)), models)
in_mse  <- setNames(rep(NA_real_, length(models)), models)
in_mae  <- setNames(rep(NA_real_, length(models)), models)

out_bias <- setNames(rep(NA_real_, length(models)), models)
out_var  <- setNames(rep(NA_real_, length(models)), models)
out_mse  <- setNames(rep(NA_real_, length(models)), models)
out_mae  <- setNames(rep(NA_real_, length(models)), models)


# ---------- IN-SAMPLE ----------
# POLS
b <- coef(POLS)
yhat <- as.numeric(b["(Intercept)"] +
                     b["L1_lev"]*in.sample$L1_lev +
                     b["size"]*in.sample$size +
                     b["prof"]*in.sample$prof +
                     b["tang"]*in.sample$tang)
m <- eval_metrics_mc(in.sample$lev_book, yhat)
in_bias["POLS"] <- m["bias"]; in_var["POLS"] <- m["var"]; in_mse["POLS"] <- m["mse"]; in_mae["POLS"]  <- m["mae"]

# FE
fe_i <- fixef(FE, effect="individual")
b <- coef(FE)
yhat <- as.numeric(fe_i[as.character(in.sample$gvkey)] +
                     b["L1_lev"]*in.sample$L1_lev +
                     b["size"]*in.sample$size +
                     b["prof"]*in.sample$prof +
                     b["tang"]*in.sample$tang)
m <- eval_metrics_mc(in.sample$lev_book, yhat)
in_bias["FE"] <- m["bias"]; in_var["FE"] <- m["var"]; in_mse["FE"] <- m["mse"]; in_mae["FE"]  <- m["mae"]

# AH (levels): yhat = lagY + DYhat
df_AH_IS$DYhat <- as.numeric(fitted(AH))
df_AH_IS$yhat  <- as.numeric(df_AH_IS$lagY) + df_AH_IS$DYhat
m <- eval_metrics_mc(df_AH_IS$lev_book, df_AH_IS$yhat)
in_bias["AH"] <- m["bias"]; in_var["AH"] <- m["var"]; in_mse["AH"] <- m["mse"]; in_mae["AH"]  <- m["mae"]

# GMM (MC-style): yhat = lagY + (g*LDY + b'DX)
res1 <- predict_gmm_levels_mcstyle_WRDS(GMM1, df_panel_bal, years_keep = yrs_is)
m <- eval_metrics_mc(res1$y_true, res1$yhat)
in_bias["GMM1"] <- m["bias"]; in_var["GMM1"] <- m["var"]; in_mse["GMM1"] <- m["mse"]; in_mae["GMM1"]  <- m["mae"]

res2 <- predict_gmm_levels_mcstyle_WRDS(GMM2, df_panel_bal, years_keep = yrs_is)
m <- eval_metrics_mc(res2$y_true, res2$yhat)
in_bias["GMM2"] <- m["bias"]; in_var["GMM2"] <- m["var"]; in_mse["GMM2"] <- m["mse"]; in_mae["GMM2"]  <- m["mae"]

# EBMLE
b <- coef(EBMLE$fit_lm)

yhat <- as.numeric(
  EBMLE$alpha_shrunk[as.character(in.sample$gvkey)] +
    b["L1_lev"] * in.sample$L1_lev +
    b["size"]   * in.sample$size +
    b["prof"]   * in.sample$prof +
    b["tang"]   * in.sample$tang
)

m <- eval_metrics_mc(in.sample$lev_book, yhat)
in_bias["EBMLE"] <- m["bias"]; in_var["EBMLE"] <- m["var"]; in_mse["EBMLE"] <- m["mse"]; in_mae["EBMLE"]  <- m["mae"]

# URE
b <- coef(URE$fit_lm)

yhat <- as.numeric(
  URE$alpha_shrunk[as.character(in.sample$gvkey)] +
    b["L1_lev"] * in.sample$L1_lev +
    b["size"]   * in.sample$size +
    b["prof"]   * in.sample$prof +
    b["tang"]   * in.sample$tang
)

m <- eval_metrics_mc(in.sample$lev_book, yhat)
in_bias["URE"] <- m["bias"]; in_var["URE"] <- m["var"]; in_mse["URE"] <- m["mse"]; in_mae["URE"]  <- m["mae"]

# Shrinkage IS (aligned)
dat_in_eval <- build_Xy_glmnet(in.sample, firm_levels = firm_levels)
X_in_eval <- dat_in_eval$X
y_in_eval <- dat_in_eval$y

yhat <- as.numeric(predict(LASSO, newx = X_in_eval, s = lam_lasso))
m <- eval_metrics_mc(y_in_eval, yhat)
in_bias["LASSO"] <- m["bias"]; in_var["LASSO"] <- m["var"]; in_mse["LASSO"] <- m["mse"]; in_mae["LASSO"]  <- m["mae"]

yhat <- as.numeric(predict(RIDGE, newx = X_in_eval, s = lam_ridge))
m <- eval_metrics_mc(y_in_eval, yhat)
in_bias["RIDGE"] <- m["bias"]; in_var["RIDGE"] <- m["var"]; in_mse["RIDGE"] <- m["mse"]; in_mae["RIDGE"]  <- m["mae"]

yhat <- as.numeric(predict(ELASTIC_NET, newx = X_in_eval, s = lam_en))
m <- eval_metrics_mc(y_in_eval, yhat)
in_bias["ELASTIC_NET"] <- m["bias"]; in_var["ELASTIC_NET"] <- m["var"]; in_mse["ELASTIC_NET"] <- m["mse"]; in_mae["ELASTIC_NET"]  <- m["mae"]


# ---------- OUT-OF-SAMPLE ----------
# POLS
b <- coef(POLS)
yhat <- as.numeric(b["(Intercept)"] +
                     b["L1_lev"]*out.sample$L1_lev +
                     b["size"]*out.sample$size +
                     b["prof"]*out.sample$prof +
                     b["tang"]*out.sample$tang)
m <- eval_metrics_mc(out.sample$lev_book, yhat)
out_bias["POLS"] <- m["bias"]; out_var["POLS"] <- m["var"]; out_mse["POLS"] <- m["mse"]; out_mae["POLS"]  <- m["mae"]

# FE (plug-in alpha_i from IS)
fe_i <- fixef(FE, effect="individual")
b <- coef(FE)
yhat <- as.numeric(fe_i[as.character(out.sample$gvkey)] +
                     b["L1_lev"]*out.sample$L1_lev +
                     b["size"]*out.sample$size +
                     b["prof"]*out.sample$prof +
                     b["tang"]*out.sample$tang)
m <- eval_metrics_mc(out.sample$lev_book, yhat)
out_bias["FE"] <- m["bias"]; out_var["FE"] <- m["var"]; out_mse["FE"] <- m["mse"]; out_mae["FE"]  <- m["mae"]

# AH OOS (levels): yhat = lagY + (b'Z)
df_AH_OOS <- df_diffs_full %>% filter(fyear >= OOS_START, fyear <= OOS_END)
bA <- coef(AH)

df_AH_OOS$DYhat <- as.numeric(
  bA["LDY"]     * df_AH_OOS$LDY +
    bA["D_size"] * df_AH_OOS$D_size +
    bA["D_prof"] * df_AH_OOS$D_prof +
    bA["D_tang"] * df_AH_OOS$D_tang
)
df_AH_OOS$yhat <- df_AH_OOS$lagY + df_AH_OOS$DYhat

m <- eval_metrics_mc(df_AH_OOS$lev_book, df_AH_OOS$yhat)
out_bias["AH"] <- m["bias"]; out_var["AH"] <- m["var"]; out_mse["AH"] <- m["mse"]; out_mae["AH"]  <- m["mae"]

# GMM OOS (MC-style reconstruction)
yrs_out <- sort(unique(out.sample$fyear))

res1 <- predict_gmm_levels_mcstyle_WRDS(GMM1, df_panel_bal, years_keep = yrs_out)
m <- eval_metrics_mc(res1$y_true, res1$yhat)
out_bias["GMM1"] <- m["bias"]; out_var["GMM1"] <- m["var"]; out_mse["GMM1"] <- m["mse"]; out_mae["GMM1"]  <- m["mae"]

res2 <- predict_gmm_levels_mcstyle_WRDS(GMM2, df_panel_bal, years_keep = yrs_out)
m <- eval_metrics_mc(res2$y_true, res2$yhat)
out_bias["GMM2"] <- m["bias"]; out_var["GMM2"] <- m["var"]; out_mse["GMM2"] <- m["mse"]; out_mae["GMM2"]  <- m["mae"]

# EBMLE OOS
b <- coef(EBMLE$fit_lm)

yhat <- as.numeric(
  EBMLE$alpha_shrunk[as.character(out.sample$gvkey)] +
    b["L1_lev"] * out.sample$L1_lev +
    b["size"]   * out.sample$size +
    b["prof"]   * out.sample$prof +
    b["tang"]   * out.sample$tang
)

m <- eval_metrics_mc(out.sample$lev_book, yhat)
out_bias["EBMLE"] <- m["bias"]; out_var["EBMLE"] <- m["var"]; out_mse["EBMLE"] <- m["mse"]; out_mae["EBMLE"]  <- m["mae"]

# URE OOS
b <- coef(URE$fit_lm)

yhat <- as.numeric(
  URE$alpha_shrunk[as.character(out.sample$gvkey)] +
    b["L1_lev"] * out.sample$L1_lev +
    b["size"]   * out.sample$size +
    b["prof"]   * out.sample$prof +
    b["tang"]   * out.sample$tang
)

m <- eval_metrics_mc(out.sample$lev_book, yhat)
out_bias["URE"] <- m["bias"]; out_var["URE"] <- m["var"]; out_mse["URE"] <- m["mse"]; out_mae["URE"]  <- m["mae"]

# Shrinkage OOS (same firm_levels)
dat_out_eval <- build_Xy_glmnet(out.sample, firm_levels = firm_levels)
X_out_eval <- dat_out_eval$X
y_out_eval <- dat_out_eval$y

yhat <- as.numeric(predict(LASSO, newx = X_out_eval, s = lam_lasso))
m <- eval_metrics_mc(y_out_eval, yhat)
out_bias["LASSO"] <- m["bias"]; out_var["LASSO"] <- m["var"]; out_mse["LASSO"] <- m["mse"]; out_mae["LASSO"]  <- m["mae"]

yhat <- as.numeric(predict(RIDGE, newx = X_out_eval, s = lam_ridge))
m <- eval_metrics_mc(y_out_eval, yhat)
out_bias["RIDGE"] <- m["bias"]; out_var["RIDGE"] <- m["var"]; out_mse["RIDGE"] <- m["mse"]; out_mae["RIDGE"]  <- m["mae"]

yhat <- as.numeric(predict(ELASTIC_NET, newx = X_out_eval, s = lam_en))
m <- eval_metrics_mc(y_out_eval, yhat)
out_bias["ELASTIC_NET"] <- m["bias"]; out_var["ELASTIC_NET"] <- m["var"]; out_mse["ELASTIC_NET"] <- m["mse"]; out_mae["ELASTIC_NET"]  <- m["mae"]

# ============================================================
# 9) Print tables (MC style)
# ============================================================
table_IS <- data.frame(
  Estimator = models,
  Bias = as.numeric(in_bias[models]),
  Var  = as.numeric(in_var[models]),
  RMSE = sqrt(as.numeric(in_mse[models])),
  MAE  = as.numeric(in_mae[models])
) %>%
  mutate(across(c(Bias, Var, RMSE,MAE), ~ round(.x, 6)))

table_OOS <- data.frame(
  Estimator = models,
  Bias = as.numeric(out_bias[models]),
  Var  = as.numeric(out_var[models]),
  RMSE = sqrt(as.numeric(out_mse[models])),
  MAE  = as.numeric(out_mae[models])
) %>%
  mutate(across(c(Bias, Var, RMSE,MAE), ~ round(.x, 6)))

cat("\n--- In-sample table (rotated) ---\n")
print(table_IS)

cat("\n--- Out-of-sample table (rotated) ---\n")
print(table_OOS)



# ============================================================
# 10) Save outputs
# ============================================================
write.csv(table_IS,
          file = file.path(table_dir, "table_IS_1991_2010.csv"),
          row.names = FALSE)

write.csv(table_OOS,
          file = file.path(table_dir, "table_OOS_2011_2024.csv"),
          row.names = FALSE)


# ============================================================
# 11) LaTeX tables
# ============================================================

latex_table <- function(tab, caption, label = NULL){
  
  tab2 <- tab
  
  # Bold minimum RMSE
  min_rmse <- which.min(as.numeric(tab2$RMSE))
  tab2$RMSE[min_rmse] <- paste0(
    "\\textbf{",
    sprintf("%.4f", as.numeric(tab2$RMSE[min_rmse])),
    "}"
  )
  
  # Bold minimum MAE
  min_mae <- which.min(as.numeric(tab2$MAE))
  tab2$MAE[min_mae] <- paste0(
    "\\textbf{",
    sprintf("%.4f", as.numeric(tab2$MAE[min_mae])),
    "}"
  )
  
  tab2 %>%
    kbl(
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      caption = caption,
      label = label,
      align = "lcccc",
      col.names = c("Estimator", "Bias", "Variance", "RMSE", "MAE")
    ) %>%
    kable_styling(latex_options = c("hold_position"))
}

latex_table(
  table_IS,
  caption = "In-sample forecasting performance",
  label = "insample_compustat"
)

latex_table(
  table_OOS,
  caption = "Out-of-sample forecasting performance",
  label = "oos_compustat"
)


# ============================================================
# 12) Figure: relative out-of-sample performance
# ============================================================
# table_OOS must contain: Estimator, Bias, Variance, RMSE, MAE

plot_data <- table_OOS %>%
  mutate(
    RMSE = as.numeric(RMSE),
    MAE  = as.numeric(MAE)
  ) %>%
  mutate(
    rmse_fe = RMSE[Estimator == "FE"],
    mae_fe  = MAE[Estimator == "FE"],
    rel_RMSE = 100 * (RMSE / rmse_fe - 1),
    rel_MAE  = 100 * (MAE  / mae_fe  - 1)
  ) %>%
  select(Estimator, rel_RMSE, rel_MAE) %>%
  pivot_longer(
    cols = c(rel_RMSE, rel_MAE),
    names_to = "Metric",
    values_to = "Relative"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      rel_RMSE = "Relative RMSE",
      rel_MAE  = "Relative MAE"
    ),
    Group = case_when(
      Estimator == "FE" ~ "FE benchmark",
      
      Estimator %in% c("POLS", "OLS", "LSDV") ~ "OLS-type",
      
      Estimator %in% c("AH", "GMM1", "GMM2") ~ "IV-type",
      
      Estimator %in% c("EBMLE", "URE") ~ "Model-based shrinkage",
      
      Estimator %in% c("LASSO", "RIDGE", "ELASTIC_NET", "Ridge", "Elastic Net") ~ "ML shrinkage",
      
      TRUE ~ "Other"
    ),
    Estimator = factor(
      Estimator,
      levels = c("GMM2", "GMM1", "AH", "FE", "URE", "EBMLE",
                 "POLS", "ELASTIC_NET", "Elastic Net", "LASSO", "RIDGE", "Ridge")
    )
  )



p <- ggplot(plot_data, aes(x = Relative, y = Estimator, fill = Group)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_text(
    aes(
      label = sprintf("%.1f", Relative),
      hjust = ifelse(Relative >= 0, -0.15, 1.15)
    ),
    size = 3
  ) +
  facet_wrap(~ Metric, ncol = 2, scales = "free_x") +
  scale_fill_manual(
    values = c(
      "OLS-type" = "#C44E52",              # rojo
      "IV-type" = "#55A868",               # verde
      "Model-based shrinkage" = "#1F4E79",# azul oscuro
      "ML shrinkage" = "#8ECAE6",          # celeste
      "FE benchmark" = "black"
    )
) +
  labs(
    #title = "Relative Out-of-Sample Forecasting Performance",
    #subtitle = "Benchmark: Fixed Effects",
    x = "Relative difference with respect to FE (%)",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_dir, "figure_relative_rmse_mae.eps"),
  plot = p,
  device = cairo_ps,
  width = 8,
  height = 5
)
