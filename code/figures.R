# =============================================================================
# Figures 1 and 2
# "Shrinkage and Persistence in Dynamic Panel Forecasting"
# Magdalena Cornejo and Walter Sosa-Escudero
#
# Figure 1: limiting mean squared forecast error along the ridge path
#           (sigma2_alpha = sigma2_eps = 1), T = 4 and 10.
# Figure 2: Monte Carlo mean squared forecast error as a function of the degree
#           of shrinkage (N = 100, sigma2_eps = 1, phi = 1, 1,000 replications),
#           with the limits implied by Proposition 2 and Theorem 1.
#
# Output: results/figures/fig_ridge_risk.pdf, results/figures/fig_shrinkage_path.pdf
# Only base R is used.
# =============================================================================

rm(list = ls())
OUT_DIR <- "results/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Closed-form limits (Section 3)
# -----------------------------------------------------------------------------
consts_T <- function(g, T) {
  k <- seq_len(T - 1)
  list(v = 1 + (2 / T) * sum((T - k) * g^k), h = sum((1 - g^(T - k)) / (1 - g)),
       c = g * (1 - g^T) / (T * (1 - g)), p = (1 - g^T) / (T * (1 - g)))
}
# Asymptotic bias b_T(omega) of the ridge estimator of gamma (Proposition 2)
bias_T <- function(om, g, T, s2 = 1, sa2 = 1) {
  k <- consts_T(g, T); su2 <- s2 / (1 - g^2)
  (-s2 * k$h / T + om * (sa2 / (1 - g) + s2 * k$h / T^2)) /
    (su2 * (T - k$v) + om * (sa2 / (1 - g)^2 + su2 * k$v / T))
}
# Limiting forecast risk R_T(kappa, g) (Theorem 1)
risk_T <- function(kappa, gf, g, T, s2 = 1, sa2 = 1) {
  k <- consts_T(g, T); su2 <- s2 / (1 - g^2); d <- g - gf
  s2 + kappa^2 * s2 / T + (1 - kappa)^2 * (1 - gf)^2 * sa2 / (1 - g)^2 +
    d^2 * su2 * (1 - 2 * kappa * k$c + kappa^2 * k$v / T) -
    2 * d * kappa * s2 * (k$p - kappa * k$h / T^2)
}
# Known-persistence benchmark R^K_T (Proposition 1)
r_known <- function(T, s2 = 1, sa2 = 1) s2 + sa2 * s2 / (T * sa2 + s2)
# Best linear forecast R^*_T (Proposition 4)
r_star <- function(g, T, s2 = 1, sa2 = 1) {
  phi <- sa2 / s2; om2 <- sa2 / (1 + phi * (1 + g) / (1 - g))
  s2 + om2 * (s2 / T) / (om2 + s2 / T)
}

GAMMAS <- c(0.2, 0.5, 0.8)
LTY    <- c(1, 2, 3)          # solid, dashed, dotted

# -----------------------------------------------------------------------------
# Figure 1
# -----------------------------------------------------------------------------
lams <- exp(seq(log(1e-3), log(1e4), length.out = 3000))
pdf(file.path(OUT_DIR, "fig_ridge_risk.pdf"), width = 10, height = 3.9)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
curves <- lapply(c(4, 10), function(T) sapply(GAMMAS, function(g)
  sapply(lams, function(l) risk_T(T / (T + l), g + bias_T(T * l / (T + l), g, T), g, T))))
ylim <- range(unlist(curves), sapply(c(4, 10), r_known))
for (j in 1:2) {
  T <- c(4, 10)[j]
  matplot(lams, curves[[j]], type = "l", lty = LTY, col = "black", lwd = 1.2, log = "x",
          ylim = ylim, xlab = expression(lambda), ylab = if (j == 1) "Limiting MSFE" else "",
          main = bquote(T == .(T)))
  abline(h = r_known(T), col = "grey55", lwd = 1.8)
  abline(v = 1, col = "grey60", lty = 4)
  if (j == 2) legend("topleft", bty = "n",
                     legend = c(expression(gamma == 0.2), expression(gamma == 0.5),
                                expression(gamma == 0.8), expression(R[T]^K)),
                     lty = c(LTY, 1), lwd = c(1.2, 1.2, 1.2, 1.8),
                     col = c("black", "black", "black", "grey55"))
}
dev.off()

# -----------------------------------------------------------------------------
# Figure 2
# -----------------------------------------------------------------------------
N <- 100; REPS <- 1000
kg <- seq(0, 1, by = 0.05)        # simulated grid of kappa
kf <- seq(0, 1, by = 0.005)       # grid for the limits
set.seed(2026)

sim_path <- function(T, g) {
  R_r <- R_w <- R_k <- numeric(length(kg))
  for (rep in 1:REPS) {
    a <- rnorm(N)
    y <- matrix(0, N, T + 2)
    y[, 1] <- a / (1 - g) + rnorm(N, sd = 1 / sqrt(1 - g^2))     # stationary start
    for (t in 2:(T + 2)) y[, t] <- a + g * y[, t - 1] + rnorm(N)
    x  <- y[, 1:T]; yy <- y[, 2:(T + 1)]                         # lag and outcome, t = 1..T
    xb <- rowMeans(x); yb <- rowMeans(yy)
    W10 <- sum((x - xb) * (yy - yb)); W11 <- sum((x - xb)^2)
    B10 <- sum((xb - mean(xb)) * (yb - mean(yb))); B11 <- sum((xb - mean(xb))^2)
    gW <- W10 / W11
    target <- y[, T + 2]; last <- y[, T + 1]
    for (j in seq_along(kg)) {
      k  <- kg[j]; om <- T * (1 - k)
      gr <- (W10 + om * B10) / (W11 + om * B11)                    # ridge (Lemma 1)
      q  <- (last - mean(xb)) - k * (xb - mean(xb))
      base <- mean(yb) + k * (yb - mean(yb))                       # Lemma 2
      R_r[j] <- R_r[j] + mean((target - base - gr * q)^2)          # ridge path
      R_w[j] <- R_w[j] + mean((target - base - gW * q)^2)          # post-within
      R_k[j] <- R_k[j] + mean((target - base - g  * q)^2)          # known gamma
    }
  }
  list(r = R_r / REPS, w = R_w / REPS, k = R_k / REPS)
}

pdf(file.path(OUT_DIR, "fig_shrinkage_path.pdf"), width = 11, height = 6.6)
layout(rbind(matrix(1:6, 2, 3, byrow = TRUE), 7), heights = c(1, 1, 0.18))
par(mar = c(4, 4.2, 2, 1))
for (T in c(5, 10)) for (g in GAMMAS) {
  s <- sim_path(T, g)
  lim_r <- sapply(kf, function(k) risk_T(k, g + bias_T(T * (1 - k), g, T), g, T))
  lim_w <- sapply(kf, function(k) risk_T(k, g + bias_T(0, g, T), g, T))
  lim_k <- sapply(kf, function(k) risk_T(k, g, g, T))
  Rs <- r_star(g, T); RF <- risk_T(1, g + bias_T(0, g, T), g, T)
  ylim <- c(min(Rs, s$r, s$k) - 0.03, 1.3 * RF)                    # truncated at 1.3 x R_FE
  plot(1 - kf, lim_r, type = "n", ylim = ylim,
       xlab = if (T == 10) expression("Degree of shrinkage, " * 1 - kappa) else "",
       ylab = if (g == 0.2) "Mean squared forecast error" else "",
       main = bquote(list(T == .(T), gamma == .(g))), cex.main = 1)
  grid(col = "grey90")
  lines(1 - kf, lim_r, lwd = 1.4)
  points(1 - kg, s$r, pch = 16, cex = 0.6)
  lines(1 - kf, lim_w, lty = 2, col = "grey35", lwd = 1.2)
  points(1 - kg, s$w, pch = 22, cex = 0.6, col = "grey35", bg = "white")
  lines(1 - kf, lim_k, lty = 3, lwd = 1.4)
  points(1 - kg, s$k, pch = 24, cex = 0.6, bg = "white")
  abline(h = Rs, lty = 4, col = "grey50")
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", horiz = TRUE, bty = "n", cex = 0.95,
       legend = c("Ridge path", expression("Post-within (" * hat(gamma)[W] * ")"),
                  expression("Known " * gamma), expression("Best linear forecast " * R[T]^"*")),
       lty = c(1, 2, 3, 4), lwd = c(1.4, 1.2, 1.4, 1), col = c("black", "grey35", "black", "grey50"),
       pch = c(16, 22, 24, NA))
dev.off()
