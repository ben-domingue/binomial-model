# src/diag_flex03_crfs.R
#
# Single-rep diagnostic at flex = 0.3.
# Generates NP DGP data, fits GPCM (mirt) and 1PL (GLM with EAP thetas),
# saves true + estimated params, and plots true vs estimated CRFs.
#
# Uses mirt for GPCM to match the actual analysis pipeline.
#
# Outputs:
#   src/diag_flex03_crfs.pdf   -- 10-item x 2-model CRF overlay
#   src/diag_flex03_data.rds   -- true params, fits, X_mat, thetas

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(mirt)
})

# ── Parameters ────────────────────────────────────────────────────────────────
set.seed(42)
N_PERSONS       <- 1000L
N_ITEMS         <- 10L
K               <- 5L
FLEX            <- 0.3
LNORM_SDLOG     <- 0.5
THRESH_NOISE_SD <- 0.3
theta_g         <- seq(-3.5, 3.5, length.out = 301L)

# ── Model helpers ─────────────────────────────────────────────────────────────

gpcm_probs <- function(theta, a, b) {
  # b: length K-1 vector of thresholds
  K  <- length(b) + 1L
  ln <- matrix(0, nrow = length(theta), ncol = K)
  for (k in seq_len(K - 1L))
    ln[, k + 1L] <- ln[, k] + a * (theta - b[k])
  ln <- ln - apply(ln, 1, max)
  ex <- exp(ln); ex / rowSums(ex)
}

binom_probs <- function(theta, a, b, K) {
  p <- pmin(pmax(plogis(a * theta - b), 1e-10), 1 - 1e-10)
  outer(p, 0:(K - 1L), function(pp, k) dbinom(k, K - 1L, pp))
}

rand_mono_spline <- function(n_knots = 7L) {
  x_k <- seq(-4, 4, length.out = n_knots)
  y_k <- sort(runif(n_knots))
  splinefun(x_k, y_k, method = "monoH.FC")
}

np_probs <- function(theta, a_j, b_j, K, flex, splines_j) {
  n       <- K - 1L
  theta_c <- pmin(pmax(theta, -4), 4)
  gpcm_cat <- gpcm_probs(theta, a_j, b_j)
  gpcm_cum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n))
    gpcm_cum[, k] <- 1 - rowSums(gpcm_cat[, seq_len(k), drop = FALSE])
  rand_cum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n))
    rand_cum[, k] <- pmin(pmax(splines_j[[k]](theta_c), 0), 1)
  cum <- (1 - flex) * gpcm_cum + flex * rand_cum
  if (n > 1L) cum <- t(apply(cum, 1, sort, decreasing = TRUE))
  probs <- cbind(
    1 - cum[, 1],
    if (n > 1L) cum[, seq_len(n - 1L), drop = FALSE] -
                cum[, seq_len(n - 1L) + 1L, drop = FALSE],
    cum[, n]
  )
  pmax(probs, 1e-10)
}

# ── Generate data ──────────────────────────────────────────────────────────────
cat("Generating NP DGP data (flex =", FLEX, ")...\n")
theta_true <- rnorm(N_PERSONS)
shifts     <- runif(N_ITEMS, -1, 1)

true_params <- vector("list", N_ITEMS)
X_mat       <- matrix(NA_integer_, N_PERSONS, N_ITEMS)

for (j in seq_len(N_ITEMS)) {
  a_j       <- rlnorm(1L, meanlog = 0, sdlog = LNORM_SDLOG)
  b_j       <- seq(-1.5, 1.5, length.out = K - 1L) + shifts[j] +
               rnorm(K - 1L, 0, THRESH_NOISE_SD)
  splines_j <- lapply(seq_len(K - 1L), function(k) rand_mono_spline())

  probs_j <- np_probs(theta_true, a_j, b_j, K, FLEX, splines_j)
  probs_j <- probs_j / rowSums(probs_j)
  X_j     <- apply(probs_j, 1, function(p) sample.int(K, 1L, prob = p) - 1L)

  tries <- 0L
  while (length(unique(X_j)) < K && tries < 5L) {
    shifts[j]  <- runif(1, -1, 1)
    b_j        <- seq(-1.5, 1.5, length.out = K - 1L) + shifts[j] +
                  rnorm(K - 1L, 0, THRESH_NOISE_SD)
    splines_j  <- lapply(seq_len(K - 1L), function(k) rand_mono_spline())
    probs_j    <- np_probs(theta_true, a_j, b_j, K, FLEX, splines_j)
    probs_j    <- probs_j / rowSums(probs_j)
    X_j        <- apply(probs_j, 1, function(p) sample.int(K, 1L, prob = p) - 1L)
    tries      <- tries + 1L
  }

  true_params[[j]] <- list(a = a_j, b = b_j, splines = splines_j)
  X_mat[, j]       <- X_j
}

cat("Response distributions per item:\n")
for (j in 1:N_ITEMS)
  cat(sprintf("  Item %2d: %s\n", j,
              paste(table(factor(X_mat[,j], 0:(K-1))), collapse = "  ")))

cat("\nTrue discriminations (a_j):\n")
a_true <- sapply(true_params, `[[`, "a")
print(round(a_true, 3))

# ── Fit GPCM with mirt ────────────────────────────────────────────────────────
cat("\nFitting GPCM with mirt...\n")
X_df     <- as.data.frame(X_mat)
mod_gpcm <- mirt(X_df, 1, itemtype = "gpcm", verbose = FALSE)
theta_eap <- as.numeric(fscores(mod_gpcm, method = "EAP")[, 1])

# Extract per-item GPCM parameters (IRT parameterization)
gpcm_params <- lapply(seq_len(N_ITEMS), function(j) {
  p <- coef(mod_gpcm, IRTpars = TRUE)[[j]]
  list(a = as.numeric(p[1, "a"]),
       b = as.numeric(p[1, paste0("b", seq_len(K - 1L))]))
})

cat("GPCM discrimination estimates:\n")
print(round(sapply(gpcm_params, `[[`, "a"), 3))

# ── Fit 1PL with GLM (using EAP thetas) ──────────────────────────────────────
cat("\nFitting 1PL (GLM with EAP thetas)...\n")
binom_params <- lapply(seq_len(N_ITEMS), function(j) {
  df  <- data.frame(X = X_mat[, j], theta = theta_eap)
  mod <- glm(cbind(X, K - 1L - X) ~ offset(theta),
             data = df, family = binomial("logit"))
  list(a = 1.0, b = -coef(mod)[["(Intercept)"]])
})

cat("1PL b estimates:\n")
print(round(sapply(binom_params, `[[`, "b"), 3))

# ── Print parameter comparison ────────────────────────────────────────────────
cat("\nParameter comparison (true vs GPCM vs 1PL):\n")
cat(sprintf("%-4s  %-6s  %-6s  |  %-6s  %-6s  |  %-6s  %-6s  |  %-6s\n",
            "item","a_true","a_gpcm","b1_true","b1_gpcm","b4_true","b4_gpcm","b_1pl"))
for (j in seq_len(N_ITEMS)) {
  tp <- true_params[[j]]; gp <- gpcm_params[[j]]; bp <- binom_params[[j]]
  cat(sprintf("%-4d  %-6.3f  %-6.3f  |  %-6.3f  %-6.3f  |  %-6.3f  %-6.3f  |  %-6.3f\n",
              j, tp$a, gp$a, tp$b[1], gp$b[1], tp$b[K-1], gp$b[K-1], bp$b))
}

# ── Build CRF data ─────────────────────────────────────────────────────────────
cat("\nBuilding CRF data...\n")
crf_rows <- lapply(seq_len(N_ITEMS), function(j) {
  tp <- true_params[[j]]
  gp <- gpcm_params[[j]]
  bp <- binom_params[[j]]

  true_mat <- np_probs(theta_g, tp$a, tp$b, K, FLEX, tp$splines)
  gpcm_mat <- gpcm_probs(theta_g, gp$a, gp$b)
  b1pl_mat <- binom_probs(theta_g, bp$a, bp$b, K)

  bind_rows(
    as.data.frame(true_mat) |>
      setNames(paste0("k", 0:(K-1))) |>
      mutate(theta = theta_g, source = "True", item = j),
    as.data.frame(gpcm_mat) |>
      setNames(paste0("k", 0:(K-1))) |>
      mutate(theta = theta_g, source = "GPCM", item = j),
    as.data.frame(b1pl_mat) |>
      setNames(paste0("k", 0:(K-1))) |>
      mutate(theta = theta_g, source = "Binom 1PL", item = j)
  )
})

crf_df <- bind_rows(crf_rows) |>
  pivot_longer(starts_with("k"), names_to = "category", values_to = "prob") |>
  mutate(
    source   = factor(source, levels = c("True", "GPCM", "Binom 1PL")),
    item_lab = factor(paste0("Item ", item), levels = paste0("Item ", 1:N_ITEMS)),
    category = factor(category, levels = paste0("k", 0:(K-1)),
                      labels = paste0("k=", 0:(K-1)))
  )

# 3-column layout: True | GPCM | Binom 1PL — one clean set of curves per panel
plot_df <- crf_df |>
  mutate(panel = factor(source, levels = c("True", "GPCM", "Binom 1PL")))

# ── Plot ──────────────────────────────────────────────────────────────────────
cat_colors <- c("k=0" = "#1b7837", "k=1" = "#2166ac", "k=2" = "#d6604d",
                "k=3" = "#762a83", "k=4" = "#e08214")

p <- ggplot(plot_df,
            aes(x = theta, y = prob, colour = category, group = category)) +
  geom_line(linewidth = 0.75) +
  scale_colour_manual(values = cat_colors, name = "Category") +
  facet_grid(item_lab ~ panel) +
  labs(x = expression(theta), y = expression(P(X == k ~"|"~ theta)),
       title = sprintf(
         "CRFs: True | GPCM fit | Binom 1PL fit  —  flex = %.1f, N = %d, K = %d",
         FLEX, N_PERSONS, K)) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92"),
        panel.spacing.x  = unit(0.5, "lines"),
        panel.spacing.y  = unit(0.2, "lines"))

ggsave("src/diag_flex03_crfs.pdf", p, width = 10, height = 18)
message("Saved src/diag_flex03_crfs.pdf")

# ── Save data ──────────────────────────────────────────────────────────────────
saveRDS(
  list(true_params  = true_params,
       gpcm_params  = gpcm_params,
       binom_params = binom_params,
       theta_true   = theta_true,
       theta_eap    = theta_eap,
       X_mat        = X_mat,
       flex         = FLEX, K = K, N = N_PERSONS),
  "src/diag_flex03_data.rds"
)
message("Saved src/diag_flex03_data.rds")
