## diag_threshold_monotonicity.R
## Checks whether the NP DGP produces non-monotone threshold functions P(X>k|theta)
## at high flex — the proposed mechanism behind the non-monotone imv_t_oracle pattern.
##
## Output: src/diag_threshold_monotonicity.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUT <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model/src"

## ── NP DGP functions (copied from gpcm_binomial_sweep_np.R) ──────────────────

gpcm_probs <- function(theta, a, b) {
  K  <- length(b) + 1L
  ln <- matrix(0, nrow = length(theta), ncol = K)
  for (k in seq_len(K - 1L))
    ln[, k + 1L] <- ln[, k] + a * (theta - b[k])
  ln <- ln - apply(ln, 1, max)
  ex <- exp(ln); ex / rowSums(ex)
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

## ── Diagnostic ───────────────────────────────────────────────────────────────

K       <- 5L
n_bound <- K - 1L
theta_g <- seq(-4, 4, length.out = 300)
n_items <- 9L
flex_vals <- c(0.1, 0.5, 0.9)

set.seed(42)
rows <- lapply(flex_vals, function(flex) {
  lapply(seq_len(n_items), function(j) {
    a_j       <- rlnorm(1, 0, 0.5)
    b_j       <- seq(-1.5, 1.5, length.out = n_bound) + rnorm(n_bound, 0, 0.3)
    splines_j <- lapply(seq_len(n_bound), function(k) rand_mono_spline())
    probs     <- np_probs(theta_g, a_j, b_j, K, flex, splines_j)
    ## cumulative P(X > k) for k = 1 ... K-1
    cum <- t(apply(probs, 1, function(p) rev(cumsum(rev(p)))[-1]))
    data.frame(
      theta     = rep(theta_g, n_bound),
      prob      = as.vector(cum),
      threshold = factor(rep(paste0("k=", seq_len(n_bound)), each = length(theta_g))),
      item      = paste0("item ", j),
      flex      = flex
    )
  })
})
df <- do.call(rbind, do.call(c, rows))

## Flag non-monotone curves: P(X>k|theta) should be non-decreasing in theta
non_mono <- df |>
  group_by(flex, item, threshold) |>
  arrange(theta) |>
  summarize(
    is_non_monotone = any(diff(prob) < -1e-6),
    max_decrease    = if (any(diff(prob) < 0)) min(diff(prob)) else 0,
    .groups = "drop"
  )

cat("Non-monotone threshold curves by flex level:\n")
non_mono |>
  group_by(flex) |>
  summarize(
    n_curves        = n(),
    n_non_monotone  = sum(is_non_monotone),
    pct             = round(100 * mean(is_non_monotone), 1),
    worst_decrease  = round(min(max_decrease), 4),
    .groups = "drop"
  ) |> print()

cat("\nPer-item breakdown at flex=0.9:\n")
non_mono |> filter(flex == 0.9) |>
  select(item, threshold, is_non_monotone, max_decrease) |>
  arrange(item, threshold) |> print(n = 40)

## ── Plot: threshold curves coloured by monotonicity ──────────────────────────
df_flag <- df |>
  left_join(non_mono |> select(flex, item, threshold, is_non_monotone),
            by = c("flex", "item", "threshold"))

p <- ggplot(df_flag, aes(x = theta, y = prob,
                         colour = threshold,
                         linetype = is_non_monotone,
                         group = interaction(threshold, item))) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  facet_grid(flex ~ item, labeller = labeller(flex = label_both)) +
  scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"),
                        labels = c("FALSE" = "monotone", "TRUE" = "NON-MONOTONE"),
                        name = NULL) +
  scale_colour_brewer(palette = "Set1", name = "threshold") +
  labs(
    title = "P(X > k | θ) threshold functions from NP DGP  —  dashed = non-monotone in θ",
    subtitle = sprintf("K = %d, %d items, 3 flex levels", K, n_items),
    x = "θ", y = "P(X > k | θ)"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92"),
        legend.position  = "bottom")

ggsave(file.path(OUT, "diag_threshold_monotonicity.pdf"),
       plot = p, width = 14, height = 7, device = cairo_pdf)
message("Saved src/diag_threshold_monotonicity.pdf")
