library(ggplot2)
library(patchwork)

# Illustrates the NP DGP boundary-function construction (Section "Study 2"):
# g_k(theta) = (1-delta) * g_k^GPCM(theta) + delta * s_k(theta),
# for a fixed set of random monotone splines s_k, shown at several delta.

set.seed(42)
K      <- 5L
a      <- 1.5
b      <- seq(-1.5, 1.5, length.out = K - 1L)   # symmetric GPCM baseline thresholds
DELTAS <- c(0, 0.3, 0.7, 1)
theta_grid <- seq(-4, 4, by = 0.02)

COL_K <- c("#185FA5", "#1D9E75", "#BA7517", "#C53A35")

# ── GPCM cumulative boundary functions P(X >= k | theta) ─────────────────────
gpcm_pr <- function(theta, a, b) {
  K  <- length(b) + 1L
  ln <- matrix(0, length(theta), K)
  for (k in seq_len(K - 1L)) ln[, k + 1] <- ln[, k] + a * (theta - b[k])
  ln <- ln - apply(ln, 1, max)
  ex <- exp(ln); ex / rowSums(ex)
}
gp_cat <- gpcm_pr(theta_grid, a, b)
gp_cum <- sapply(seq_len(K - 1L), function(k) 1 - rowSums(gp_cat[, seq_len(k), drop = FALSE]))

# ── Random monotone Hermite splines, one per boundary (fixed draw) ───────────
rand_mono_spline <- function(n_knots = 7L) {
  x_k <- seq(-4, 4, length.out = n_knots)
  y_k <- sort(runif(n_knots))
  splinefun(x_k, y_k, method = "monoH.FC")
}
splines_j <- lapply(seq_len(K - 1L), function(k) rand_mono_spline())
rand_cum  <- sapply(splines_j, function(s) pmin(pmax(s(theta_grid), 0), 1))

# ── Mixture at each delta, re-sorted to enforce g_1 > g_2 > ... > g_{K-1} ────
mix_df <- do.call(rbind, lapply(DELTAS, function(delta) {
  cum <- (1 - delta) * gp_cum + delta * rand_cum
  cum <- t(apply(cum, 1, sort, decreasing = TRUE))
  do.call(rbind, lapply(seq_len(K - 1L), function(k) {
    data.frame(theta = theta_grid, g = cum[, k], k = factor(k),
               delta = factor(paste0("delta == ", delta), levels = paste0("delta == ", DELTAS)))
  }))
}))

p <- ggplot(mix_df, aes(x = theta, y = g, colour = k)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ delta, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = COL_K, name = "k") +
  scale_x_continuous(breaks = -2:2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(theta), y = expression(g[k](theta) == P(X >= k*"|"*theta))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

ggsave("/home/ben/Dropbox/Apps/Overleaf/binomial_model/fig_np_boundary_mixture.pdf",
       p, width = 11, height = 3.2, device = cairo_pdf)
cat("saved fig_np_boundary_mixture.pdf\n")
