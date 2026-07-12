library(ggplot2)
library(patchwork)

# Illustrates a mechanism cited in si.tex (sec:si-oracle): as flex increases
# under the NP DGP, responses concentrate in the two extreme categories
# (X=0 and X=K-1) at the expense of the interior categories, which drives
# the omega_t weighting artefact (non-monotonic gap vs. GPCM/oracle in
# Figure fig:np-imv of the main text).
# Reuses the exact DGP mechanism (helper functions, parameter distributions)
# from gpcm_binomial_sweep_np.R so this is a faithful diagnostic.

set.seed(1)

K       <- 5L
n_bound <- K - 1L
LNORM_SDLOG     <- 0.5
THRESH_NOISE_SD <- 0.3
theta_grid <- seq(-4, 4, by = 0.02)
CAT_COLS   <- c("#185FA5","#1D9E75","#BA7517","#C53A35","#7F77DD")

gpcm_probs <- function(theta, a, b) {
  K  <- length(b) + 1L
  ln <- matrix(0, nrow = length(theta), ncol = K)
  for (k in seq_len(K - 1L)) ln[, k + 1L] <- ln[, k] + a * (theta - b[k])
  ln <- ln - apply(ln, 1, max)
  ex <- exp(ln); ex / rowSums(ex)
}
rand_mono_spline <- function(n_knots = 7L) {
  x_k <- seq(-4, 4, length.out = n_knots)
  y_k <- sort(runif(n_knots))
  splinefun(x_k, y_k, method = "monoH.FC")
}
np_probs <- function(theta, a_j, b_j, K, flex, splines_j) {
  n <- K - 1L
  theta_c <- pmin(pmax(theta, -4), 4)
  gpcm_cat <- gpcm_probs(theta, a_j, b_j)
  gpcm_cum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n)) gpcm_cum[, k] <- 1 - rowSums(gpcm_cat[, seq_len(k), drop = FALSE])
  rand_cum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n)) rand_cum[, k] <- pmin(pmax(splines_j[[k]](theta_c), 0), 1)
  cum <- (1 - flex) * gpcm_cum + flex * rand_cum
  if (n > 1L) cum <- t(apply(cum, 1, sort, decreasing = TRUE))
  probs <- cbind(1 - cum[, 1],
                  if (n > 1L) cum[, seq_len(n - 1L), drop = FALSE] - cum[, seq_len(n - 1L) + 1L, drop = FALSE],
                  cum[, n])
  pmax(probs, 1e-10)
}
draw_item <- function() {
  a_j <- rlnorm(1L, meanlog = 0, sdlog = LNORM_SDLOG)
  shift <- runif(1, -1, 1)
  b_j <- seq(-1.5, 1.5, length.out = n_bound) + shift + rnorm(n_bound, 0, THRESH_NOISE_SD)
  splines_j <- lapply(seq_len(n_bound), function(k) rand_mono_spline())
  list(a = a_j, b = b_j, splines = splines_j)
}

## ---- Panel A: one example item's CRFs at low vs. high flex ----------------
item <- draw_item()
crf_df <- do.call(rbind, lapply(c(0.1, 0.9), function(flex) {
  pr <- np_probs(theta_grid, item$a, item$b, K, flex, item$splines)
  do.call(rbind, lapply(0:(K - 1), function(k) {
    data.frame(theta = theta_grid, prob = pr[, k + 1], k = factor(k),
               flex = factor(paste0("delta == ", flex), levels = c("delta == 0.1", "delta == 0.9")))
  }))
}))
pA <- ggplot(crf_df, aes(theta, prob, colour = k)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~flex, labeller = label_parsed) +
  scale_colour_manual(values = CAT_COLS, name = "k") +
  labs(x = expression(theta), y = expression(P(X==k*"|"*theta)),
       title = "A. Example item: category response functions at low vs. high flex") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        plot.title = element_text(size = 10, face = "bold"))

## ---- Panel B: extreme- vs. interior-category share vs. flex ---------------
## Integrates over theta ~ N(0,1) and many items at each flex value: the
## empirical quantity behind the omega_t weighting artefact.
N_PERSONS_MC <- 3000L
N_ITEMS_MC   <- 60L
flex_vals    <- seq(0, 1, by = 0.1)
theta_mc     <- rnorm(N_PERSONS_MC)

marg_df <- do.call(rbind, lapply(flex_vals, function(flex) {
  cat_counts <- rep(0, K)
  for (i in seq_len(N_ITEMS_MC)) {
    it <- draw_item()
    pr <- np_probs(theta_mc, it$a, it$b, K, flex, it$splines)
    cat_counts <- cat_counts + colSums(pr)
  }
  data.frame(k = factor(0:(K - 1)), share = cat_counts / sum(cat_counts), flex = flex)
}))

share_df <- do.call(rbind, lapply(flex_vals, function(fv) {
  sub <- marg_df[marg_df$flex == fv, ]
  x0      <- sub$share[sub$k == "0"]
  xK1     <- sub$share[sub$k == as.character(K - 1)]
  exterior <- x0 + xK1
  data.frame(flex = fv,
             group = c("Exterior (X = 0 or X = K-1)", "Interior (X = 1, ..., K-2)"),
             share = c(exterior, 1 - exterior))
}))
share_df$group <- factor(share_df$group,
                          levels = c("Exterior (X = 0 or X = K-1)", "Interior (X = 1, ..., K-2)"))

pB <- ggplot(share_df, aes(x = flex, y = share, colour = group)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = c("Exterior (X = 0 or X = K-1)" = "#185FA5",
                                  "Interior (X = 1, ..., K-2)" = "#BA7517"),
                       name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(delta~"(flex)"), y = "Proportion of responses",
       title = "B. Exterior- vs. interior-category share as flex increases") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

p_final <- pA / pB + plot_layout(heights = c(1, 1))

ggsave("/home/ben/Dropbox/Apps/Overleaf/binomial_model/fig_np_extreme_categories.pdf",
       p_final, width = 8.5, height = 8, device = cairo_pdf)
cat("saved fig_np_extreme_categories.pdf\n")
cat("\nExtreme- vs. interior-category share by flex:\n")
print(share_df)
