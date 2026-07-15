## Extends Figure 1 (fig_isoerf_distance2_altcontour.pdf / iso.grm.fig) from two
## discrimination snapshots (a'=0.7, a'=1.5) into a continuous sweep over a',
## to make the "broad at low a', compact at high a'" claim in SI Section 1.4
## (si:isoerf.summary) directly visible and quantitative rather than inferred
## from two static panels.
##
## Same K=3 reference (a=1, b=(0,1)) and grid ranges as isoerf_distance2_altcontour.R,
## and the same common near-isoERF threshold (RMISE=0.15) used there, so the
## two figures are read together.

library(ggplot2)
library(patchwork)

logis <- function(x) 1 / (1 + exp(-x))

erf_grm <- function(theta, a, b) rowSums(outer(theta, b, function(t, bk) logis(a * (t - bk))))
erf_gpcm <- function(theta, a, b) {
  sapply(theta, function(t) {
    u1 <- exp(a * (t - b[1])); u2 <- exp(a * (t - b[1]) + a * (t - b[2]))
    (u1 + 2 * u2) / (1 + u1 + u2)
  })
}
erf_binom <- function(theta, a, b) 2 * logis(a * (theta - b[1]))

refs <- list(
  GRM      = list(a = 1, b = c(0, 1),    erf_fn = erf_grm),
  GPCM     = list(a = 1, b = c(0, 1),    erf_fn = erf_gpcm),
  Binomial = list(a = 1, b = c(0.5, NA), erf_fn = erf_binom)
)

erf_dist <- function(erf_ref_fn, erf_fn, a_prime, b_prime) {
  integrate(function(theta) {
    (erf_ref_fn(theta) - erf_fn(theta, a_prime, b_prime))^2 * dnorm(theta)
  }, lower = -5, upper = 5, subdivisions = 200)$value
}

b1_vals <- seq(-2, 2, length.out = 60)
b2_vals <- seq(-1, 3, length.out = 60)
b1_vals_binom <- seq(-2, 3, length.out = 400)

grids_2d <- list(
  GRM      = subset(expand.grid(b1 = b1_vals, b2 = b2_vals), b1 < b2),
  GPCM     = expand.grid(b1 = b1_vals, b2 = b2_vals)
)

a_seq <- seq(0.4, 3.0, by = 0.1)
contour_level_common <- 0.15   # shared with Figure 1 / iso.grm.fig

sweep_model <- function(model) {
  ref <- refs[[model]]
  efn <- ref$erf_fn
  erf_ref_fn <- function(theta) efn(theta, ref$a, ref$b)

  if (model == "Binomial") {
    n_total <- length(b1_vals_binom)
    out <- lapply(a_seq, function(ap) {
      dist <- sapply(b1_vals_binom, function(b1) erf_dist(erf_ref_fn, efn, ap, c(b1, NA)))
      rmise <- sqrt(dist)
      data.frame(a_prime = ap, min_rmise = min(rmise),
                 frac_below = mean(rmise <= contour_level_common))
    })
  } else {
    g <- grids_2d[[model]]
    n_total <- nrow(g)
    out <- lapply(a_seq, function(ap) {
      dist <- mapply(function(b1, b2) erf_dist(erf_ref_fn, efn, ap, c(b1, b2)), g$b1, g$b2)
      rmise <- sqrt(dist)
      data.frame(a_prime = ap, min_rmise = min(rmise),
                 frac_below = mean(rmise <= contour_level_common))
    })
  }
  res <- do.call(rbind, out)
  res$model <- model
  cat(sprintf("done %s\n", model))
  res
}

res <- do.call(rbind, lapply(c("GRM", "GPCM", "Binomial"), sweep_model))
res$model <- factor(res$model, levels = c("GRM", "GPCM", "Binomial"))

model_cols <- c(GRM = "#1D9E75", GPCM = "#D85A30", Binomial = "#378ADD")

theme_iso <- theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 10.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = "bottom")

wrap <- function(s, w = 95) paste(strwrap(s, width = w), collapse = "\n")

# Panel A: size of the near-isoERF neighborhood (fraction of the sampled
# admissible domain with RMISE <= 0.15) as a function of candidate a'.
# This is the direct, continuous version of "broad at a'=0.7, compact at
# a'=1.5" -- Figure 1 shows two vertical slices through this curve.
# Zero-fraction stretches (no admissible candidate within threshold at all)
# are left as gaps rather than floored, so an absent line segment is itself
# the signal, not an artifact of the log scale.
res$frac_plot <- ifelse(res$frac_below > 0, res$frac_below, NA)
p_frac <- ggplot(res, aes(x = a_prime, y = frac_plot, color = model)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = c(0.7, 1.5), linetype = "dashed", color = "grey75") +
  geom_line(linewidth = 1.1) +
  scale_y_log10(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = model_cols, name = NULL) +
  labs(title = "Size of the near-isoERF neighborhood vs. candidate discrimination",
       subtitle = wrap("Fraction of the sampled admissible (b1', b2') domain with RMISE <= 0.15 (log scale); gaps mean no candidate at that a' falls within threshold at all. Vertical dashed lines mark a'=0.7 and a'=1.5 from Figure 1 (main text); dotted line marks a=1 (truth)."),
       x = expression(a*"'"), y = "Neighborhood fraction") +
  theme_iso

# Panel B: minimum achievable RMISE (best-compensating candidate) vs a' --
# shows *why* the binomial cannot exhibit the same transition: with only one
# free threshold it cannot compensate for a mismatched a', so even its best
# candidate stays far from the reference ERF once a' departs from 1.
p_min <- ggplot(res, aes(x = a_prime, y = min_rmise, color = model)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = contour_level_common, linetype = "dashed", color = "grey75") +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = model_cols, name = NULL) +
  labs(title = "Best-case ERF match vs. candidate discrimination",
       subtitle = wrap("Minimum RMISE over candidate thresholds b' at each a'; dashed line marks the RMISE=0.15 threshold used throughout."),
       x = expression(a*"'"), y = "min RMISE") +
  theme_iso

p_final <- (p_frac / p_min) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("/home/ben/Dropbox/Apps/Overleaf/binomial_model/fig_isoerf_neighborhood_vs_a.pdf",
       p_final, width = 9.5, height = 9)
cat("saved\n")
