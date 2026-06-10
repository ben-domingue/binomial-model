library(ggplot2)
library(patchwork)

# ── Parameters ────────────────────────────────────────────────────────────────
K      <- 5L
a      <- 1.5
s      <- log(8)          # upper-wider asymmetry
THETAS <- c(-1, 1)      # two highlighted ability values
theta_grid <- seq(-3, 3, by = 0.05)

COL_GPCM <- "#333333"
COL_NP   <- "#C53A35"
CAT_COLS <- c("#185FA5","#1D9E75","#BA7517","#C53A35","#7F77DD")

# ── Helpers ───────────────────────────────────────────────────────────────────
make_thresholds <- function(K, s) {
  n <- K - 1L; rng <- 3.0
  if (abs(s) < 1e-9 || n <= 2L)
    return(seq(-rng/2, rng/2, length.out = n))
  wi   <- exp(s * seq(0, n-2) / (n - 2))
  gaps <- wi / sum(wi) * rng
  b    <- cumsum(c(-rng/2, gaps[-length(gaps)]))
  c(b, b[n-1] + gaps[n-1])
}

gpcm_pr <- function(theta, a, b) {
  K  <- length(b) + 1L
  ln <- matrix(0, length(theta), K)
  for (k in seq_len(K-1)) ln[, k+1] <- ln[, k] + a*(theta - b[k])
  ln <- ln - apply(ln, 1, max)
  ex <- exp(ln); ex / rowSums(ex)
}

dbinom_mat <- function(p, K)
  outer(p, 0:(K-1), function(pp, k) dbinom(k, K-1, pp))

np_binom_pr <- function(theta, a, b, K) {
  gp   <- gpcm_pr(theta, a, b)
  emn  <- as.numeric(gp %*% 0:(K-1))
  phat <- pmin(pmax(emn / (K-1), 1e-10), 1-1e-10)
  dbinom_mat(phat, K)
}

erf <- function(pr, K) as.numeric(pr %*% 0:(K-1))

# ── Compute grid data ─────────────────────────────────────────────────────────
b  <- make_thresholds(K, s)
gp <- gpcm_pr(theta_grid, a, b)
np <- np_binom_pr(theta_grid, a, b, K)

grid_df <- data.frame(
  th   = theta_grid,
  erfG = erf(gp, K),
  erfN = erf(np, K)
)
for (k in 0:(K-1)) {
  grid_df[[paste0("crfG", k)]] <- gp[, k+1]
  grid_df[[paste0("crfN", k)]] <- np[, k+1]
}

# ── Point data at the two θ values ────────────────────────────────────────────
pt_rows <- lapply(THETAS, function(th) {
  gp_pt <- as.numeric(gpcm_pr(th, a, b))
  np_pt <- as.numeric(np_binom_pr(th, a, b, K))
  data.frame(
    theta = th,
    k     = factor(0:(K-1)),
    prob  = c(gp_pt, np_pt),
    model = rep(c("GPCM", "NP spline"), each = K)
  )
})
pt_df <- do.call(rbind, pt_rows)
pt_df$model <- factor(pt_df$model, levels = c("GPCM", "NP spline"))
pt_df$theta_lab <- factor(
  paste0("θ = ", pt_df$theta),
  levels = paste0("θ = ", THETAS)
)

# ── vline helper ──────────────────────────────────────────────────────────────
vlines <- data.frame(th = THETAS)

# ── Panel A: ERF ──────────────────────────────────────────────────────────────
pA <- ggplot(grid_df, aes(x = th)) +
  geom_vline(data = vlines, aes(xintercept = th),
             linetype = "dashed", colour = "grey70", linewidth = 0.5) +
  geom_line(aes(y = erfG, colour = "GPCM"),    linewidth = 1.1) +
  geom_line(aes(y = erfN, colour = "NP spline"), linewidth = 1.0,
            linetype = "dashed") +
  scale_colour_manual(values = c("GPCM" = COL_GPCM, "NP spline" = COL_NP),
                      name = NULL) +
  scale_x_continuous(breaks = -2:2) +
  scale_y_continuous(limits = c(0, K-1), breaks = 0:(K-1)) +
  labs(x = expression(theta), y = expression(E*"["*X*"|"*theta*"]"),
       title = "A  Expected response function") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10, face = "bold"))

# ── Panel B: CRF ──────────────────────────────────────────────────────────────
crf_long <- do.call(rbind, lapply(0:(K-1), function(k) {
  rbind(
    data.frame(th = theta_grid, prob = grid_df[[paste0("crfG", k)]],
               k = factor(k), model = "GPCM"),
    data.frame(th = theta_grid, prob = grid_df[[paste0("crfN", k)]],
               k = factor(k), model = "NP spline")
  )
}))
crf_long$model <- factor(crf_long$model, levels = c("GPCM", "NP spline"))

pB <- ggplot(crf_long, aes(x = th, y = prob, colour = k)) +
  geom_vline(data = vlines, aes(xintercept = th),
             linetype = "dashed", colour = "grey70", linewidth = 0.5,
             inherit.aes = FALSE) +
  geom_line(aes(linetype = model), linewidth = 0.9) +
  scale_colour_manual(values = CAT_COLS[seq_len(K)], name = "k") +
  scale_linetype_manual(values = c("GPCM" = "solid", "NP spline" = "dashed"),
                        name = NULL) +
  scale_x_continuous(breaks = -2:2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(theta), y = expression(P*"("*X == k*"|"*theta*")"),
       title = "B  Category response functions") +
  guides(colour   = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2)) +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10, face = "bold"))

# ── Panel C: bar chart of P(X=k|θ) at two θ values ───────────────────────────
pC <- ggplot(pt_df, aes(x = k, y = prob, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  facet_wrap(~ theta_lab, ncol = 2) +
  scale_fill_manual(values = c("GPCM" = COL_GPCM, "NP spline" = COL_NP),
                    name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Response category k",
       y = expression(P*"("*X == k*"|"*theta*")"),
       title = "C  Category probabilities at marked θ values") +
  theme_bw(base_size = 11) +
  theme(legend.position   = "bottom",
        panel.grid.minor  = element_blank(),
        strip.background  = element_rect(fill = "grey92"),
        plot.title        = element_text(size = 10, face = "bold"))

# ── Assemble ──────────────────────────────────────────────────────────────────
p_final <- (pA | pB | pC) +
  plot_annotation(
    title    = paste0("GPCM vs. NP spline  ·  K = ", K,
                      ",  a = ", a, ",  s = ln8  (upper-wider thresholds)"),
    subtitle = paste0("Thresholds: ", paste(round(b, 2), collapse = ", "),
                      "   ·   Vertical lines mark θ = ",
                      paste(THETAS, collapse = " and ")),
    theme = theme(plot.title    = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey40"))
  ) &
  theme(legend.justification = "center")

ggsave("/home/ben/Dropbox/Apps/Overleaf/binomial_model/fig_kl_category.pdf",
       p_final, width = 11, height = 4, device = cairo_pdf)
cat("saved fig_kl_category.pdf\n")
