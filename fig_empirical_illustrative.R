# src/fig_empirical_illustrative.R
#
# Two real IRW datasets illustrating when the binomial approximation works
# well vs. poorly, for the new Section 4.3 ("Two illustrative datasets").
#
# Selection (from irw/poly_compare_data/multi_*.rds pairwise IMV cache,
# Binom 2PL (logit) -> GPCM comparison, restricted to K=5 datasets with a
# similar number of items/respondents so the two rows are comparable):
#
#   GOOD FIT: assessment_time_fournier_2025_bscs
#     omega_c ~ -0.009 (binomial ~ties/edges GPCM), mean log gap-ratio ~ 0.84
#     (near-symmetric GPCM thresholds)
#
#   POOR FIT: LSHS-E_Quidaja_2022
#     omega_c ~ 0.20, omega_t ~ 0.10 (GPCM clearly better), mean log
#     gap-ratio ~ 1.9 (more asymmetric GPCM thresholds)
#
# Exact omega_c/gap-ratio values vary slightly run to run (random 80/20
# holdout, no fixed seed guarantee across mirt/EM internals); rerun and
# read the console output for the numbers to cite in prose.
#
# NB: florida_twins_class was tried first and rejected -- 2 of its 10 items
# had non-identified GPCM threshold estimates (e.g. -13 to +16, sparse-
# category degeneracy, not genuine asymmetry). chakraborty2026_BIDR_IRW was
# also tried and rejected -- half its items had negative discrimination
# (unrecoded reverse-keyed BIDR items), which is a scale-coding artifact,
# not a threshold-structure story. LSHS-E_Quidaja_2022 was verified clean:
# all 16 items have positive discrimination and sane, properly-ordered
# thresholds (max |b| = 5.4, vs. e.g. florida_twins_class's 16). See
# scratch/diag_illustrative_thresholds.R for the full per-item threshold
# diagnostic this was checked against.
#
# Both are refit here from scratch (not read from the cache, which only
# stores summary statistics, not fitted objects) using the same
# preprocessing/model-fitting logic as irw/poly_compare_multi_compute.R:
# 80/20 cell holdout, GPCM vs. custom Binom 2PL (logit) via mirt, EAP theta.
#
# For each dataset we plot the item with the largest GPCM-vs-Binom-2PL CRF
# divergence (same convention as fig_single_item_example.R), with the GPCM
# threshold estimates marked as rug ticks on the ERF panel to show threshold
# spacing directly.
#
# Output: fig_empirical_illustrative.pdf (Overleaf project root)

suppressPackageStartupMessages({
  library(irw)
  library(imv)
  library(mirt)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})
source("/home/ben/Dropbox/projects/binomial_model/src/irw_fetch_cached.R")

set.seed(20260717)
OUT <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model/figs/fig_empirical_illustrative.pdf"

MIN_N     <- 200
MAX_ITEMS <- 60L
MAX_N     <- 5000L

# ── Custom binomial 2PL item (matches poly_compare_multi_compute.R) ──────────
binom_P <- function(prob, ncat) {
  size <- ncat - 1L
  prob <- pmin(pmax(prob, 1e-10), 1 - 1e-10)
  outer(prob, 0:size, function(p, k) dbinom(k, size, p))
}
binom_2pl <- createItem("binom_2pl",
  par = c(a = 1, b = 0), est = c(TRUE, TRUE),
  P = function(par, Theta, ncat) binom_P(plogis(par[1] * Theta[, 1] - par[2]), ncat)
)

# ── IMV helpers (identical to poly_compare_multi_compute.R) ──────────────────
imv_c_fn <- function(y, pctt.tab, p1, p2) {
  nn  <- length(pctt.tab); om <- numeric(nn); iis <- 0:(nn - 1)
  for (ii in iis) {
    ns <- om.tmp <- numeric()
    for (jj in iis[-match(ii, iis)]) {
      y2   <- y[y$resp %in% c(ii, jj), ]
      resp <- ifelse(y2$resp == ii, 1, 0)
      p1v  <- y2[[paste0(p1, ii)]] / (y2[[paste0(p1, ii)]] + y2[[paste0(p1, jj)]])
      p2v  <- y2[[paste0(p2, ii)]] / (y2[[paste0(p2, ii)]] + y2[[paste0(p2, jj)]])
      om.tmp[as.character(jj)] <- imv.binary(resp, p1v, p2v)
      ns[as.character(jj)]     <- nrow(y2)
    }
    om[ii + 1] <- sum(om.tmp * ns) / sum(ns)
  }
  sum(om * pctt.tab) / sum(pctt.tab)
}
imv_t_fn <- function(y, pctt.tab, p1, p2) {
  nn <- length(pctt.tab); om <- numeric(nn - 1)
  for (ii in 0:(nn - 2)) {
    resp <- ifelse(y$resp <= ii, 1, 0)
    pr1  <- rowSums(y[, paste0(p1, 0:ii), drop = FALSE])
    pr2  <- rowSums(y[, paste0(p2, 0:ii), drop = FALSE])
    om[ii + 1] <- imv.binary(resp, pr1, pr2)
  }
  wts <- pctt.tab[1:(nn - 1)] / (1 - pctt.tab[nn])
  sum(om * wts) / sum(wts)
}

# ── Fetch + fit one dataset ───────────────────────────────────────────────────
fit_dataset <- function(tab) {
  message("=== ", tab, " ===")
  df <- irw_fetch_cached(tab)
  K_min    <- min(df$resp, na.rm = TRUE)
  df$resp0 <- df$resp - K_min
  resp_wide <- irw_long2resp(df %>% select(-resp) %>% rename(resp = resp0))
  item_cols <- setdiff(names(resp_wide), "id")
  resp_mat  <- as.matrix(resp_wide[, item_cols])

  keep     <- apply(resp_mat, 2, function(x) length(unique(na.omit(x))) >= 2)
  resp_mat <- resp_mat[, keep, drop = FALSE]

  item_K   <- apply(resp_mat, 2, function(x) length(unique(na.omit(x))))
  target_K <- as.integer(names(sort(table(item_K), decreasing = TRUE))[1])
  resp_mat <- resp_mat[, item_K == target_K, drop = FALSE]
  K <- target_K

  if (ncol(resp_mat) > MAX_ITEMS) resp_mat <- resp_mat[, sample(ncol(resp_mat), MAX_ITEMS), drop = FALSE]
  if (nrow(resp_mat) > MAX_N)     resp_mat <- resp_mat[sample(nrow(resp_mat), MAX_N), , drop = FALSE]

  ni <- ncol(resp_mat); N <- nrow(resp_mat)
  message("  N=", N, "  items=", ni, "  K=", K)

  non_missing <- which(!is.na(resp_mat), arr.ind = TRUE)
  holdout_idx <- sample(nrow(non_missing), floor(0.2 * nrow(non_missing)))
  train_mat   <- resp_mat
  train_mat[non_missing[holdout_idx, , drop = FALSE]] <- NA
  test_cells  <- non_missing[holdout_idx, , drop = FALSE]
  obs_vals    <- resp_mat[test_cells]

  mod_gpcm <- mirt(train_mat, 1, itemtype = rep("gpcm", ni), verbose = FALSE,
                    technical = list(NCYCLES = 2000))
  mod_b2pl <- mirt(train_mat, 1, itemtype = rep("binom_2pl", ni), verbose = FALSE,
                    customItems = list(binom_2pl = binom_2pl),
                    technical = list(NCYCLES = 2000))

  # Each model scores persons on its own EAP theta (matches
  # irw/poly_compare_multi_compute.R). Using one model's theta to evaluate
  # the other's item parameters is invalid: the two latent scales are
  # separately identified and need not agree in sign/scale, even though in
  # practice they are highly correlated.
  theta_gpcm <- fscores(mod_gpcm, response.pattern = train_mat, method = "EAP")[, 1]
  theta_b2pl <- fscores(mod_b2pl, response.pattern = train_mat, method = "EAP")[, 1]

  predict_probs <- function(mod, theta_eap) {
    probs <- matrix(NA_real_, nrow(test_cells), K)
    for (j in seq_len(ni)) {
      rows_j <- which(test_cells[, 2] == j)
      if (!length(rows_j)) next
      Theta_j <- matrix(theta_eap[test_cells[rows_j, 1]], ncol = 1)
      pr <- probtrace(extract.item(mod, j), Theta_j)
      if (ncol(pr) == K) probs[rows_j, ] <- pr
    }
    probs
  }
  probs_gpcm <- predict_probs(mod_gpcm, theta_gpcm)
  probs_b2pl <- predict_probs(mod_b2pl, theta_b2pl)

  train_obs <- as.integer(train_mat[!is.na(train_mat)])
  pctt.tab  <- as.numeric(table(factor(train_obs, levels = 0:(K - 1)))) / length(train_obs)

  y_imv <- as.data.frame(probs_b2pl); colnames(y_imv) <- paste0("p1", 0:(K - 1))
  y_imv$resp <- obs_vals
  for (k in 0:(K - 1)) y_imv[[paste0("p2", k)]] <- probs_gpcm[, k + 1]
  y_imv <- y_imv[complete.cases(y_imv), ]
  omega_c <- imv_c_fn(y_imv, pctt.tab, "p1", "p2")
  omega_t <- imv_t_fn(y_imv, pctt.tab, "p1", "p2")
  message(sprintf("  Binom 2PL -> GPCM:  omega_c = %.3f  omega_t = %.3f", omega_c, omega_t))

  # GPCM threshold (b) estimates, per item, for the asymmetry annotation
  item_b <- lapply(seq_len(ni), function(j) {
    p <- coef(mod_gpcm, IRTpars = TRUE)[[j]]
    as.numeric(p[1, grep("^b", colnames(p))])
  })
  item_asym <- sapply(item_b, function(b) {
    gaps <- diff(sort(b)); gaps <- gaps[gaps > 0]
    if (length(gaps) < 2) return(NA_real_)
    log(max(gaps) / min(gaps))
  })

  # The two models' theta scales are separately identified (arbitrary sign);
  # align them before plotting on a shared theta axis. Magnitude (~N(0,1))
  # is comparable by construction of EAP scoring, so only sign needs fixing.
  sign_flip <- sign(cor(theta_gpcm, theta_b2pl, use = "complete.obs"))
  message("  cor(theta_gpcm, theta_b2pl) = ", round(cor(theta_gpcm, theta_b2pl, use = "complete.obs"), 3))

  # Per-item CRF divergence between fitted GPCM and Binom 2PL, over a theta grid
  theta_g <- seq(-3, 3, length.out = 201)
  divergence <- sapply(seq_len(ni), function(j) {
    pg <- probtrace(extract.item(mod_gpcm, j), matrix(theta_g, ncol = 1))
    pb <- probtrace(extract.item(mod_b2pl, j), matrix(sign_flip * theta_g, ncol = 1))
    mean(abs(pg - pb))
  })
  # Guard against picking a non-identified/degenerate GPCM fit (sparse
  # categories can drive thresholds to implausible extremes) as the
  # "most divergent" item -- restrict to items with sane threshold estimates.
  well_identified <- sapply(item_b, function(b) max(abs(b)) < 8)
  if (!any(well_identified)) well_identified <- rep(TRUE, ni)
  divergence[!well_identified] <- -Inf
  item_pick <- which.max(divergence)

  pg <- probtrace(extract.item(mod_gpcm, item_pick), matrix(theta_g, ncol = 1))
  pb <- probtrace(extract.item(mod_b2pl, item_pick), matrix(sign_flip * theta_g, ncol = 1))

  list(tab = tab, N = N, ni = ni, K = K, omega_c = omega_c, omega_t = omega_t,
       item_pick = item_pick, item_b = item_b[[item_pick]],
       mean_log_gap_ratio = mean(item_asym, na.rm = TRUE),
       theta_g = theta_g, probs_gpcm = pg, probs_b2pl = pb)
}

# ── Fit both datasets ─────────────────────────────────────────────────────────
GOOD <- fit_dataset("assessment_time_fournier_2025_bscs")
POOR <- fit_dataset("LSHS-E_Quidaja_2022")

# ── Build long data frame for plotting ────────────────────────────────────────
erf_fn <- function(probs, K) as.numeric(probs %*% 0:(K - 1)) / (K - 1)

build_df <- function(res, row_label) {
  K <- res$K
  bind_rows(
    bind_rows(
      data.frame(theta = res$theta_g, prob = erf_fn(res$probs_gpcm, K), model = "GPCM",      panel = "ERF"),
      data.frame(theta = res$theta_g, prob = erf_fn(res$probs_b2pl, K), model = "Binom 2PL",  panel = "ERF")
    ),
    bind_rows(lapply(0:(K - 1), function(k)
      bind_rows(
        data.frame(theta = res$theta_g, prob = res$probs_gpcm[, k + 1], model = "GPCM",     panel = paste0("k = ", k)),
        data.frame(theta = res$theta_g, prob = res$probs_b2pl[, k + 1], model = "Binom 2PL", panel = paste0("k = ", k))
      )))
  ) |> mutate(row = row_label)
}

plot_df <- bind_rows(
  build_df(GOOD, sprintf("Good fit (omega_c = %.3f)", GOOD$omega_c)),
  build_df(POOR, sprintf("Poor fit (omega_c = %.3f)", POOR$omega_c))
) |>
  mutate(
    row   = factor(row, levels = unique(row)),
    model = factor(model, levels = c("GPCM", "Binom 2PL")),
    panel = factor(panel, levels = c("ERF", paste0("k = ", 0:6)))
  )

thresh_df <- bind_rows(
  data.frame(theta = GOOD$item_b, row = levels(plot_df$row)[1]),
  data.frame(theta = POOR$item_b, row = levels(plot_df$row)[2])
) |> mutate(row = factor(row, levels = levels(plot_df$row)), panel = factor("ERF", levels = levels(plot_df$panel)))

# ── Plot ──────────────────────────────────────────────────────────────────────
# ERF is drawn solid for both models (the two ERFs alone can look deceptively
# similar even when category-level fit diverges sharply -- that's part of
# the point of this figure); category panels keep dashed/dotted so the two
# models stay distinguishable when the curves nearly overlap (top row).
model_colors <- c("GPCM" = "#2166ac", "Binom 2PL" = "#d6604d")
model_ltypes <- c("GPCM" = "dashed",  "Binom 2PL" = "dotted")

plot_df <- plot_df |>
  mutate(ltype = ifelse(panel == "ERF", "solid", as.character(model_ltypes[as.character(model)])))

p <- ggplot(plot_df, aes(x = theta, y = prob, colour = model, linetype = ltype)) +
  geom_line(linewidth = 0.8) +
  geom_rug(data = thresh_df, aes(x = theta), inherit.aes = FALSE,
           sides = "b", colour = "grey30", linewidth = 0.5) +
  scale_colour_manual(values = model_colors, name = NULL) +
  scale_linetype_identity(guide = "none") +
  scale_x_continuous(breaks = c(-2, 0, 2), limits = c(-3, 3), oob = scales::squish) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  facet_grid(row ~ panel, scales = "free_y") +
  labs(x = expression(theta), y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92"),
        strip.text.x     = element_text(size = 9),
        strip.text.y     = element_text(size = 8),
        panel.spacing    = unit(0.5, "lines"))

ggsave(OUT, p, width = 10, height = 5.5)
message("Saved ", OUT)

cat(sprintf("\nGOOD: %s  omega_c=%.4f  omega_t=%.4f  mean_log_gap_ratio=%.3f  item=%d  K=%d\n",
            GOOD$tab, GOOD$omega_c, GOOD$omega_t, GOOD$mean_log_gap_ratio, GOOD$item_pick, GOOD$K))
cat("GOOD item_b (thresholds):", GOOD$item_b, "\n")
cat(sprintf("POOR: %s  omega_c=%.4f  omega_t=%.4f  mean_log_gap_ratio=%.3f  item=%d  K=%d\n",
            POOR$tab, POOR$omega_c, POOR$omega_t, POOR$mean_log_gap_ratio, POOR$item_pick, POOR$K))
cat("POOR item_b (thresholds):", POOR$item_b, "\n")
