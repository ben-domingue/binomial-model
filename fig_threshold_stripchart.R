# src/fig_threshold_stripchart.R
#
# Manuscript figure for Section 4.3 ("Two illustrative datasets"): GPCM
# threshold structure for every item in the two datasets used to illustrate
# good vs. poor binomial fit. Promoted from
# scratch/diag_illustrative_thresholds.R after this stripchart proved
# clearer than the original CRF/ERF comparison (fig_empirical_illustrative.R,
# kept but no longer referenced in main.tex) at showing *why* the binomial
# approximation succeeds/fails: the good-fit dataset's items all have tight,
# regularly-spaced thresholds; the poor-fit dataset's items show a
# consistent asymmetric pattern (small gap at the low end, widening gaps
# toward the top).
#
# Dataset choice (see fig_empirical_illustrative.R header for the two
# rejected candidates and why): good fit = assessment_time_fournier_2025_bscs
# (omega_c ~ -0.009 vs. GPCM); poor fit = LSHS-E_Quidaja_2022 (omega_c ~
# 0.21 vs. GPCM), chosen after verifying all items in both datasets have
# positive discrimination and non-degenerate, properly-ordered thresholds
# (no sparse-category estimation artifacts).
#
# Layout: two panels side by side (not stacked) to minimize vertical space;
# items within each panel ordered by threshold asymmetry (ascending), so the
# eye reads left-to-right, top-to-bottom as "most regular" to "most
# asymmetric."
#
# Output: fig_threshold_stripchart.pdf (Overleaf project root)

suppressPackageStartupMessages({
  library(irw)
  library(mirt)
  library(dplyr)
  library(ggplot2)
})
source("/home/ben/Dropbox/projects/binomial_model/src/irw_fetch_cached.R")

set.seed(20260717)
OUT <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model/figs/fig_threshold_stripchart.pdf"

fit_thresholds <- function(tab) {
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
  ni <- ncol(resp_mat)
  message("  N=", nrow(resp_mat), "  items=", ni, "  K=", K)

  mod_gpcm <- mirt(resp_mat, 1, itemtype = rep("gpcm", ni), verbose = FALSE,
                    technical = list(NCYCLES = 2000))

  item_b <- lapply(seq_len(ni), function(j) {
    p <- coef(mod_gpcm, IRTpars = TRUE)[[j]]
    as.numeric(p[1, grep("^b", colnames(p))])
  })

  thresholds_long <- bind_rows(lapply(seq_len(ni), function(j) {
    b <- sort(item_b[[j]])
    gaps <- diff(b)
    asym <- if (length(gaps) >= 2 && all(gaps > 0)) log(max(gaps) / min(gaps)) else NA_real_
    data.frame(item = j, threshold = b, k = seq_along(b), asym = asym)
  }))
  list(N = nrow(resp_mat), ni = ni, K = K, df = thresholds_long)
}

GOOD <- fit_thresholds("assessment_time_fournier_2025_bscs")
POOR <- fit_thresholds("LSHS-E_Quidaja_2022")

mean_asym <- function(res) res$df |> distinct(item, asym) |> pull(asym) |> mean(na.rm = TRUE)
cat(sprintf("GOOD mean log gap-ratio = %.3f\n", mean_asym(GOOD)))
cat(sprintf("POOR mean log gap-ratio = %.3f\n", mean_asym(POOR)))

GOOD_label <- "Good fit: assessment_time_fournier_2025_bscs"
POOR_label <- "Poor fit: LSHS-E_Quidaja_2022"

all_df <- bind_rows(
  GOOD$df |> mutate(row = GOOD_label),
  POOR$df |> mutate(row = POOR_label)
)

# Order items within each panel by asymmetry (ascending): most regular first.
item_order <- all_df |>
  distinct(row, item, asym) |>
  arrange(row, asym) |>
  mutate(item_rank = row_number(), .by = row)

all_df <- all_df |>
  left_join(item_order, by = c("row", "item", "asym"))

p <- ggplot(all_df, aes(x = threshold, y = item_rank, colour = factor(k))) +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 1.6) +
  scale_colour_brewer(palette = "Set1", name = "threshold\nindex") +
  scale_y_continuous(breaks = NULL) +
  facet_wrap(~row, scales = "free_y", ncol = 2) +
  labs(x = "GPCM threshold (b)", y = "item") +
  theme_bw(base_size = 10) +
  theme(strip.text       = element_text(size = 9),
        legend.position  = "bottom",
        panel.spacing    = unit(1, "lines"))

ggsave(OUT, p, width = 8, height = 3.4)
message("Saved ", OUT)
