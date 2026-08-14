# src/scratch/diag_illustrative_thresholds.R
#
# Diagnostic (not a manuscript figure): sanity-check the GPCM threshold
# structure for ALL items in the two datasets used in
# fig_empirical_illustrative.R (good fit: assessment_time_fournier_2025_bscs;
# poor fit: LSHS-E_Quidaja_2022), before trusting the "threshold asymmetry"
# story the main figure is meant to illustrate.
#
# Fits GPCM on the FULL data (no holdout -- we just want a clean read of the
# threshold structure, not a predictive check) and plots every item's
# K-1 thresholds as a stripplot, one row per item, faceted by dataset.
#
# Output: diag_illustrative_thresholds.pdf (this directory)

suppressPackageStartupMessages({
  library(irw)
  library(mirt)
  library(dplyr)
  library(ggplot2)
})
source("/home/ben/Dropbox/projects/binomial_model/src/irw_fetch_cached.R")

set.seed(20260717)

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
  item_a <- sapply(seq_len(ni), function(j) {
    p <- coef(mod_gpcm, IRTpars = TRUE)[[j]]
    as.numeric(p[1, "a"])
  })

  bind_rows(lapply(seq_len(ni), function(j) {
    b <- sort(item_b[[j]])
    gaps <- diff(b)
    asym <- if (length(gaps) >= 2 && all(gaps > 0)) log(max(gaps) / min(gaps)) else NA_real_
    data.frame(tab = tab, item = j, a = item_a[j], threshold = b,
               k = seq_along(b), asym = asym)
  }))
}

GOOD_df <- fit_thresholds("assessment_time_fournier_2025_bscs")
POOR_df <- fit_thresholds("LSHS-E_Quidaja_2022")

all_df <- bind_rows(
  GOOD_df |> mutate(row = sprintf("Good fit: assessment_time_fournier_2025_bscs")),
  POOR_df |> mutate(row = sprintf("Poor fit: LSHS-E_Quidaja_2022"))
)

# Order items within each dataset by mean threshold (just for readability)
item_order <- all_df |>
  group_by(row, item) |>
  summarise(mean_b = mean(threshold), asym = first(asym), .groups = "drop") |>
  arrange(row, mean_b)

all_df <- all_df |>
  mutate(item_label = factor(paste(row, item), levels = paste(item_order$row, item_order$item)))

cat("\n--- Per-item mean log gap-ratio (asymmetry) ---\n")
print(item_order |> arrange(row, desc(asym)), n = 30)

p <- ggplot(all_df, aes(x = threshold, y = item_label, colour = factor(k))) +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Set1", name = "threshold\nindex") +
  facet_wrap(~row, scales = "free_y", ncol = 1) +
  labs(x = "GPCM threshold (b)", y = "item (sorted by mean threshold)") +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 7))

ggsave("/home/ben/Dropbox/projects/binomial_model/src/scratch/diag_illustrative_thresholds.pdf",
       p, width = 8, height = 9)
message("Saved diag_illustrative_thresholds.pdf")
