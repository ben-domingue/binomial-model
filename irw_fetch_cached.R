# src/irw_fetch_cached.R
#
# Thin wrapper around irw::irw_fetch() that caches the raw fetched
# data frame to /tmp so repeated runs (e.g. while iterating on
# fig_empirical_illustrative.R / scratch/diag_illustrative_thresholds.R)
# don't re-hit Redivis. Temporary convenience, not part of the
# reproducibility pipeline -- /tmp is not durable storage.

irw_fetch_cached <- function(tab, cache_dir = "/tmp/irw_cache") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_path <- file.path(cache_dir, paste0(tab, ".rds"))
  if (file.exists(cache_path)) {
    message("  [cache] ", tab, " <- ", cache_path)
    return(readRDS(cache_path))
  }
  df <- irw::irw_fetch(tab)
  saveRDS(df, cache_path)
  message("  [cache] ", tab, " -> ", cache_path)
  df
}
