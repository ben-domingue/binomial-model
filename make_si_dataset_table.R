suppressPackageStartupMessages({ library(dplyr) })

OUT      <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model"
poly_dir <- "/home/ben/Dropbox/projects/binomial_model/src/irw/poly_compare_data"

normalise_model <- function(m) {
  m <- as.character(m)
  m[m == "1PL+logit"]   <- "Binom 1PL (logit)"
  m[m == "2PL+logit"]   <- "Binom 2PL (logit)"
  m[m == "1PL+cloglog"] <- "Binom 1PL (cloglog)"
  m
}

fs  <- list.files(poly_dir, pattern="^multi_.+\\.rds$", full.names=TRUE)
fs  <- fs[!grepl("multi_(results|summary)\\.rds$", fs)]
raw <- lapply(fs, readRDS)
nms <- sub("^multi_","",sub("\\.rds$","",basename(fs)))
names(raw) <- nms
raw <- Filter(Negate(is.null), raw)
raw <- lapply(raw, function(r) { if(is.null(r[["K"]])) r$K <- r$K_range+1L; r })
results <- Filter(function(r) r$K >= 3, raw)

meta <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(dataset=nm, K=r$K, ni=r$ni, N=r$N, stringsAsFactors=FALSE)
}))

all_rmse <- do.call(rbind, lapply(results, `[[`, "rmse")) |>
  mutate(model = normalise_model(model))
rmse_wide <- all_rmse |>
  filter(model %in% c("GPCM","Binom 2PL (logit)")) |>
  mutate(model = ifelse(model=="GPCM","rmse_gpcm","rmse_binom2pl")) |>
  tidyr::pivot_wider(id_cols=dataset, names_from=model, values_from=rmse)

all_pair <- do.call(rbind, Filter(Negate(is.null), lapply(results, `[[`, "pair_imv"))) |>
  mutate(comparison = sub("1PL\\+logit", "Binom 1PL (logit)",
                     sub("2PL\\+logit", "Binom 2PL (logit)", comparison)))
pair_gpcm <- all_pair |>
  filter(comparison == "Binom 2PL (logit)→GPCM") |>
  select(dataset, imv_c, imv_t)

all_asym <- do.call(rbind, Filter(Negate(is.null), lapply(results, `[[`, "asym")))
if (is.null(all_asym)) all_asym <- data.frame(dataset=character(0), mean_log_gap_ratio=numeric(0))
asym_sel <- all_asym |> select(dataset, mean_log_gap_ratio)

tab <- meta |>
  left_join(rmse_wide, by="dataset") |>
  left_join(pair_gpcm,  by="dataset") |>
  left_join(asym_sel,   by="dataset") |>
  arrange(dataset)

cat("rows:", nrow(tab), "\n")

fmt <- function(x, d=3) ifelse(is.na(x), "--", formatC(x, format="f", digits=d))
esc <- function(x) gsub("_", "\\\\_", x, fixed=FALSE)

lines <- sprintf(
  "%s & %d & %d & %d & %s & %s & %s & %s & %s \\\\",
  esc(tab$dataset), tab$K, tab$ni, tab$N,
  fmt(tab$rmse_gpcm), fmt(tab$rmse_binom2pl),
  fmt(tab$imv_c), fmt(tab$imv_t), fmt(tab$mean_log_gap_ratio, 2)
)

header <- c(
  "\\begin{footnotesize}",
  "\\begin{longtable}{l r r r r r r r r}",
  "\\caption{Per-dataset results for all 436 IRW datasets used in the empirical evaluation. $K$: number of response categories; $n_i$: number of items; $N$: number of respondents. RMSE(GPCM) and RMSE(Binom 2PL): cross-validated RMSE on the rescaled response. $\\omega_c$, $\\omega_t$: pairwise IMV of GPCM over the Binom 2PL (logit), category- and threshold-level (main text Figure~\\ref{fig:pair-imv}). Asym: mean log gap-ratio from fitted GPCM thresholds, an unsigned analogue of the simulation's $|s|$ (main text Figure~\\ref{fig:asym-imv}); `--' where GPCM threshold parameters were not available. Datasets ordered alphabetically.} \\label{tab:si-all-datasets} \\\\",
  "\\toprule",
  "Dataset & $K$ & $n_i$ & $N$ & RMSE(GPCM) & RMSE(Binom 2PL) & $\\omega_c$ & $\\omega_t$ & Asym \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{9}{l}{\\textit{(continued from previous page)}} \\\\",
  "\\toprule",
  "Dataset & $K$ & $n_i$ & $N$ & RMSE(GPCM) & RMSE(Binom 2PL) & $\\omega_c$ & $\\omega_t$ & Asym \\\\",
  "\\midrule",
  "\\endhead",
  "\\midrule",
  "\\multicolumn{9}{r}{\\textit{(continued on next page)}} \\\\",
  "\\endfoot",
  "\\bottomrule",
  "\\endlastfoot"
)
footer <- c("\\end{longtable}", "\\end{footnotesize}")

writeLines(c(header, lines, footer), file.path(OUT, "si_dataset_table.tex"))
cat("wrote si_dataset_table.tex\n")
