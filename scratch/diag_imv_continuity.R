# scratch/diag_imv_continuity.R
#
# Section 4.3 datasets (Binom 2PL -> GPCM) scored with the imv v0.4.0 polytomous
# IMVs (klintkanopka/imv, branch poly-weights-aggregate), on identical predictions:
#   paper      the manuscript's own imv_c/imv_t (poly_compare_multi_compute.R),
#              test cells pooled, weights from training proportions
#   legacy     imv_c/imv_t(weights = "legacy"), pooled (check: should ~ paper;
#              differs only in taking weights from held-out rather than training counts)
#   cont_pool  weights = "continuity", pooled
#   cont_item  weights = "continuity", one IMV per item, mean weighted by held-out
#              responses (package default for polytomous mirt: aggregate = "items",
#              item.weights = "responses")  <- adopted 2026-11-11
#   cont_eq    as cont_item with equal item weights
#
# Usage: Rscript diag_imv_continuity.R <imv checkout on poly-weights-aggregate> <outdir>
# Needs irw/Redivis access (IRW v252 pinned), mirt, dplyr; Linux/Mac for mclapply.

suppressPackageStartupMessages({library(irw); library(mirt); library(dplyr); library(parallel)})
args <- commandArgs(TRUE); IMV_DIR <- args[1]; OUTDIR <- args[2]
new <- new.env(parent = globalenv())
for (f in c("imv.R", "imv_binary.R", "imv_poly.R")) sys.source(file.path(IMV_DIR, "R", f), envir = new)
imv.binary <- new$imv.binary
irw::irw_use_version(version = 252)
TABS <- c(good = "assessment_time_fournier_2025_bscs", poor = "LSHS-E_Quidaja_2022")

binom_P <- function(prob, ncat) {
  size <- ncat - 1L; prob <- pmin(pmax(prob, 1e-10), 1 - 1e-10)
  outer(prob, 0:size, function(p, k) dbinom(k, size, p))
}
binom_2pl <- createItem("binom_2pl", par = c(a = 1, b = 0), est = c(TRUE, TRUE),
  P = function(par, Theta, ncat) binom_P(plogis(par[1] * Theta[, 1] - par[2]), ncat))

paper_c <- function(x, pct, P0, P1) {
  nn <- length(pct); om <- numeric(nn); iis <- 0:(nn - 1)
  for (ii in iis) {
    ns <- o <- numeric()
    for (jj in iis[-match(ii, iis)]) {
      i <- x %in% c(ii, jj)
      q <- function(P) P[i, ii + 1] / (P[i, ii + 1] + P[i, jj + 1])
      o[as.character(jj)]  <- imv.binary(as.numeric(x[i] == ii), q(P0), q(P1))
      ns[as.character(jj)] <- sum(i)
    }
    om[ii + 1] <- sum(o * ns) / sum(ns)
  }
  sum(om * pct) / sum(pct)
}
paper_t <- function(x, pct, P0, P1) {
  nn <- length(pct); om <- numeric(nn - 1)
  for (ii in 0:(nn - 2)) {
    om[ii + 1] <- imv.binary(as.numeric(x <= ii),
                             rowSums(P0[, 1:(ii + 1), drop = FALSE]),
                             rowSums(P1[, 1:(ii + 1), drop = FALSE]))
  }
  w <- pct[1:(nn - 1)] / (1 - pct[nn])
  sum(om * w) / sum(w)
}

load_mat <- function(tab) {
  df <- irw::irw_fetch(tab)
  df$resp0 <- df$resp - min(df$resp, na.rm = TRUE)
  w <- irw_long2resp(df %>% select(-resp) %>% rename(resp = resp0))
  m <- as.matrix(w[, setdiff(names(w), "id")])
  m <- m[, apply(m, 2, function(x) length(unique(na.omit(x))) >= 2), drop = FALSE]
  iK <- apply(m, 2, function(x) length(unique(na.omit(x))))
  m[, iK == as.integer(names(sort(table(iK), decreasing = TRUE))[1]), drop = FALSE]
}

q <- function(expr) suppressWarnings(as.numeric(expr))   # fair-coin floor warnings

one_split <- function(tab, m) {
  ni <- ncol(m); K <- length(unique(na.omit(as.vector(m))))
  nm <- which(!is.na(m), arr.ind = TRUE)
  h  <- sample(nrow(nm), floor(0.2 * nrow(nm)))
  tr <- m; tr[nm[h, , drop = FALSE]] <- NA
  tc <- nm[h, , drop = FALSE]; x <- m[tc]
  fit <- function(type, ...) mirt(tr, 1, itemtype = rep(type, ni), verbose = FALSE,
                                  technical = list(NCYCLES = 2000), ...)
  probs <- function(mod) {
    th <- fscores(mod, response.pattern = tr, method = "EAP")[, 1]
    P <- matrix(NA_real_, nrow(tc), K)
    for (j in seq_len(ni)) {
      r <- which(tc[, 2] == j); if (!length(r)) next
      pr <- probtrace(extract.item(mod, j), matrix(th[tc[r, 1]], ncol = 1))
      if (ncol(pr) == K) P[r, ] <- pr
    }
    P
  }
  P1 <- probs(fit("gpcm"))
  P0 <- probs(fit("binom_2pl", customItems = list(binom_2pl = binom_2pl)))
  ok <- complete.cases(P0, P1); x <- x[ok]; P0 <- P0[ok, ]; P1 <- P1[ok, ]; it <- tc[ok, 2]
  to <- as.integer(tr[!is.na(tr)])
  pct <- as.numeric(table(factor(to, levels = 0:(K - 1)))) / length(to)
  byitem <- function(f, equal = FALSE) {
    js <- unique(it)
    v <- sapply(js, function(j) { i <- it == j; q(f(x[i], P0[i, , drop = FALSE], P1[i, , drop = FALSE], K = K)) })
    n <- sapply(js, function(j) sum(it == j)); if (equal) n[] <- 1
    ok <- !is.na(v); sum(n[ok] * v[ok]) / sum(n[ok])
  }
  cc <- function(...) new$imv_c(..., weights = "continuity")
  ct <- function(...) new$imv_t(..., weights = "continuity")
  data.frame(tab = names(TABS)[TABS == tab],
    paper_c = paper_c(x, pct, P0, P1),
    legacy_c = q(new$imv_c(x, P0, P1, K = K, weights = "legacy")),
    cont_pool_c = q(cc(x, P0, P1, K = K)),
    cont_item_c = byitem(cc), cont_eq_c = byitem(cc, TRUE),
    paper_t = paper_t(x, pct, P0, P1),
    legacy_t = q(new$imv_t(x, P0, P1, K = K, weights = "legacy")),
    cont_pool_t = q(ct(x, P0, P1, K = K)),
    cont_item_t = byitem(ct), cont_eq_t = byitem(ct, TRUE))
}

mats <- lapply(TABS, load_mat)
suppressWarnings(new$.warn_legacy_weights())
set.seed(20260717)                    # same seed/order as fig_empirical_illustrative.R
paper <- do.call(rbind, lapply(TABS, function(t) one_split(t, mats[[names(TABS)[TABS == t]]])))
saveRDS(paper, file.path(OUTDIR, "cont_paper_split.rds"))
cat("seed 20260717\n"); print(paper, digits = 3, row.names = FALSE)

RNGkind("L'Ecuyer-CMRG"); set.seed(1)
reps <- mclapply(1:20, function(s) do.call(rbind, lapply(TABS, function(t)
  cbind(split = s, one_split(t, mats[[names(TABS)[TABS == t]]])))), mc.cores = 10, mc.preschedule = FALSE)
reps <- do.call(rbind, reps); saveRDS(reps, file.path(OUTDIR, "cont_splits.rds"))
v <- setdiff(names(reps), c("tab", "split"))
cat("\nmean over 20 splits\n"); print(aggregate(reps[v], reps["tab"], mean), digits = 3)
cat("\nSD over 20 splits\n");   print(aggregate(reps[v], reps["tab"], sd), digits = 2)
