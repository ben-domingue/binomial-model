# scratch/diag_imv_new_vs_old.R
#
# Old (poly_compare_multi_compute.R) vs. new (imv branch polytomous-imv,
# R/imv_poly.R) omega_c / omega_t for Binom 2PL -> GPCM in the two Section 4.3
# datasets, on identical predictions. "new_pooled" swaps in the new functions
# with test cells pooled across items, as the manuscript pipeline does;
# "new_byitem" is the unweighted item mean used by imv() for mirt models.
#
# TO REPLICATE
#   1. Get the new IMV code:
#        git clone --branch polytomous-imv https://github.com/ben-domingue/imv.git
#   2. Install R packages irw, imv (CRAN version is fine; only imv.binary is used
#      from it), mirt, dplyr. parallel ships with R.
#   3. Set up IRW/Redivis access so irw::irw_fetch() works. The script pins IRW
#      v252, so everyone gets the same data.
#   4. Make an empty output directory, then run:
#        Rscript diag_imv_new_vs_old.R <path to imv checkout> <outdir>
#   Needs Linux or Mac: the 20 splits use mclapply with mc.cores = 10 (lower it
#   in the last block on a smaller machine). Under 10 minutes on 10 cores.
#
#   Output: prints the seed-20260717 results with per-task weights, then the mean
#   and SD over 20 splits. Also writes paper_split.rds and splits.rds to <outdir>.
#
#   Expect the LSHS seed-20260717 values to match the paper (0.209, 0.111). The
#   Fournier values do not (paper: -0.009, -0.004; here about 0.000); this is a
#   known mismatch within holdout noise (SD about 0.004 for omega_c).
#   The workaround for the integer overflow in .omega_t below can stay after
#   the branch is fixed; it is harmless.

suppressPackageStartupMessages({library(irw); library(imv); library(mirt); library(dplyr); library(parallel)})
args <- commandArgs(TRUE); IMV_DIR <- args[1]; OUTDIR <- args[2]
new <- new.env(parent = globalenv())
sys.source(file.path(IMV_DIR, "R/imv_binary.R"), envir = new)
sys.source(file.path(IMV_DIR, "R/imv_poly.R"), envir = new)
# Branch bug workaround: nlo * (N - nlo) * (...) in .omega_t overflows integers
# for large N. Make the counts double; the weight formula is unchanged.
local({
  orig <- new$.poly_setup
  new$.poly_setup <- function(...) { s <- orig(...); s$Nk <- as.numeric(s$Nk); s }
})
IRW_VERSION <- 252
irw::irw_use_version(version = IRW_VERSION)
TABS <- c(good = "assessment_time_fournier_2025_bscs", poor = "LSHS-E_Quidaja_2022")

binom_P <- function(prob, ncat) {
  size <- ncat - 1L; prob <- pmin(pmax(prob, 1e-10), 1 - 1e-10)
  outer(prob, 0:size, function(p, k) dbinom(k, size, p))
}
binom_2pl <- createItem("binom_2pl", par = c(a = 1, b = 0), est = c(TRUE, TRUE),
  P = function(par, Theta, ncat) binom_P(plogis(par[1] * Theta[, 1] - par[2]), ncat))

old_c <- function(x, pct, P0, P1) {
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
old_t <- function(x, pct, P0, P1) {
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
  byitem <- function(f) mean(sapply(unique(it), function(j) {
    i <- it == j; as.numeric(f(x[i], P0[i, , drop = FALSE], P1[i, , drop = FALSE], K = K))
  }), na.rm = TRUE)
  nc <- new$imv_c(x, P0, P1, K = K); nt <- new$imv_t(x, P0, P1, K = K)
  list(row = data.frame(tab = tab,
         old_c = old_c(x, pct, P0, P1), new_pooled_c = as.numeric(nc), new_byitem_c = byitem(new$imv_c),
         old_t = old_t(x, pct, P0, P1), new_pooled_t = as.numeric(nt), new_byitem_t = byitem(new$imv_t)),
       comp_c = attr(nc, "components"), comp_t = attr(nt, "components"), pct = pct)
}

mats <- lapply(TABS, load_mat)
set.seed(20260717)                    # same seed/order as fig_empirical_illustrative.R
paper <- lapply(names(TABS), function(n) one_split(TABS[[n]], mats[[n]])); names(paper) <- names(TABS)
saveRDS(paper, file.path(OUTDIR, "paper_split.rds"))
print(do.call(rbind, lapply(paper, `[[`, "row"))[, -1], digits = 3)
for (n in names(paper)) { cat("\n", n, "\n"); print(paper[[n]]$comp_c, digits = 3); print(paper[[n]]$comp_t, digits = 3) }

RNGkind("L'Ecuyer-CMRG"); set.seed(1)
reps <- mclapply(1:20, function(s) do.call(rbind, lapply(names(TABS), function(n)
  cbind(split = s, one_split(TABS[[n]], mats[[n]])$row))), mc.cores = 10, mc.preschedule = FALSE)
reps <- do.call(rbind, reps); saveRDS(reps, file.path(OUTDIR, "splits.rds"))
v <- setdiff(names(reps), c("tab", "split"))
cat("\nmean over 20 splits\n"); print(aggregate(reps[v], reps["tab"], mean), digits = 3)
cat("\nSD over 20 splits\n");   print(aggregate(reps[v], reps["tab"], sd), digits = 2)
