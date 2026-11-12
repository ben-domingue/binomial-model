## kl_decomposition.R
## Exact population decomposition of the binomial approximation's log loss
## into a structural part (loss that remains when the binomial's mean equals
## the true expected score) and a link part (extra loss from an incorrect mean).
##
## For a true category distribution P on {0,...,n} at theta, with mean mu,
## m = mu/n, and a binomial B_f = Binom(n, f):
##
##   KL(P || B_f) = KL(P || B_m) + n * KL(Bern(m) || Bern(f)).
##
## KL(P || B_m) is the smallest loss any binomial can achieve at that theta.
## Integrating over theta ~ N(0,1) gives the population decomposition.
##
## Items are drawn from the Study 1 (GPCM) and Study 2 (NP mixture) item
## distributions. For each item we take the Binom 2PL that minimizes
## expected KL (equivalently, by the identity, the link term) and report the
## share of its loss that is structural.
##
## Output: tab_kl_decomp.tex (Overleaf project root), \input by main.tex.
## Run from the repository root:  Rscript kl_decomposition.R

set.seed(20261105)
OUT     <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model/tab_kl_decomp.tex"
N_ITEMS <- 200L                       # items drawn per condition
th      <- seq(-6, 6, length.out = 1201L)
w       <- dnorm(th); w <- w / sum(w) # theta ~ N(0,1)

## ── DGP pieces (copied from sim/gpcm_binomial_sweep*.R) ─────────────────────
make_thresholds <- function(K, asym = 0) {
  n <- K - 1L; total_range <- 3.0
  sym <- seq(-total_range / 2, total_range / 2, length.out = n)
  if (asym == 0 || n <= 2L) return(sym)
  k_idx <- 0:(n - 2L)
  gw <- exp(asym * k_idx / (n - 2L)); gaps <- gw / sum(gw) * total_range
  b <- cumsum(c(-total_range / 2, gaps[-length(gaps)]))
  c(b, b[n - 1L] + gaps[n - 1L])
}
gpcm_probs <- function(theta, a, b) {
  K <- length(b) + 1L; ln <- matrix(0, length(theta), K)
  for (k in seq_len(K - 1L)) ln[, k + 1L] <- ln[, k] + a * (theta - b[k])
  ln <- ln - apply(ln, 1, max); ex <- exp(ln); ex / rowSums(ex)
}
rand_mono_spline <- function(n_knots = 7L) {
  x_k <- seq(-4, 4, length.out = n_knots); y_k <- sort(runif(n_knots))
  splinefun(x_k, y_k, method = "monoH.FC")
}
np_probs <- function(theta, a_j, b_j, K, flex, splines_j) {
  n <- K - 1L; theta_c <- pmin(pmax(theta, -4), 4)
  gc <- gpcm_probs(theta, a_j, b_j)
  gcum <- sapply(seq_len(n), function(k) 1 - rowSums(gc[, seq_len(k), drop = FALSE]))
  rcum <- sapply(seq_len(n), function(k) pmin(pmax(splines_j[[k]](theta_c), 0), 1))
  cum <- (1 - flex) * gcum + flex * rcum
  if (n > 1L) cum <- t(apply(cum, 1, sort, decreasing = TRUE))
  probs <- cbind(1 - cum[, 1],
                 if (n > 1L) cum[, 1:(n - 1L), drop = FALSE] - cum[, 2:n, drop = FALSE],
                 cum[, n])
  pmax(probs, 1e-10)
}

## ── decomposition for one item ─────────────────────────────────────────────
kl_rows <- function(P, Q) rowSums(P * (log(P) - log(pmax(Q, 1e-300))))
bern_kl <- function(m, f) m * log(m / f) + (1 - m) * log((1 - m) / (1 - f))
binom_mat <- function(p, n) outer(p, 0:n, function(p, k) dbinom(k, n, p))

decompose_item <- function(P) {
  P <- P / rowSums(P); n <- ncol(P) - 1L
  m <- pmin(pmax(as.numeric(P %*% (0:n)) / n, 1e-9), 1 - 1e-9)
  structural <- sum(w * kl_rows(P, binom_mat(m, n)))
  ## best Binom 2PL: minimize the link term (equivalent to total KL)
  link_obj <- function(q) {
    f <- pmin(pmax(plogis(exp(q[1]) * (th - q[2])), 1e-12), 1 - 1e-12)
    sum(w * n * bern_kl(m, f))
  }
  start <- c(0, th[which.min(abs(m - 0.5))])
  q <- optim(start, link_obj, method = "BFGS")$par
  f <- pmin(pmax(plogis(exp(q[1]) * (th - q[2])), 1e-12), 1 - 1e-12)
  link  <- sum(w * n * bern_kl(m, f))
  total <- sum(w * kl_rows(P, binom_mat(f, n)))
  stopifnot(abs(total - (structural + link)) < 1e-6 * max(1, total))  # identity check
  c(structural = structural, link = link, total = total)
}

summarise_cond <- function(mat) {
  tot <- mean(mat[, "total"]); str <- mean(mat[, "structural"])
  c(total = tot, share = 100 * str / tot)
}

## ── Study 1: GPCM DGP ──────────────────────────────────────────────────────
s_vals <- c(-log(8), -1, 0, 1, log(8))
study1 <- do.call(rbind, lapply(4:7, function(K) do.call(rbind, lapply(s_vals, function(s) {
  b0 <- make_thresholds(K, s)
  res <- t(replicate(N_ITEMS, {
    a <- rlnorm(1, 0, 0.5); b <- b0 + runif(1, -1, 1) + rnorm(K - 1L, 0, 0.3)
    decompose_item(gpcm_probs(th, a, b))
  }))
  data.frame(K = K, s = s, t(summarise_cond(res)))
}))))

## ── Study 2: NP DGP ────────────────────────────────────────────────────────
flex_vals <- c(0, 0.25, 0.5, 0.75, 1)
study2 <- do.call(rbind, lapply(c(5L, 7L), function(K) do.call(rbind, lapply(flex_vals, function(fl) {
  res <- t(replicate(N_ITEMS, {
    a <- rlnorm(1, 0, 0.5)
    b <- seq(-1.5, 1.5, length.out = K - 1L) + runif(1, -1, 1) + rnorm(K - 1L, 0, 0.3)
    spl <- lapply(seq_len(K - 1L), function(k) rand_mono_spline())
    decompose_item(np_probs(th, a, b, K, fl, spl))
  }))
  data.frame(K = K, flex = fl, t(summarise_cond(res)))
}))))

print(study1, digits = 3); print(study2, digits = 3)
message(sprintf("Study 1 structural share: range %.0f-%.0f%%; Study 2: %.0f-%.0f%%",
                min(study1$share), max(study1$share), min(study2$share), max(study2$share)))

## ── LaTeX table ────────────────────────────────────────────────────────────
cell <- function(sh, tot) sprintf("%.1f\\%% (%.3f)", sh, tot)
s_lab <- c("$-\\ln 8$", "$-1$", "$0$", "$1$", "$\\ln 8$")
rows1 <- sapply(4:7, function(K) {
  d <- study1[study1$K == K, ]
  paste0("$K=", K, "$ & ", paste(cell(d$share, d$total), collapse = " & "), " \\\\")
})
rows2 <- sapply(c(5L, 7L), function(K) {
  d <- study2[study2$K == K, ]
  paste0("$K=", K, "$ & ", paste(cell(d$share, d$total), collapse = " & "), " \\\\")
})
tex <- c(
  "% Generated by kl_decomposition.R -- do not edit by hand.",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Exact decomposition of the binomial approximation's excess expected log loss, its KL divergence from the true category distribution (Equation~\\ref{eq:kl-decomp}). Each cell gives the structural share---the percentage of the best-fitting Binom 2PL's loss that remains when the binomial's mean equals the true expected score---followed in parentheses by that model's total expected loss, $\\E_\\theta[D_{\\mathrm{KL}}(P_\\theta \\Vert B_{f(\\theta)})]$, in nats per response. Population values with $\\theta \\sim \\mathcal{N}(0,1)$, averaged over ",
  paste0("  ", N_ITEMS, " items drawn from each study's item distribution; no data are simulated and no models are estimated.}"),
  "\\label{tab:kl-decomp}",
  "\\small",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "\\multicolumn{6}{l}{\\textit{Study 1 (GPCM DGP), by gap asymmetry $s$}} \\\\",
  paste0(" & ", paste(s_lab, collapse = " & "), " \\\\"),
  "\\midrule", rows1,
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Study 2 (NP DGP), by flex $\\delta$}} \\\\",
  paste0(" & ", paste(sprintf("$%s$", c("0", "0.25", "0.5", "0.75", "1")), collapse = " & "), " \\\\"),
  "\\midrule", rows2,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}")
writeLines(tex, OUT)
message("wrote ", OUT)
