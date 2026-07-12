## diag_imvt_decomp.R
## Decomposes ω_t oracle advantage for the NP spline by threshold k, at a range
## of flex values.  Goal: locate where the non-monotone ω_t pattern comes from.
##
## For each (flex, rep): compute per-threshold binary IMV(NP spline, oracle)
## and the marginal P(X > k) from the oracle, to see if the log-score gap
## saturates or changes character at high flex.
##
## Output: src/diag_imvt_decomp.pdf

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scam); library(imv)
})

OUT <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model/src"

## ── helpers copied from gpcm_binomial_sweep_np.R ─────────────────────────────

gpcm_probs <- function(theta, a, b) {
  K  <- length(b) + 1L
  ln <- matrix(0, nrow = length(theta), ncol = K)
  for (k in seq_len(K - 1L))
    ln[, k + 1L] <- ln[, k] + a * (theta - b[k])
  ln <- ln - apply(ln, 1, max)
  ex <- exp(ln); ex / rowSums(ex)
}

rand_mono_spline <- function(n_knots = 7L) {
  x_k <- seq(-4, 4, length.out = n_knots)
  y_k <- sort(runif(n_knots))
  splinefun(x_k, y_k, method = "monoH.FC")
}

np_probs <- function(theta, a_j, b_j, K, flex, splines_j) {
  n       <- K - 1L; theta_c <- pmin(pmax(theta, -4), 4)
  gpcm_cat <- gpcm_probs(theta, a_j, b_j)
  gpcm_cum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n))
    gpcm_cum[, k] <- 1 - rowSums(gpcm_cat[, seq_len(k), drop = FALSE])
  rand_cum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n))
    rand_cum[, k] <- pmin(pmax(splines_j[[k]](theta_c), 0), 1)
  cum <- (1 - flex) * gpcm_cum + flex * rand_cum
  if (n > 1L) cum <- t(apply(cum, 1, sort, decreasing = TRUE))
  probs <- cbind(1 - cum[, 1],
                 if (n > 1L) cum[, seq_len(n-1L), drop=FALSE] -
                              cum[, seq_len(n-1L)+1L, drop=FALSE],
                 cum[, n])
  pmax(probs, 1e-10)
}

## per-threshold binary IMV: positive = oracle better than model
imv_t_by_k <- function(pr_model, pr_oracle, X, K) {
  sapply(0:(K - 2L), function(ii) {
    resp <- as.integer(X <= ii)
    p1   <- rowSums(pr_model [, 1:(ii + 1L), drop = FALSE])
    p2   <- rowSums(pr_oracle[, 1:(ii + 1L), drop = FALSE])
    tryCatch(imv.binary(resp, p1, p2), error = function(e) NA_real_)
  })
}

## mean oracle P(X <= k) and NP spline P(X <= k) across holdout observations
mean_cum_probs <- function(pr, K)
  colMeans(apply(pr, 1, function(p) cumsum(p)[1:(K-1)]))

## ── mini simulation ───────────────────────────────────────────────────────────

set.seed(42)
K        <- 5L
N        <- 1000L
N_ITEMS  <- 10L
N_REPS   <- 30L
flex_vals <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)

results <- do.call(rbind, lapply(flex_vals, function(flex) {
  do.call(rbind, lapply(seq_len(N_REPS), function(rep_idx) {
    set.seed(rep_idx * 997L + as.integer(round(flex * 1000)) + 10000L)

    theta    <- rnorm(N)
    n_bound  <- K - 1L
    X_mat    <- matrix(NA_integer_, N, N_ITEMS)
    item_par <- vector("list", N_ITEMS)

    for (j in seq_len(N_ITEMS)) {
      a_j <- rlnorm(1, 0, 0.5)
      b_j <- seq(-1.5, 1.5, length.out = n_bound) + rnorm(n_bound, 0, 0.3)
      spl <- lapply(seq_len(n_bound), function(k) rand_mono_spline())
      pr  <- np_probs(theta, a_j, b_j, K, flex, spl)
      pr  <- pr / rowSums(pr)
      X_mat[, j] <- apply(pr, 1, function(p) sample.int(K, 1L, prob = p) - 1L)
      item_par[[j]] <- list(a = a_j, b = b_j, flex = flex, splines = spl)
    }
    if (any(sapply(seq_len(N_ITEMS), function(j) length(unique(X_mat[, j])) < K)))
      return(NULL)

    ## 20% holdout
    n_ho    <- as.integer(0.2 * N * N_ITEMS)
    ho_lin  <- sample(N * N_ITEMS, n_ho)
    ho_mat  <- matrix(FALSE, N, N_ITEMS); ho_mat[ho_lin] <- TRUE
    ho_row  <- row(ho_mat)[ho_lin]; ho_col <- col(ho_mat)[ho_lin]
    X_ho    <- X_mat[ho_lin]
    pctt    <- as.numeric(table(factor(X_ho, levels = 0L:(K-1L)))) / n_ho

    X_tr_df <- as.data.frame(X_mat); X_tr_df[ho_mat] <- NA_integer_

    ## oracle predictions (true theta + true params)
    pr_oracle <- matrix(NA_real_, n_ho, K)
    for (j in seq_len(N_ITEMS)) {
      cells_j <- which(ho_col == j)
      if (!length(cells_j)) next
      tp <- item_par[[j]]
      pr_oracle[cells_j, ] <- np_probs(theta[ho_row[cells_j]], tp$a, tp$b, K,
                                        tp$flex, tp$splines)
    }

    ## NP spline: fit per item on training, predict on holdout using true theta
    pr_np <- matrix(NA_real_, n_ho, K)
    ok <- TRUE
    for (j in seq_len(N_ITEMS)) {
      obs_j <- which(!ho_mat[, j])
      df_j  <- data.frame(X = X_mat[obs_j, j], th = theta[obs_j])
      if (length(unique(df_j$X)) < K) { ok <- FALSE; break }
      fit_j <- tryCatch(
        scam(cbind(X, K - 1L - X) ~ s(th, bs = "mpi"),
             data = df_j, family = binomial()),
        error = function(e) NULL
      )
      if (is.null(fit_j)) { ok <- FALSE; break }
      cells_j <- which(ho_col == j)
      if (!length(cells_j)) next
      p_hat <- pmin(pmax(
        as.numeric(predict(fit_j, newdata = data.frame(th = theta[ho_row[cells_j]]),
                           type = "response")), 1e-10), 1 - 1e-10)
      pr_np[cells_j, ] <- outer(p_hat, 0:(K-1L), function(p, k) dbinom(k, K-1L, p))
    }
    if (!ok || anyNA(pr_np) || anyNA(pr_oracle)) return(NULL)

    ## per-threshold IMV (positive = oracle beats NP spline)
    imvt_k <- imv_t_by_k(pr_np, pr_oracle, X_ho, K)

    ## mean cumulative probs from oracle and NP spline across holdout obs
    oracle_cum <- mean_cum_probs(pr_oracle, K)
    np_cum     <- mean_cum_probs(pr_np, K)

    data.frame(
      flex     = flex,
      rep      = rep_idx,
      threshold = paste0("k<=", 0:(K-2L)),
      imvt_k   = imvt_k,
      oracle_cum = oracle_cum,
      np_cum     = np_cum,
      wt         = (pctt[1:(K-1L)] / (1 - pctt[K]))[1:(K-1L)]
    )
  }))
}))

## ── aggregate and plot ────────────────────────────────────────────────────────

agg <- results |>
  group_by(flex, threshold) |>
  summarize(
    mean_imvt   = mean(imvt_k,    na.rm = TRUE),
    mean_oracle = mean(oracle_cum, na.rm = TRUE),
    mean_np     = mean(np_cum,     na.rm = TRUE),
    mean_wt     = mean(wt,         na.rm = TRUE),
    .groups = "drop"
  )

## weighted ω_t (should reproduce main results)
wt_imvt <- agg |>
  group_by(flex) |>
  summarize(
    omega_t = sum(mean_imvt * mean_wt) / sum(mean_wt),
    .groups = "drop"
  )
cat("Weighted ω_t by flex (should peak ~0.45 then decrease):\n")
print(wt_imvt)

## Panel A: per-threshold IMV by flex
pA <- ggplot(agg, aes(x = flex, y = mean_imvt, colour = threshold, group = threshold)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0,
           linetype = "dashed", colour = "grey50") +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "A: per-threshold binary IMV (oracle vs NP spline)",
       subtitle = "positive = oracle wins at that threshold",
       x = "flex", y = "IMV(NP spline, oracle)") +
  theme_bw(base_size = 11)

## Panel B: mean oracle P(X <= k) by flex
pB <- ggplot(agg, aes(x = flex, y = mean_oracle, colour = threshold, group = threshold)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "B: mean oracle cumulative probability P(X ≤ k)",
       subtitle = "approaches 0.5 as flex -> 1 (random splines are uniform on [0,1])",
       x = "flex", y = "mean P(X ≤ k) from oracle") +
  theme_bw(base_size = 11)

## Panel C: mean NP spline P(X <= k) by flex
pC <- ggplot(agg, aes(x = flex, y = mean_np, colour = threshold, group = threshold)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "C: mean NP spline cumulative probability P(X ≤ k)",
       x = "flex", y = "mean P(X ≤ k) from NP spline") +
  theme_bw(base_size = 11)

## Panel D: |oracle - NP spline| by threshold and flex
pD <- agg |>
  mutate(gap = abs(mean_oracle - mean_np)) |>
  ggplot(aes(x = flex, y = gap, colour = threshold, group = threshold)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "D: |oracle mean P(X ≤ k) - NP spline mean P(X ≤ k)|",
       subtitle = "shrinking gap = oracle and NP spline converging in mean prediction",
       x = "flex", y = "|oracle - NP spline| mean cum. prob") +
  theme_bw(base_size = 11)

library(patchwork)
p <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "ω_t decomposition: NP spline vs oracle by threshold and flex",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )

ggsave(file.path(OUT, "diag_imvt_decomp.pdf"), plot = p,
       width = 14, height = 9, device = cairo_pdf)
message("Saved src/diag_imvt_decomp.pdf")
