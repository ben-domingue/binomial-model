# src/diag_grid_resolution.R
#
# Tests whether the 101-point theta grid causes artifical IMV identity
# across binomial models at high flex.
#
# Approach: take 10 reps from the existing RDS (5 high-flex, 5 low-flex),
# re-run the theta-scoring step for 1PL-logit vs 2PL-logit with 101 and
# 501 grid points, and compare theta_est and holdout IMV.
#
# Key question: at high flex, does a finer grid break the exact identity,
# or do the models genuinely converge?

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scam)
  library(imv)
})

SIM <- "/home/ben/Dropbox/projects/irw/irw_site/sims"

# ── Shared model helpers ──────────────────────────────────────────────────────

gpcm_probs <- function(theta, a, b) {
  K  <- length(b) + 1L
  ln <- matrix(0, nrow = length(theta), ncol = K)
  for (k in seq_len(K - 1L))
    ln[, k + 1L] <- ln[, k] + a * (theta - b[k])
  ln <- ln - apply(ln, 1, max)
  ex <- exp(ln); ex / rowSums(ex)
}

binom_probs_fn <- function(theta, a, b, K) {
  p <- pmin(pmax(plogis(a * theta - b), 1e-10), 1 - 1e-10)
  outer(p, 0:(K - 1L), function(pp, k) dbinom(k, K - 1L, pp))
}

item_ll <- function(probs, X)
  sum(log(pmax(probs[cbind(seq_along(X), X + 1L)], 1e-15)))

fit_binom_1pl <- function(X, theta) {
  K <- max(X) + 1L
  mod <- glm(cbind(X, K - 1L - X) ~ offset(theta), family = binomial("logit"))
  list(a = 1.0, b = -coef(mod)[["(Intercept)"]])
}

fit_binom_2pl <- function(X, theta) {
  mod <- glm(cbind(X, max(X) + 1L - 1L - X) ~ theta, family = binomial("logit"))
  list(a = coef(mod)[["theta"]], b = -coef(mod)[["(Intercept)"]])
}

rand_mono_spline <- function(n_knots = 7L) {
  x_k <- seq(-4, 4, length.out = n_knots)
  y_k <- sort(runif(n_knots))
  splinefun(x_k, y_k, method = "monoH.FC")
}

np_probs_fn <- function(theta, a_j, b_j, K, flex, splines_j) {
  n       <- K - 1L
  theta_c <- pmin(pmax(theta, -4), 4)
  gpcm_cat <- gpcm_probs(theta, a_j, b_j)
  gpcm_cum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n))
    gpcm_cum[, k] <- 1 - rowSums(gpcm_cat[, seq_len(k), drop = FALSE])
  rand_cum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n))
    rand_cum[, k] <- pmin(pmax(splines_j[[k]](theta_c), 0), 1)
  cum <- (1 - flex) * gpcm_cum + flex * rand_cum
  if (n > 1L) cum <- t(apply(cum, 1, sort, decreasing = TRUE))
  probs <- cbind(
    1 - cum[, 1],
    if (n > 1L) cum[, seq_len(n - 1L), drop = FALSE] -
                cum[, seq_len(n - 1L) + 1L, drop = FALSE],
    cum[, n]
  )
  pmax(probs, 1e-10)
}

imv_c_fn <- function(y, pctt.tab, p1, p2) {
  nn  <- length(pctt.tab); iis <- 0:(nn - 1L); om <- numeric(nn)
  for (ii in iis) {
    ns <- om.tmp <- numeric()
    for (jj in iis[-match(ii, iis)]) {
      y2   <- y[y$resp %in% c(ii, jj), ]
      resp <- ifelse(y2$resp == ii, 1, 0)
      p1ii <- y2[[paste0(p1, ii)]]; p1jj <- y2[[paste0(p1, jj)]]
      p2ii <- y2[[paste0(p2, ii)]]; p2jj <- y2[[paste0(p2, jj)]]
      z    <- data.frame(resp = resp,
                         p1   = p1ii / (p1ii + p1jj),
                         p2   = p2ii / (p2ii + p2jj))
      om.tmp[as.character(jj)] <- imv.binary(z$resp, z$p1, z$p2)
      ns[as.character(jj)]     <- nrow(z)
    }
    om[ii + 1L] <- sum(om.tmp * ns) / sum(ns)
  }
  sum(om * pctt.tab) / sum(pctt.tab)
}

# ── Score theta from grid for one model ───────────────────────────────────────
score_theta <- function(fits, predict_fn, X_mat, ho_mat, N, K, n_grid) {
  theta_grid  <- seq(-4, 4, length.out = n_grid)
  N_ITEMS     <- ncol(X_mat)
  total_ll    <- matrix(0, N, length(theta_grid))
  for (j in seq_len(N_ITEMS)) {
    obs_j   <- which(!ho_mat[, j])
    pr_grid <- predict_fn(fits[[j]], theta_grid, K)
    total_ll[obs_j, ] <- total_ll[obs_j, ] +
      t(log(pmax(pr_grid, 1e-15)))[X_mat[obs_j, j] + 1L, ]
  }
  theta_grid[apply(total_ll, 1, which.max)]
}

# ── One simulated rep ─────────────────────────────────────────────────────────
run_grid_check <- function(flex, rep_seed, n_grids = c(101L, 501L),
                           N = 1000L, N_ITEMS = 10L, K = 5L) {
  set.seed(rep_seed)
  theta_true <- rnorm(N)
  shifts     <- runif(N_ITEMS, -1, 1)

  X_mat       <- matrix(NA_integer_, N, N_ITEMS)
  item_params <- vector("list", N_ITEMS)
  for (j in seq_len(N_ITEMS)) {
    a_j       <- rlnorm(1L, meanlog = 0, sdlog = 0.5)
    b_j       <- seq(-1.5, 1.5, length.out = K - 1L) + shifts[j] +
                 rnorm(K - 1L, 0, 0.3)
    splines_j <- lapply(seq_len(K - 1L), function(k) rand_mono_spline())
    probs_j   <- np_probs_fn(theta_true, a_j, b_j, K, flex, splines_j)
    probs_j   <- probs_j / rowSums(probs_j)
    X_j       <- apply(probs_j, 1, function(p) sample.int(K, 1L, prob = p) - 1L)
    tries <- 0L
    while (length(unique(X_j)) < K && tries < 5L) {
      shifts[j]  <- runif(1, -1, 1)
      b_j        <- seq(-1.5, 1.5, length.out = K - 1L) + shifts[j] + rnorm(K-1L,0,0.3)
      splines_j  <- lapply(seq_len(K - 1L), function(k) rand_mono_spline())
      probs_j    <- np_probs_fn(theta_true, a_j, b_j, K, flex, splines_j)
      probs_j    <- probs_j / rowSums(probs_j)
      X_j        <- apply(probs_j, 1, function(p) sample.int(K, 1L, prob = p) - 1L)
      tries      <- tries + 1L
    }
    X_mat[, j]       <- X_j
    item_params[[j]] <- list(a = a_j, b = b_j, splines = splines_j)
  }

  # 20% cell holdout
  n_ho   <- as.integer(0.2 * N * N_ITEMS)
  ho_lin <- sample(N * N_ITEMS, n_ho)
  ho_mat <- matrix(FALSE, N, N_ITEMS)
  ho_mat[ho_lin] <- TRUE
  ho_row <- row(ho_mat)[ho_lin]
  ho_col <- col(ho_mat)[ho_lin]
  all_X_ho <- X_mat[ho_lin]
  pctt.tab <- as.numeric(table(factor(all_X_ho, 0:(K-1)))) / length(all_X_ho)

  # Fit 1PL and 2PL using true theta (on training cells only)
  fits_1pl <- lapply(seq_len(N_ITEMS), function(j) {
    obs_j <- !ho_mat[, j]
    fit_binom_1pl(X_mat[obs_j, j], theta_true[obs_j])
  })
  fits_2pl <- lapply(seq_len(N_ITEMS), function(j) {
    obs_j <- !ho_mat[, j]
    fit_binom_2pl(X_mat[obs_j, j], theta_true[obs_j])
  })

  pred_1pl <- function(f, th, K) binom_probs_fn(th, f$a, f$b, K)
  pred_2pl <- function(f, th, K) binom_probs_fn(th, f$a, f$b, K)

  results <- lapply(n_grids, function(ng) {
    th_1pl <- score_theta(fits_1pl, pred_1pl, X_mat, ho_mat, N, K, ng)
    th_2pl <- score_theta(fits_2pl, pred_2pl, X_mat, ho_mat, N, K, ng)

    # Fraction of persons where theta_est is identical
    frac_same <- mean(th_1pl == th_2pl)

    # Holdout predictions
    pr_1pl <- matrix(NA_real_, n_ho, K)
    pr_2pl <- matrix(NA_real_, n_ho, K)
    for (j in seq_len(N_ITEMS)) {
      cells_j <- which(ho_col == j)
      if (length(cells_j) == 0L) next
      pr_1pl[cells_j, ] <- pred_1pl(fits_1pl[[j]], th_1pl[ho_row[cells_j]], K)
      pr_2pl[cells_j, ] <- pred_2pl(fits_2pl[[j]], th_2pl[ho_row[cells_j]], K)
    }

    # IMV: 2PL vs 1PL (positive = 2PL wins)
    df <- as.data.frame(pr_1pl); colnames(df) <- paste0("p1", 0:(K-1))
    df$resp <- all_X_ho
    for (k in 0:(K-1)) df[[paste0("p2", k)]] <- pr_2pl[, k+1L]
    imv_val <- tryCatch(imv_c_fn(df, pctt.tab, "p1", "p2"), error = function(e) NA_real_)

    data.frame(flex = flex, seed = rep_seed, n_grid = ng,
               frac_theta_same = frac_same,
               imv_2pl_vs_1pl  = imv_val,
               a_1pl_mean = mean(sapply(fits_1pl, `[[`, "a")),
               a_2pl_mean = mean(sapply(fits_2pl, `[[`, "a")))
  })
  do.call(rbind, results)
}

# ── Run checks ────────────────────────────────────────────────────────────────
cat("Running grid resolution check (low flex = 0.3, high flex = 0.85)...\n")

low_flex_seeds  <- 101:105
high_flex_seeds <- 201:205

low_results <- do.call(rbind, lapply(low_flex_seeds,  function(s) run_grid_check(0.3,  s)))
high_results <- do.call(rbind, lapply(high_flex_seeds, function(s) run_grid_check(0.85, s)))

all_results <- bind_rows(low_results, high_results)

cat("\n=== Results: fraction of persons with identical theta_est (1PL vs 2PL) ===\n")
print(
  all_results |>
    group_by(flex, n_grid) |>
    summarise(frac_same_mean = mean(frac_theta_same),
              imv_mean       = mean(imv_2pl_vs_1pl, na.rm=TRUE),
              a_2pl_mean     = mean(a_2pl_mean),
              .groups = "drop"),
  digits = 4
)

cat("\n=== Per-rep detail ===\n")
print(all_results, digits = 4)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(all_results |> mutate(flex_lab = paste0("flex = ", flex)),
            aes(x = factor(n_grid), y = frac_theta_same, colour = factor(flex_lab))) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_colour_manual(values = c("flex = 0.3" = "#2166ac", "flex = 0.85" = "#d6604d"),
                      name = NULL) +
  labs(x = "Theta grid size (points)", y = "Fraction of persons: theta_est identical (1PL = 2PL)",
       title = "Grid resolution check: does a finer grid break the binomial identity?") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("src/diag_grid_resolution.pdf", p, width = 6, height = 4)
message("Saved src/diag_grid_resolution.pdf")
