# src/fig_single_item_example.R
#
# Single-item figure for the NP DGP section.
# Picks the item with the highest True-vs-1PL CRF divergence
# from the 20-item flex=0.3 dataset (same seed as fig_np_crf_example.R).
#
# Layout: 2 rows x 3 columns
#   Row 1: ERF  |  k=0  |  k=1
#   Row 2: k=2  |  k=3  |  k=4
#
# Each panel overlays True (black solid), GPCM (blue dashed), 1PL (red dotted).
#
# Outputs:
#   fig_single_item_example.pdf (Overleaf project root)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(mirt)
  library(imv)
})

# ── Parameters (must match fig_np_crf_example.R) ─────────────────────────────
set.seed(123)
N       <- 2000L
N_ITEMS <- 10L
K       <- 5L
FLEX    <- 0.3
theta_g <- seq(-3.5, 3.5, length.out = 301L)

# ── Helpers ───────────────────────────────────────────────────────────────────
gpcm_probs <- function(theta, a, b) {
  Kk <- length(b) + 1L
  ln <- matrix(0, length(theta), Kk)
  for (k in seq_len(Kk - 1L)) ln[, k+1L] <- ln[, k] + a*(theta - b[k])
  ln <- ln - apply(ln, 1, max)
  ex <- exp(ln); ex / rowSums(ex)
}
binom_probs <- function(theta, a, b, K) {
  p <- pmin(pmax(plogis(a*theta - b), 1e-10), 1-1e-10)
  outer(p, 0:(K-1L), function(pp, k) dbinom(k, K-1L, pp))
}
rand_mono_spline <- function(n_knots = 7L) {
  x_k <- seq(-4, 4, length.out = n_knots)
  y_k <- sort(runif(n_knots))
  splinefun(x_k, y_k, method = "monoH.FC")
}
np_probs <- function(theta, a_j, b_j, K, flex, splines_j) {
  n <- K - 1L; theta_c <- pmin(pmax(theta, -4), 4)
  gc <- gpcm_probs(theta, a_j, b_j)
  gcum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n)) gcum[,k] <- 1 - rowSums(gc[, seq_len(k), drop=FALSE])
  rcum <- matrix(NA_real_, length(theta), n)
  for (k in seq_len(n)) rcum[,k] <- pmin(pmax(splines_j[[k]](theta_c), 0), 1)
  cum <- (1-flex)*gcum + flex*rcum
  if (n > 1L) cum <- t(apply(cum, 1, sort, decreasing=TRUE))
  probs <- cbind(1-cum[,1],
                 if(n>1L) cum[,seq_len(n-1L),drop=FALSE]-cum[,seq_len(n-1L)+1L,drop=FALSE],
                 cum[,n])
  pmax(probs, 1e-10)
}
erf_fn <- function(probs, K) as.numeric(probs %*% 0:(K-1L)) / (K - 1L)

# ── Generate data ─────────────────────────────────────────────────────────────
cat("Generating data...\n")
theta_true <- rnorm(N)
shifts     <- runif(N_ITEMS, -1, 1)

true_params <- vector("list", N_ITEMS)
X_mat       <- matrix(NA_integer_, N, N_ITEMS)

for (j in seq_len(N_ITEMS)) {
  a_j <- rlnorm(1L, meanlog=0, sdlog=0.5)
  b_j <- seq(-1.5, 1.5, length.out=K-1L) + shifts[j] + rnorm(K-1L, 0, 0.3)
  spl <- lapply(seq_len(K-1L), function(k) rand_mono_spline())
  pr  <- np_probs(theta_true, a_j, b_j, K, FLEX, spl)
  pr  <- pr / rowSums(pr)
  X_j <- apply(pr, 1, function(p) sample.int(K, 1L, prob=p) - 1L)
  tries <- 0L
  while (length(unique(X_j)) < K && tries < 5L) {
    shifts[j] <- runif(1, -1, 1)
    b_j <- seq(-1.5,1.5,length.out=K-1L) + shifts[j] + rnorm(K-1L,0,0.3)
    spl <- lapply(seq_len(K-1L), function(k) rand_mono_spline())
    pr  <- np_probs(theta_true, a_j, b_j, K, FLEX, spl)
    pr  <- pr / rowSums(pr); X_j <- apply(pr, 1, function(p) sample.int(K,1L,prob=p)-1L)
    tries <- tries + 1L
  }
  true_params[[j]] <- list(a=a_j, b=b_j, splines=spl)
  X_mat[,j]        <- X_j
}

# ── Fit models ────────────────────────────────────────────────────────────────
cat("Fitting GPCM...\n")
mod_gpcm  <- suppressMessages(mirt(as.data.frame(X_mat), 1, itemtype="gpcm", verbose=FALSE))
theta_eap <- as.numeric(fscores(mod_gpcm, method="EAP")[,1])

gpcm_params <- lapply(seq_len(N_ITEMS), function(j) {
  p <- coef(mod_gpcm, IRTpars=TRUE)[[j]]
  list(a=as.numeric(p[1,"a"]), b=as.numeric(p[1, paste0("b",seq_len(K-1L))]))
})
binom_params <- lapply(seq_len(N_ITEMS), function(j) {
  df  <- data.frame(X=X_mat[,j], theta=theta_eap)
  mod <- glm(cbind(X, K-1L-X) ~ offset(theta), data=df, family=binomial("logit"))
  list(a=1.0, b=-coef(mod)[["(Intercept)"]])
})

# ── Pick item with highest True-vs-1PL CRF divergence ────────────────────────
divergence <- sapply(seq_len(N_ITEMS), function(j) {
  tr <- np_probs(theta_g, true_params[[j]]$a, true_params[[j]]$b, K, FLEX, true_params[[j]]$splines)
  b1 <- binom_probs(theta_g, binom_params[[j]]$a, binom_params[[j]]$b, K)
  mean(abs(tr - b1))
})
item_j <- which.max(divergence)
cat("Using item", item_j, "(CRF divergence =", round(divergence[item_j], 4), ")\n")

# ── Compute curves for chosen item ────────────────────────────────────────────
tr <- np_probs(theta_g, true_params[[item_j]]$a, true_params[[item_j]]$b, K, FLEX,
               true_params[[item_j]]$splines)
gp <- gpcm_probs(theta_g, gpcm_params[[item_j]]$a, gpcm_params[[item_j]]$b)
b1 <- binom_probs(theta_g, binom_params[[item_j]]$a, binom_params[[item_j]]$b, K)

models     <- c("True", "GPCM", "Binom 1PL")
probs_list <- list("True"=tr, "GPCM"=gp, "Binom 1PL"=b1)

# ── Build long data frame ─────────────────────────────────────────────────────
plot_df <- bind_rows(
  # ERF panel
  bind_rows(lapply(models, function(m)
    data.frame(theta=theta_g, prob=erf_fn(probs_list[[m]], K),
               model=m, panel="ERF"))),
  # One panel per category
  bind_rows(lapply(0:(K-1L), function(k)
    bind_rows(lapply(models, function(m)
      data.frame(theta=theta_g, prob=probs_list[[m]][, k+1L],
                 model=m, panel=paste0("k = ", k))))))
) |>
  mutate(
    model = factor(model, levels=c("True","GPCM","Binom 1PL")),
    panel = factor(panel, levels=c("ERF", paste0("k = ", 0:(K-1L))))
  )

# ── IMV (20% cell holdout, all items) ─────────────────────────────────────────
set.seed(77)
n_ho   <- as.integer(0.2 * N * N_ITEMS)
ho_lin <- sample(N * N_ITEMS, n_ho)
ho_mat <- matrix(FALSE, N, N_ITEMS); ho_mat[ho_lin] <- TRUE
ho_row <- row(ho_mat)[ho_lin]; ho_col <- col(ho_mat)[ho_lin]
X_ho   <- X_mat[ho_lin]
pctt   <- as.numeric(table(factor(X_ho, 0:(K-1)))) / length(X_ho)

pred_ho <- function(fn) {
  out <- matrix(NA_real_, n_ho, K)
  for (j in seq_len(N_ITEMS)) {
    cc <- which(ho_col == j); if (!length(cc)) next
    out[cc,] <- fn(j, theta_eap[ho_row[cc]])
  }
  out
}
pr_g <- pred_ho(function(j,th) gpcm_probs(th, gpcm_params[[j]]$a, gpcm_params[[j]]$b))
pr_b <- pred_ho(function(j,th) binom_probs(th, binom_params[[j]]$a, binom_params[[j]]$b, K))

df_gb <- as.data.frame(pr_b); names(df_gb) <- paste0("p1",0:(K-1))
df_gb$resp <- X_ho
for (k in 0:(K-1)) df_gb[[paste0("p2",k)]] <- pr_g[,k+1L]

imv_c_val <- {
  iis <- 0:(K-1L); om <- numeric(K)
  for (ii in iis) {
    ns <- om.tmp <- numeric()
    for (jj in iis[-match(ii,iis)]) {
      y2   <- df_gb[df_gb$resp %in% c(ii,jj),]
      resp <- ifelse(y2$resp==ii,1,0)
      p1   <- y2[[paste0("p1",ii)]]/(y2[[paste0("p1",ii)]]+y2[[paste0("p1",jj)]])
      p2   <- y2[[paste0("p2",ii)]]/(y2[[paste0("p2",ii)]]+y2[[paste0("p2",jj)]])
      om.tmp[as.character(jj)] <- imv.binary(resp,p1,p2)
      ns[as.character(jj)]     <- nrow(y2)
    }
    om[ii+1L] <- sum(om.tmp*ns)/sum(ns)
  }
  sum(om*pctt)/sum(pctt)
}
imv_t_val <- {
  om <- numeric(K-1L)
  for (ii in 0:(K-2L)) {
    resp <- ifelse(df_gb$resp<=ii,1,0)
    p1   <- rowSums(df_gb[,paste0("p1",0:ii),drop=FALSE])
    p2   <- rowSums(df_gb[,paste0("p2",0:ii),drop=FALSE])
    om[ii+1L] <- imv.binary(resp,p1,p2)
  }
  wts <- pctt[1:(K-1L)]/(1-pctt[K])
  sum(om*wts)/sum(wts)
}
cat(sprintf("IMV (GPCM vs 1PL, all items):  omega_c = %.3f  omega_t = %.3f\n",
            imv_c_val, imv_t_val))

# ── Plot ──────────────────────────────────────────────────────────────────────
model_colors <- c("True"="#000000","GPCM"="#2166ac","Binom 1PL"="#d6604d")
model_ltypes <- c("True"="solid",  "GPCM"="dashed", "Binom 1PL"="dotted")
model_widths <- c("True"=1.0,      "GPCM"=0.85,     "Binom 1PL"=0.85)

sub_label <- sprintf("GPCM vs Binom 1PL  (all %d items):  omega_c = %.3f,  omega_t = %.3f",
                     N_ITEMS, imv_c_val, imv_t_val)

p <- ggplot(plot_df,
            aes(x=theta, y=prob, colour=model, linetype=model, linewidth=model)) +
  geom_line() +
  scale_colour_manual(values=model_colors, name=NULL) +
  scale_linetype_manual(values=model_ltypes, name=NULL) +
  scale_linewidth_manual(values=model_widths, name=NULL, guide="none") +
  scale_x_continuous(breaks=c(-2, 0, 2)) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~panel, ncol=3, scales="free_y") +
  labs(x=expression(theta), y=NULL, subtitle=sub_label) +
  theme_bw(base_size=10) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill="grey92"),
        strip.text       = element_text(size=9),
        panel.spacing    = unit(0.5, "lines"),
        plot.subtitle    = element_text(size=8, hjust=0.5, colour="grey30"))

ggsave("/home/ben/Dropbox/Apps/Overleaf/binomial_model/fig_single_item_example.pdf", p, width=6.5, height=4.5)
message("Saved fig_single_item_example.pdf")
