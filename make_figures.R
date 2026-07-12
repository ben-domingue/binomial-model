## make_figures.R
## Generates all manuscript figure PDFs into the Overleaf working directory.
## Run from: /home/ben/Dropbox/Apps/Overleaf/binomial_model/
## Rscript make_figures.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(mgcv)
  library(RColorBrewer)
})

OUT  <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model"
SIM  <- "/home/ben/Dropbox/projects/binomial_model/src/sim"
VIG  <- "/home/ben/Dropbox/projects/binomial_model/src/irw"

W <- 7   # default figure width  (inches)
H <- 4.5 # default figure height (inches)

save_pdf <- function(p, name, w = W, h = H) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), plot = p,
         width = w, height = h, device = cairo_pdf)
  message("  saved ", name, ".pdf")
}

## ── shared helpers ────────────────────────────────────────────────────────────

gpcm_pr_fn <- function(theta, a, b) {
  K <- length(b) + 1L; ln <- matrix(0, nrow = length(theta), ncol = K)
  for (k in seq_len(K - 1L)) ln[, k + 1L] <- ln[, k] + a * (theta - b[k])
  ln <- ln - apply(ln, 1, max); ex <- exp(ln); ex / rowSums(ex)
}
grm_pr_fn <- function(theta, a, b) {
  K <- length(b) + 1L
  cum <- cbind(1, plogis(outer(theta, b, function(th, bk) a * (th - bk))), 0)
  pmax(cum[, 1:K] - cum[, 2:(K + 1L)], 1e-15)
}
binom_pr_fn <- function(theta, a, b, K, link) {
  eta <- a * theta - b
  p <- switch(link, logit = plogis(eta), probit = pnorm(eta),
              cloglog = 1 - exp(-exp(eta)), cauchit = 0.5 + atan(eta) / pi)
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  outer(p, 0:(K - 1L), function(pp, k) dbinom(k, K - 1L, pp))
}
es_fn <- function(pr, K) as.numeric(pr %*% 0:(K - 1L))

make_thresholds <- function(K, asym = 0) {
  n <- K - 1L; total_range <- 3.0
  if (asym == 0 || n <= 2L) return(seq(-total_range/2, total_range/2, length.out = n))
  k_idx <- 0:(n - 2L); gaps <- exp(asym * k_idx / (n - 2L))
  gaps  <- gaps / sum(gaps) * total_range
  b <- cumsum(c(-total_range / 2, gaps[-length(gaps)]))
  c(b, b[n - 1L] + gaps[n - 1L])
}

model_labels <- c(
  gpcm              = "GPCM",       pcm               = "PCM",
  grm               = "GRM",        tutz              = "Tutz",
  binom_2pl_logit   = "Binom 2PL (logit)",
  binom_1pl_logit   = "Binom 1PL (logit)",
  binom_1pl_scobit  = "Binom 1PL (scobit)",
  binom_1pl_cloglog = "Binom 1PL (cloglog)",
  binom_1pl_cauchit = "Binom 1PL (cauchit)",
  np_spline         = "NP spline"
)
model_colors <- c(
  "GPCM"               = "#1b7837", "PCM"                = "#762a83",
  "GRM"                = "#2166ac", "Tutz"               = "#543005",
  "Binom 2PL (logit)"  = "#d6604d", "Binom 1PL (logit)"  = "#f4a582",
  "Binom 1PL (scobit)" = "#f1a340", "Binom 1PL (cloglog)"= "#c2a5cf",
  "Binom 1PL (cauchit)"= "#a6dba0", "NP spline"          = "#e7298a"
)

## normalise old empirical model names (pre-June-2026 compute runs)
normalise_model <- function(m) {
  m <- as.character(m)
  m[m == "1PL+logit"]   <- "Binom 1PL (logit)"
  m[m == "2PL+logit"]   <- "Binom 2PL (logit)"
  m[m == "1PL+cloglog"] <- "Binom 1PL (cloglog)"
  m
}

theme_set(
  theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92"),
          legend.position  = "none",
          plot.margin      = margin(5, 5, 5, 5))
)

## ═══════════════════════════════════════════════════════════════════════════
## SECTION 1: GPCM DGP SIMULATION FIGURES
## ═══════════════════════════════════════════════════════════════════════════
message("\n── Loading simulation data ──")
res <- readRDS(file.path(SIM, "gpcm_binomial_results.rds")) |>
  filter(K == 5) |>
  mutate(model_label = factor(model_labels[as.character(model)], levels = model_labels))

lpos_fn <- function(data, y_var) {
  grps <- split(data, list(as.character(data$model_label)), drop = TRUE)
  do.call(rbind, lapply(grps, function(g) {
    x_anch <- mean(g$asym[g$asym >= quantile(g$asym, 0.9)])
    fit    <- tryCatch(mgcv::gam(as.formula(paste(y_var, "~ s(asym, bs='cs')")), data = g),
                       error = function(e) NULL)
    y_anch <- if (!is.null(fit))
      as.numeric(predict(fit, newdata = data.frame(asym = x_anch)))
    else mean(g[[y_var]][g$asym >= quantile(g$asym, 0.9)], na.rm = TRUE)
    data.frame(model_label = g$model_label[1], asym = x_anch, y_val = y_anch,
               stringsAsFactors = FALSE)
  }))
}

make_imv_pair <- function(data, v1, v2, l1, l2, y_lab, excl = "GPCM") {
  long <- bind_rows(
    mutate(data, y = .data[[v1]], metric = l1),
    mutate(data, y = .data[[v2]], metric = l2)
  ) |> mutate(metric = factor(metric, levels = c(l1, l2)))
  lpos <- do.call(rbind, lapply(
    split(long, list(as.character(long$model_label), as.character(long$metric)), drop = TRUE),
    function(g) {
      x_anch <- mean(g$asym[g$asym >= quantile(g$asym, 0.9)])
      fit <- tryCatch(mgcv::gam(y ~ s(asym, bs = "cs"), data = g), error = function(e) NULL)
      y_anch <- if (!is.null(fit))
        as.numeric(predict(fit, newdata = data.frame(asym = x_anch)))
      else mean(g$y[g$asym >= quantile(g$asym, 0.9)], na.rm = TRUE)
      data.frame(model_label = g$model_label[1], metric = g$metric[1],
                 asym = x_anch, y_val = y_anch, stringsAsFactors = FALSE)
    }
  ))
  lpos$metric <- factor(lpos$metric, levels = c(l1, l2))
  ggplot(long, aes(x = asym, y = y, colour = model_label, group = model_label)) +
    annotate("segment", x = -log(8), xend = log(8), y = 0, yend = 0,
             linetype = "dashed", colour = "grey60") +
    geom_smooth(data = function(d) if (is.null(excl)) d else filter(d, !model_label %in% excl),
                method = "gam", formula = y ~ s(x, bs = "cs"),
                se = TRUE, linewidth = 0.8, alpha = 0.12) +
    annotate("rect", xmin = log(8), xmax = log(8) + 3, ymin = -Inf, ymax = Inf,
             fill = "white", colour = NA) +
    geom_label_repel(
      data = lpos, aes(x = asym, y = y_val, label = model_label),
      size = 2.6, hjust = 0, nudge_x = log(8) + 0.3 - lpos$asym,
      direction = "y", segment.size = 0.3, segment.color = "grey50",
      box.padding = 0.2, label.padding = 0.15,
      show.legend = FALSE, xlim = c(log(8), log(8) + 3)
    ) +
    facet_wrap(~ metric) +
    scale_colour_manual(values = model_colors) +
    scale_x_continuous(
      name = "Gap asymmetry s  (s < 0 = wider lower;  s = 0 = symmetric;  s > 0 = wider upper)",
      limits = c(-log(8), log(8) + 3),
      breaks = round(c(-log(8), -1, 0, 1, log(8)), 2),
      labels = c("-ln8", "-1", "0", "1", "ln8"),
      expand = expansion(mult = c(0.02, 0))
    ) +
    coord_cartesian(clip = "off") + ylab(y_lab)
}

## ── fig_erf_overview ──────────────────────────────────────────────────────
message("fig_erf_overview")
theta_g <- seq(-3, 3, by = 0.05); K_ex <- 5L; a_ex <- 1.5
b_ov    <- make_thresholds(K_ex, log(8))
gpcm_es <- es_fn(gpcm_pr_fn(theta_g, a_ex, b_ov), K_ex)
p2      <- optim(c(a_ex, 0),
                 function(p) mean((gpcm_es -
                   es_fn(binom_pr_fn(theta_g, p[1], p[2], K_ex, "logit"), K_ex))^2),
                 method = "BFGS")$par
binom_es <- es_fn(binom_pr_fn(theta_g, p2[1], p2[2], K_ex, "logit"), K_ex)
erf_df <- rbind(
  data.frame(theta = theta_g, es = gpcm_es,  model = "GPCM"),
  data.frame(theta = theta_g, es = binom_es, model = "Binom 2PL (logit)")
)
erf_df$model <- factor(erf_df$model, levels = c("GPCM", "Binom 2PL (logit)"))
p <- ggplot(erf_df, aes(x = theta, y = es, colour = model, linetype = model)) +
  geom_line(linewidth = 1.0) +
  scale_colour_manual(values = c("GPCM" = "black", "Binom 2PL (logit)" = "#d6604d"), name = NULL) +
  scale_linetype_manual(values = c("GPCM" = "solid", "Binom 2PL (logit)" = "dashed"), name = NULL) +
  scale_x_continuous(name = expression(theta), breaks = -2:2) +
  scale_y_continuous(name = expression(E*"["*X*"|"*theta*"]"), limits = c(0, K_ex - 1L)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_pdf(p, "fig_erf_overview", w = 4.5, h = 3.5)

## ── fig_erf_grid (3 conditions × 8 models) ───────────────────────────────
message("fig_erf_grid")
panel_levels <- c("PCM","GRM","NP spline","Binom 2PL (logit)","Binom 1PL (logit)",
                  "Binom 1PL (scobit)","Binom 1PL (cloglog)","Binom 1PL (cauchit)")
panel_color  <- "#e6550d"
panel_colors <- setNames(c(rep(panel_color,2),"#e7298a",rep(panel_color,5)), panel_levels)

make_erf_panel <- function(s_val, label) {
  b_ex    <- make_thresholds(K_ex, s_val)
  gpcm_es <- es_fn(gpcm_pr_fn(theta_g, a_ex, b_ex), K_ex)
  mse     <- function(es) mean((gpcm_es - es)^2)
  pcm_b   <- optim(b_ex, function(b) mse(es_fn(gpcm_pr_fn(theta_g,1.0,sort(b)),K_ex)), method="BFGS")$par
  grm_p   <- optim(c(a_ex,b_ex), function(p) { if(p[1]<=0) return(1e9)
    mse(es_fn(grm_pr_fn(theta_g,p[1],sort(p[-1])),K_ex)) }, method="BFGS")$par
  p2l     <- optim(c(a_ex,0), function(p) mse(es_fn(binom_pr_fn(theta_g,p[1],p[2],K_ex,"logit"),K_ex)), method="BFGS")$par
  opt1    <- function(link) optim(0, function(b) mse(es_fn(binom_pr_fn(theta_g,1,b,K_ex,link),K_ex)),
                                   method="Brent",lower=-3,upper=3)$par
  sco_es  <- function(b,la) es_fn(outer(pmin(pmax(plogis(theta_g-b)^exp(la),1e-10),1-1e-10),
                                         0:(K_ex-1L),function(pp,k) dbinom(k,K_ex-1L,pp)),K_ex)
  sc_p    <- optim(c(0,0), function(p) mse(sco_es(p[1],p[2])), method="BFGS")$par
  df_np   <- data.frame(y=gpcm_es/(K_ex-1L), th=theta_g)
  fit_np  <- scam::scam(y~s(th,bs="mpi"), data=df_np)
  p_np    <- pmin(pmax(as.numeric(predict(fit_np)),0),1)
  es_np   <- es_fn(outer(p_np,0:(K_ex-1L),function(p,k) dbinom(k,K_ex-1L,p)),K_ex)
  m_df <- rbind(
    data.frame(panel="PCM",                 theta=theta_g,es=es_fn(gpcm_pr_fn(theta_g,1.0,sort(pcm_b)),K_ex)),
    data.frame(panel="GRM",                 theta=theta_g,es=es_fn(grm_pr_fn(theta_g,grm_p[1],sort(grm_p[-1])),K_ex)),
    data.frame(panel="NP spline",           theta=theta_g,es=es_np),
    data.frame(panel="Binom 2PL (logit)",   theta=theta_g,es=es_fn(binom_pr_fn(theta_g,p2l[1],p2l[2],K_ex,"logit"),K_ex)),
    data.frame(panel="Binom 1PL (logit)",   theta=theta_g,es=es_fn(binom_pr_fn(theta_g,1,opt1("logit"),K_ex,"logit"),K_ex)),
    data.frame(panel="Binom 1PL (scobit)",  theta=theta_g,es=sco_es(sc_p[1],sc_p[2])),
    data.frame(panel="Binom 1PL (cloglog)", theta=theta_g,es=es_fn(binom_pr_fn(theta_g,1,opt1("cloglog"),K_ex,"cloglog"),K_ex)),
    data.frame(panel="Binom 1PL (cauchit)", theta=theta_g,es=es_fn(binom_pr_fn(theta_g,1,opt1("cauchit"),K_ex,"cauchit"),K_ex))
  )
  m_df$panel <- factor(m_df$panel, levels=panel_levels)
  ref_df <- do.call(rbind, lapply(panel_levels, function(m)
    data.frame(panel=factor(m,levels=panel_levels), theta=theta_g, es=gpcm_es)))
  ggplot() +
    geom_line(data=ref_df, aes(x=theta,y=es), colour="black", linewidth=1.2, linetype="dashed") +
    geom_line(data=m_df,   aes(x=theta,y=es,colour=panel), linewidth=0.9) +
    facet_wrap(~panel, ncol=4) +
    scale_colour_manual(values=panel_colors, guide="none") +
    scale_x_continuous(name=expression(theta), breaks=-2:2) +
    scale_y_continuous(name="E[X | θ]", limits=c(0,K_ex-1)) +
    ggtitle(label) +
    theme_bw(base_size=10) + theme(strip.background=element_rect(fill="grey92"),
                                    panel.grid.minor=element_blank())
}
suppressPackageStartupMessages(library(scam))
g0 <- make_erf_panel(0,       "Symmetric (s = 0)")
gu <- make_erf_panel(log(8),  "Upper-wider (s = ln 8)")
gl <- make_erf_panel(-log(8), "Lower-wider (s = -ln 8)")
p  <- ggpubr::ggarrange(g0, gu, gl, ncol=1, nrow=3)
save_pdf(p, "fig_erf_grid", w=9, h=13)

## ── fig_imv_gpcm ──────────────────────────────────────────────────────────
message("fig_imv_gpcm")
p <- make_imv_pair(res, "imv_c", "imv_t", "ω_c (category)", "ω_t (threshold)",
                   "IMV(model, GPCM)  [positive = GPCM wins]")
save_pdf(p, "fig_imv_gpcm", w=8, h=4.5)

## ── fig_imv_floor ─────────────────────────────────────────────────────────
message("fig_imv_floor")
res_1pl <- filter(res, model_label != "Binom 1PL (logit)")
p <- make_imv_pair(res_1pl, "imv_c_1pl", "imv_t_1pl", "ω_c (category)", "ω_t (threshold)",
                   "IMV(1PL-logit, model)  [positive = model wins]", excl=NULL)
save_pdf(p, "fig_imv_floor", w=8, h=4.5)

## ── fig_rmse_sim ──────────────────────────────────────────────────────────
message("fig_rmse_sim")
lpos <- lpos_fn(res, "rmse")
p <- ggplot(res, aes(x=asym, y=rmse, colour=model_label, group=model_label)) +
  geom_smooth(data=function(d) filter(d, model_label!="GPCM"),
              method="gam", formula=y~s(x,bs="cs"), se=TRUE, linewidth=0.8, alpha=0.12) +
  annotate("rect", xmin=log(8), xmax=log(8)+3, ymin=-Inf, ymax=Inf, fill="white", colour=NA) +
  geom_label_repel(data=lpos, aes(x=asym, y=y_val, label=model_label),
                   size=2.6, hjust=0, nudge_x=log(8)+0.3-lpos$asym,
                   direction="y", segment.size=0.3, segment.color="grey50",
                   box.padding=0.2, label.padding=0.15, show.legend=FALSE, xlim=c(log(8),log(8)+3)) +
  scale_colour_manual(values=model_colors) +
  scale_x_continuous(name="Gap asymmetry s",
                     limits=c(-log(8), log(8)+3),
                     breaks=round(c(-log(8),-1,0,1,log(8)),2),
                     labels=c("-ln8","-1","0","1","ln8"),
                     expand=expansion(mult=c(0.02,0))) +
  coord_cartesian(clip="off") + ylab("RMSE (expected score, 0–1 scale)")
save_pdf(p, "fig_rmse_sim", w=7, h=4.5)

## ── fig_kl_category ── generated by src/npspline.R, not here ─────────────

## ═══════════════════════════════════════════════════════════════════════════
## SECTION 2: NP DGP FIGURES
## ═══════════════════════════════════════════════════════════════════════════
message("\n── Loading NP simulation data ──")
res_np <- readRDS(file.path(SIM, "gpcm_binomial_results_np.rds")) |>
  filter(K == 5) |>
  mutate(model_label = factor(model_labels[as.character(model)], levels = model_labels))

make_imv_pair_np <- function(data, v1, v2, l1, l2, y_lab, excl="GPCM") {
  long <- bind_rows(
    mutate(data, y=.data[[v1]], metric=l1),
    mutate(data, y=.data[[v2]], metric=l2)
  ) |> mutate(metric=factor(metric,levels=c(l1,l2)))
  lpos <- do.call(rbind, lapply(
    split(long, list(as.character(long$model_label), as.character(long$metric)), drop=TRUE),
    function(g) {
      x_anch <- mean(g$flex[g$flex >= quantile(g$flex, 0.9)])
      fit <- tryCatch(mgcv::gam(y~s(flex,bs="cs"),data=g), error=function(e) NULL)
      y_anch <- if(!is.null(fit)) as.numeric(predict(fit,newdata=data.frame(flex=x_anch)))
                else mean(g$y[g$flex >= quantile(g$flex,0.9)], na.rm=TRUE)
      data.frame(model_label=g$model_label[1], metric=g$metric[1],
                 flex=x_anch, y_val=y_anch, stringsAsFactors=FALSE)
    }
  ))
  lpos$metric <- factor(lpos$metric, levels=c(l1,l2))
  ggplot(long, aes(x=flex, y=y, colour=model_label, group=model_label)) +
    annotate("segment", x=0, xend=1, y=0, yend=0, linetype="dashed", colour="grey60") +
    geom_smooth(data=function(d) if(is.null(excl)) d else filter(d,!model_label %in% excl),
                method="gam", formula=y~s(x,bs="cs"), se=TRUE, linewidth=0.8, alpha=0.12) +
    annotate("rect", xmin=1, xmax=1.35, ymin=-Inf, ymax=Inf, fill="white", colour=NA) +
    geom_label_repel(data=lpos, aes(x=flex,y=y_val,label=model_label),
                     size=2.6, hjust=0, nudge_x=1.05-lpos$flex,
                     direction="y", segment.size=0.3, segment.color="grey50",
                     box.padding=0.2, label.padding=0.15,
                     show.legend=FALSE, xlim=c(1,1.35)) +
    facet_wrap(~metric) +
    scale_colour_manual(values=model_colors) +
    scale_x_continuous(name="flex  (0 = GPCM boundaries;  1 = fully nonparametric)",
                       limits=c(0,1.35), breaks=c(0,0.25,0.5,0.75,1),
                       labels=c("0\n(GPCM)","0.25","0.5","0.75","1\n(NP)"),
                       expand=expansion(mult=c(0.02,0))) +
    coord_cartesian(clip="off") + ylab(y_lab)
}

## ── fig_np_imv ────────────────────────────────────────────────────────────
message("fig_np_imv")
p <- make_imv_pair_np(res_np, "imv_c", "imv_t", "ω_c (category)", "ω_t (threshold)",
                      "IMV(model, GPCM)  [positive = GPCM wins;  negative = model wins]")
save_pdf(p, "fig_np_imv", w=8, h=4.5)

## ── fig_np_imv_oracle ─────────────────────────────────────────────────────
## IMV vs. the oracle (true NP DGP, true theta) — SI Figure (sec:si-oracle)
message("fig_np_imv_oracle")
res_np_oracle <- readRDS(file.path(SIM, "gpcm_binomial_results_np_oracle.rds")) |>
  filter(K == 5) |>
  mutate(model_label = factor(model_labels[as.character(model)], levels = model_labels))
p <- make_imv_pair_np(res_np_oracle, "imv_c_oracle", "imv_t_oracle",
                      "ω_c (category)", "ω_t (threshold)",
                      "IMV(oracle, model)  [positive = oracle wins]", excl = NULL)
save_pdf(p, "fig_np_imv_oracle", w=8, h=4.5)

## ── fig_np_boundary_mixture ── generated by src/fig_np_boundary_mixture.R, not here ──

## ── fig_np_extreme_categories ── generated by src/fig_np_extreme_categories.R, not here ──

## ── fig_implied_bounds ────────────────────────────────────────────────────
message("fig_implied_bounds")
K_np <- 5L
## Designed DGP: GRM backbone with asymmetric spacing — k=1,2 clustered on the left,
## k=3,4 clustered on the right, large gap in the middle. The ERF is symmetric around
## theta=0 and spans the full range (fits a logistic well), but the boundary structure
## is completely different from the evenly-spaced binomial fan.
knots_exotic <- seq(-4, 4, length.out = 13L)
b_back <- c(-1.5, -1.0, 1.0, 1.5); a_back <- 1.0
set.seed(17)
add_spline_char <- function(v, noise_sd = 0.20) {
  n <- length(v); d <- diff(v)
  d_noisy <- d * exp(rnorm(n - 1L, 0, noise_sd))
  v_new <- c(v[1], v[1] + cumsum(d_noisy))
  pmin(pmax((v_new - v_new[1]) / (v_new[n] - v_new[1]) * (v[n] - v[1]) + v[1], 0), 1)
}
exotic_vals <- lapply(b_back, function(b)
  add_spline_char(plogis(a_back * (knots_exotic - b))))
spl_exotic <- lapply(exotic_vals, function(v) splinefun(knots_exotic, v, method = "monoH.FC"))
np_cum_f1 <- do.call(cbind, lapply(spl_exotic, function(s)
  pmin(pmax(s(pmin(pmax(theta_g, -4), 4)), 0), 1)))
np_cum_f1 <- t(apply(np_cum_f1, 1, sort, decreasing = TRUE))
## Fit GRM and binomial to minimise boundary RMSE.
grm_pf1 <- optim(c(1.0, -2.0, -1.5, 1.5, 2.0), function(p) {
  if (p[1] <= 0) return(1e9)
  grm_cum <- outer(theta_g, sort(p[-1]), function(th, bk) plogis(p[1] * (th - bk)))
  mean((np_cum_f1 - grm_cum)^2)
}, method = "BFGS")$par
grm_cum_f1 <- outer(theta_g, sort(grm_pf1[-1]),
                    function(th, bk) plogis(grm_pf1[1] * (th - bk)))
b1_f1 <- optim(0, function(b) {
  p_b <- plogis(theta_g - b)
  binom_cum <- sapply(seq_len(K_np - 1L), function(k)
    pbinom(k - 1L, K_np - 1L, p_b, lower.tail = FALSE))
  mean((np_cum_f1 - binom_cum)^2)
}, method = "Brent", lower = -3, upper = 3)$par
p_f1 <- plogis(theta_g - b1_f1)
binom_cum_f1 <- sapply(seq_len(K_np - 1L), function(k)
  pbinom(k - 1L, K_np - 1L, p_f1, lower.tail = FALSE))
library(patchwork)
bound_pal4 <- brewer.pal(4,"RdYlBu")
erf_cols <- c("True NP DGP"="black", "GRM fit"="#2166ac", "Binom 1PL fit"="#d6604d")
## ERF panel
erf_df <- data.frame(
  theta = rep(theta_g, 3),
  erf   = c(rowSums(np_cum_f1), rowSums(grm_cum_f1), rowSums(binom_cum_f1)),
  model = rep(c("True NP DGP","GRM fit","Binom 1PL fit"), each=length(theta_g))
)
erf_df$model <- factor(erf_df$model, levels=names(erf_cols))
erf_df$panel <- "Expected response functions"
p_erf <- ggplot(erf_df, aes(x=theta, y=erf, colour=model, linetype=model)) +
  geom_line(linewidth=0.9) +
  facet_wrap(~panel) +
  scale_colour_manual(values=erf_cols, name=NULL) +
  scale_linetype_manual(values=c("True NP DGP"="solid","GRM fit"="dashed","Binom 1PL fit"="dotdash"), name=NULL) +
  scale_x_continuous(breaks=-2:2, name=expression(theta)) +
  scale_y_continuous(name=expression(E(X~"|"~theta))) +
  theme_bw(base_size=10) + theme(strip.background=element_rect(fill="grey92"),
                                   panel.grid.minor=element_blank(),
                                   legend.position=c(0.02, 0.98),
                                   legend.justification=c(0, 1),
                                   legend.background=element_rect(fill="white", colour=NA),
                                   legend.key.size=unit(0.9,"lines"),
                                   legend.text=element_text(size=7.5))
## Boundary panels
panel_labels <- c(np="True NP DGP\n(K-1 free splines)",
                   grm="GRM fit\n(K-1 free logistic curves)",
                   binom="Binom 1PL fit\n(1 logistic → K-1 constrained curves)")
bound_df2 <- do.call(rbind, lapply(names(panel_labels), function(src) {
  mat <- switch(src, np=np_cum_f1, grm=grm_cum_f1, binom=binom_cum_f1)
  do.call(rbind, lapply(seq_len(K_np-1L), function(k)
    data.frame(theta=theta_g, prob=mat[,k], boundary=paste0("k = ",k), src=src)))
}))
bound_df2$src <- factor(bound_df2$src, levels=names(panel_labels), labels=panel_labels)
p_bounds <- ggplot(bound_df2, aes(x=theta, y=prob, colour=boundary)) +
  geom_line(linewidth=0.9) +
  facet_wrap(~src, ncol=3) +
  scale_colour_manual(values=bound_pal4, name="Boundary") +
  scale_x_continuous(breaks=-2:2, name=expression(theta)) +
  scale_y_continuous(name=expression(P(X>=k~"|"~theta)), limits=c(0,1)) +
  theme_bw(base_size=10) + theme(strip.background=element_rect(fill="grey92"),
                                   panel.grid.minor=element_blank(),
                                   legend.position="right")
p <- p_erf + p_bounds + plot_layout(widths=c(1,3))
save_pdf(p, "fig_implied_bounds", w=12, h=3.5)

## ── fig_cross_section ─────────────────────────────────────────────────────
message("fig_cross_section")
s1c <- res    |> filter(abs(asym)>1) |> group_by(model_label) |>
  summarise(imv=mean(imv_c_1pl,na.rm=TRUE),.groups="drop") |>
  mutate(section="Section 1\n(GPCM DGP, |s|>1)", metric="ω_c (category)")
s1t <- res    |> filter(abs(asym)>1) |> group_by(model_label) |>
  summarise(imv=mean(imv_t_1pl,na.rm=TRUE),.groups="drop") |>
  mutate(section="Section 1\n(GPCM DGP, |s|>1)", metric="ω_t (threshold)")
s2c <- res_np |> filter(flex>0.67) |> group_by(model_label) |>
  summarise(imv=mean(imv_c_1pl,na.rm=TRUE),.groups="drop") |>
  mutate(section="Section 2\n(NP DGP, flex>0.67)", metric="ω_c (category)")
s2t <- res_np |> filter(flex>0.67) |> group_by(model_label) |>
  summarise(imv=mean(imv_t_1pl,na.rm=TRUE),.groups="drop") |>
  mutate(section="Section 2\n(NP DGP, flex>0.67)", metric="ω_t (threshold)")
comp_df <- rbind(s1c,s1t,s2c,s2t) |>
  filter(model_label != "Binom 1PL (logit)") |>
  mutate(section=factor(section, levels=c("Section 1\n(GPCM DGP, |s|>1)",
                                           "Section 2\n(NP DGP, flex>0.67)")),
         metric=factor(metric, levels=c("ω_c (category)","ω_t (threshold)")))
ord <- comp_df |> filter(grepl("Section 1",as.character(section)), metric=="ω_c (category)") |>
  arrange(imv) |> pull(model_label) |> as.character()
comp_df$model_label <- factor(comp_df$model_label, levels=ord)
p <- ggplot(comp_df, aes(x=model_label, y=imv, fill=section)) +
  geom_col(position=position_dodge(width=0.7), width=0.65) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
  facet_wrap(~metric, ncol=2) +
  scale_fill_manual(values=c("Section 1\n(GPCM DGP, |s|>1)"="#2166ac",
                              "Section 2\n(NP DGP, flex>0.67)"="#d6604d"), name=NULL) +
  labs(x=NULL, y="Mean IMV(1PL-logit, model)  [positive = model wins]") +
  coord_flip() +
  theme_bw(base_size=10) + theme(strip.background=element_rect(fill="grey92"),
                                   legend.position="bottom", panel.grid.minor=element_blank())
save_pdf(p, "fig_cross_section", w=8, h=4.5)

## ═══════════════════════════════════════════════════════════════════════════
## EMPIRICAL FIGURES
## ═══════════════════════════════════════════════════════════════════════════
message("\n── Loading empirical data ──")
poly_dir <- file.path(VIG, "poly_compare_data")
fs       <- list.files(poly_dir, pattern="^multi_.+\\.rds$", full.names=TRUE)
fs       <- fs[!grepl("multi_(results|summary)\\.rds$", fs)]
raw      <- lapply(fs, readRDS)
nms      <- sub("^multi_","",sub("\\.rds$","",basename(fs)))
names(raw) <- nms
raw      <- Filter(Negate(is.null), raw)
raw      <- lapply(raw, function(r) { if(is.null(r[["K"]])) r$K <- r$K_range+1L; r })
results  <- Filter(function(r) r$K >= 3, raw)
all_rmse <- do.call(rbind, lapply(results, `[[`, "rmse")) |>
  mutate(model = normalise_model(model))
all_imv  <- do.call(rbind, Filter(Negate(is.null), lapply(results, `[[`, "imv"))) |>
  mutate(model = normalise_model(model))
all_pair <- do.call(rbind, Filter(Negate(is.null), lapply(results, `[[`, "pair_imv"))) |>
  mutate(comparison = sub("1PL\\+logit", "Binom 1PL (logit)",
                     sub("2PL\\+logit", "Binom 2PL (logit)", comparison)))
meta     <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(dataset=nm, K=r$K, ni=r$ni, N=r$N, stringsAsFactors=FALSE)
}))

model_levels <- c("PCM","GPCM","GRM","Tutz",
                  "Binom 1PL (logit)","Binom 2PL (logit)","Binom 1PL (cloglog)","NP spline")
model_type   <- c(PCM="Standard",GPCM="Standard",GRM="Standard",Tutz="Standard",
                  `Binom 1PL (logit)`="Binomial",`Binom 2PL (logit)`="Binomial",
                  `Binom 1PL (cloglog)`="Binomial",`NP spline`="Binomial")
pal <- c(Standard="#2166ac", Binomial="#d6604d")

all_rmse <- all_rmse |>
  mutate(model=factor(model,levels=model_levels), mod_type=model_type[as.character(model)])
ds_med   <- all_rmse |> group_by(dataset) |> summarise(med_rmse=median(rmse),.groups="drop") |> arrange(med_rmse)
all_rmse <- all_rmse |> mutate(dataset=factor(dataset, levels=ds_med$dataset))
ranks    <- all_rmse |> group_by(dataset) |> mutate(rank=rank(rmse,ties.method="min")) |> ungroup()
all_imv  <- all_imv  |> mutate(model=factor(model,levels=model_levels),
                                mod_type=model_type[as.character(model)])
mean_imv <- all_imv  |> group_by(model,mod_type) |>
  summarise(mean_imv_c=mean(imv_c,na.rm=TRUE), se_imv_c=sd(imv_c,na.rm=TRUE)/sqrt(sum(!is.na(imv_c))),
            mean_imv_t=mean(imv_t,na.rm=TRUE), se_imv_t=sd(imv_t,na.rm=TRUE)/sqrt(sum(!is.na(imv_t))),
            .groups="drop")

## ── fig_rmse_empirical ────────────────────────────────────────────────────
message("fig_rmse_empirical")
model_ord <- all_rmse |> group_by(model) |> summarise(med=median(rmse),.groups="drop") |>
  arrange(med) |> pull(model)
all_rmse  <- all_rmse |> mutate(model=factor(model, levels=model_ord))
p <- ggplot(all_rmse, aes(x=model, y=rmse)) +
  geom_violin(aes(fill=mod_type), width=0.9, alpha=0.30, colour=NA, scale="width") +
  geom_boxplot(aes(colour=mod_type), width=0.15, fill="white", outlier.shape=1,
               outlier.size=1, outlier.alpha=0.5, linewidth=0.4) +
  scale_fill_manual(values=pal, name="Family") +
  scale_colour_manual(values=pal, guide="none") +
  coord_flip() +
  labs(x=NULL, y="CV RMSE (rescaled response)") +
  theme_minimal(base_size=12) +
  theme(legend.position="right", panel.grid.major.y=element_blank())
save_pdf(p, "fig_rmse_empirical", w=7, h=4.5)

## ── fig_ranks ─────────────────────────────────────────────────────────────
message("fig_ranks")
rank_se <- ranks |> group_by(model,mod_type) |>
  summarise(mean_rank=mean(rank), se=sd(rank)/sqrt(n()), .groups="drop")
p <- ggplot(rank_se, aes(x=reorder(model,mean_rank), y=mean_rank, fill=mod_type,
                          ymin=mean_rank-se, ymax=mean_rank+se)) +
  geom_col(width=0.6, alpha=0.85) +
  geom_errorbar(width=0.25, colour="grey40") +
  scale_fill_manual(values=pal, name="Family") +
  coord_flip() + labs(x=NULL, y="Mean rank (lower = better)") +
  theme_minimal(base_size=13) + theme(legend.position="right")
save_pdf(p, "fig_ranks", w=5.5, h=3.5)

## ── fig_imv_empirical ─────────────────────────────────────────────────────
message("fig_imv_empirical")
imv_long <- mean_imv |>
  pivot_longer(c(mean_imv_c,mean_imv_t), names_to="metric", values_to="imv") |>
  mutate(se=ifelse(metric=="mean_imv_c",se_imv_c,se_imv_t),
         metric=recode(metric,mean_imv_c="ω_c (category)",mean_imv_t="ω_t (thresholded)"))
p <- ggplot(imv_long, aes(x=reorder(model,-imv), y=imv, fill=mod_type, ymin=imv-se, ymax=imv+se)) +
  geom_col(width=0.6, alpha=0.85) +
  geom_errorbar(width=0.25, colour="grey40") +
  scale_fill_manual(values=pal, name="Family") +
  facet_wrap(~metric) + coord_flip() +
  labs(x=NULL, y="Mean IMV (higher = better)") +
  theme_minimal(base_size=13) + theme(legend.position="right")
save_pdf(p, "fig_imv_empirical", w=8, h=4)

## ── fig_pair_imv ──────────────────────────────────────────────────────────
message("fig_pair_imv")
pair_long <- all_pair |>
  filter(comparison %in% c("Binom 1PL (logit)→PCM","Binom 2PL (logit)→GPCM")) |>
  pivot_longer(c(imv_c,imv_t), names_to="metric", values_to="imv") |>
  mutate(metric=recode(metric,imv_c="ω_c (category)",imv_t="ω_t (thresholded)"),
         comparison=factor(comparison, levels=c("Binom 1PL (logit)→PCM","Binom 2PL (logit)→GPCM")))
panel_stats <- pair_long |> group_by(metric,comparison) |>
  summarise(mean_imv=mean(imv,na.rm=TRUE), .groups="drop") |>
  mutate(mean_x=pmin(pmax(mean_imv,-0.1),0.1))
pair_long <- pair_long |> group_by(metric,comparison) |>
  mutate(imv=pmin(pmax(imv,-0.1),0.1)) |> ungroup()
p <- ggplot(pair_long, aes(x=imv)) +
  geom_histogram(bins=20, fill="#2166ac", colour="white", alpha=0.85) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey40") +
  geom_vline(data=panel_stats, aes(xintercept=mean_x), colour="firebrick", linewidth=0.8) +
  geom_text(data=panel_stats, aes(label=paste0("mean = ",round(mean_imv,3))),
            x=Inf, y=Inf, hjust=1.1, vjust=1.5, size=3.3, colour="firebrick", inherit.aes=FALSE) +
  facet_grid(metric~comparison) +
  labs(x="Pairwise IMV (positive = standard polytomous model better)", y="Count") +
  theme_minimal(base_size=12)
save_pdf(p, "fig_pair_imv", w=8, h=5)

## ── fig_scatter_pairs ─────────────────────────────────────────────────────
message("fig_scatter_pairs")
make_scatter <- function(mod_x, mod_y) {
  d <- all_rmse |> filter(model %in% c(mod_x,mod_y)) |>
    select(dataset,model,rmse) |> pivot_wider(names_from=model,values_from=rmse) |>
    left_join(meta, by="dataset")
  ggplot(d, aes(x=.data[[mod_x]], y=.data[[mod_y]], size=ni)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey60") +
    geom_point(colour="#2166ac", alpha=0.75) +
    scale_size_continuous(range=c(2,7), name="# items") +
    labs(x=paste("CV RMSE —",mod_x), y=paste("CV RMSE —",mod_y)) +
    theme_minimal(base_size=11)
}
p1 <- make_scatter("Binom 2PL (logit)","GPCM")
p2 <- make_scatter("Binom 1PL (logit)","PCM")
p  <- ggpubr::ggarrange(p1, p2, ncol=2, common.legend=TRUE, legend="right")
save_pdf(p, "fig_scatter_pairs", w=9, h=4.5)

## ── fig_gap_by_k ──────────────────────────────────────────────────────────
message("fig_gap_by_k")
gap_k <- all_rmse |> filter(model %in% c("GPCM","Binom 2PL (logit)")) |>
  select(dataset,model,rmse) |> pivot_wider(names_from=model,values_from=rmse) |>
  left_join(meta |> select(dataset,K,ni,N), by="dataset") |>
  mutate(gap=`Binom 2PL (logit)`-GPCM)
p <- ggplot(gap_k, aes(x=K, y=gap)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
  geom_jitter(width=0.08, size=3, colour="#2166ac", alpha=0.75) +
  stat_summary(aes(x=K,y=gap), fun=mean, geom="point", shape=18, size=5,
               colour="black", inherit.aes=FALSE) +
  scale_x_continuous(breaks=3:8) +
  labs(x="K (number of response categories)",
       y="RMSE(Binom 2PL logit) - RMSE(GPCM)") +
  theme_minimal(base_size=13)
save_pdf(p, "fig_gap_by_k", w=6, h=4)

## ── fig_imv_1pl ───────────────────────────────────────────────────────────
## IMV vs. Binom 1PL (logit) floor — only available in new-format files
message("fig_imv_1pl")
all_imv_1pl <- do.call(rbind, Filter(Negate(is.null), lapply(results, `[[`, "imv_1pl")))
if (!is.null(all_imv_1pl) && nrow(all_imv_1pl) > 0) {
  all_imv_1pl <- all_imv_1pl |>
    mutate(model = normalise_model(model),
           model = factor(model, levels = model_levels),
           mod_type = model_type[as.character(model)])
  mean_imv_1pl <- all_imv_1pl |>
    group_by(model, mod_type) |>
    summarise(mean_imv_c = mean(imv_c, na.rm = TRUE),
              se_imv_c   = sd(imv_c, na.rm = TRUE) / sqrt(sum(!is.na(imv_c))),
              mean_imv_t = mean(imv_t, na.rm = TRUE),
              se_imv_t   = sd(imv_t, na.rm = TRUE) / sqrt(sum(!is.na(imv_t))),
              .groups = "drop")
  imv_1pl_long <- mean_imv_1pl |>
    pivot_longer(c(mean_imv_c, mean_imv_t), names_to = "metric", values_to = "imv") |>
    mutate(se     = ifelse(metric == "mean_imv_c", se_imv_c, se_imv_t),
           metric = recode(metric,
                           mean_imv_c = "ω_c (category)",
                           mean_imv_t = "ω_t (thresholded)"))
  p <- ggplot(imv_1pl_long,
              aes(x = reorder(model, -imv), y = imv, fill = mod_type,
                  ymin = imv - se, ymax = imv + se)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(width = 0.25, colour = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = pal, name = "Family") +
    facet_wrap(~metric) + coord_flip() +
    labs(x = NULL, y = "Mean IMV above Binom 1PL (logit) floor") +
    theme_minimal(base_size = 13) + theme(legend.position = "right")
  save_pdf(p, "fig_imv_1pl", w = 8, h = 4)
} else {
  message("  skipping fig_imv_1pl — no imv_1pl data in current files")
}

## ── fig_grm_vs_gpcm ───────────────────────────────────────────────────────
## Head-to-head GRM vs GPCM scatter (RMSE)
message("fig_grm_vs_gpcm")
if (all(c("GRM","GPCM") %in% levels(all_rmse$model))) {
  d_gg <- all_rmse |>
    filter(model %in% c("GRM","GPCM")) |>
    select(dataset, model, rmse) |>
    pivot_wider(names_from = model, values_from = rmse) |>
    left_join(meta |> select(dataset, K, ni, N), by = "dataset") |>
    filter(!is.na(GRM), !is.na(GPCM))
  n_grm_wins  <- sum(d_gg$GRM  < d_gg$GPCM, na.rm = TRUE)
  n_gpcm_wins <- sum(d_gg$GPCM < d_gg$GRM,  na.rm = TRUE)
  p <- ggplot(d_gg, aes(x = GPCM, y = GRM, size = ni)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(colour = "#2166ac", alpha = 0.75) +
    scale_size_continuous(range = c(2, 7), name = "# items") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 3.5,
             label = sprintf("GRM wins: %d    GPCM wins: %d", n_grm_wins, n_gpcm_wins)) +
    labs(x = "CV RMSE — GPCM", y = "CV RMSE — GRM") +
    theme_minimal(base_size = 12)
  save_pdf(p, "fig_grm_vs_gpcm", w = 5, h = 4.5)
} else {
  message("  skipping fig_grm_vs_gpcm — GRM or GPCM not in data")
}

## ── fig_asym ──────────────────────────────────────────────────────────────
## Gap asymmetry s distribution across datasets — only in new-format files
message("fig_asym")
all_asym <- do.call(rbind, Filter(Negate(is.null), lapply(results, `[[`, "asym")))
if (!is.null(all_asym) && nrow(all_asym) > 0) {
  all_asym <- all_asym |> left_join(meta |> select(dataset, K), by = "dataset")
  p <- ggplot(all_asym, aes(x = mean_log_gap_ratio)) +
    geom_histogram(bins = 30, fill = "#2166ac", colour = "white", alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = median(all_asym$mean_log_gap_ratio, na.rm = TRUE),
               colour = "firebrick", linewidth = 0.9) +
    labs(x = "Mean log gap-ratio (asymmetry s) per dataset",
         y = "Count",
         caption = sprintf("Median = %.2f  |  n = %d datasets",
                           median(all_asym$mean_log_gap_ratio, na.rm = TRUE),
                           nrow(all_asym))) +
    theme_minimal(base_size = 13)
  save_pdf(p, "fig_asym", w = 6, h = 4)
} else {
  message("  skipping fig_asym — no asym data in current files")
}

## ═══════════════════════════════════════════════════════════════════════════
## SI FIGURE: STEP 1 OF ISOERF INJECTIVITY PROOF (identification of a)
## ═══════════════════════════════════════════════════════════════════════════
message("\nstep1_isoerf")

# GPCM expected response function for K=3 (two thresholds)
gpcm_erf <- function(theta, a, b) {
  # b is a numeric vector of length K-1
  K <- length(b) + 1L
  # log-numerators: sum_{j=1}^{k} a*(theta - b_j), empty sum = 0
  ln <- matrix(0, nrow = length(theta), ncol = K)
  for (k in seq_len(K - 1L)) {
    ln[, k + 1L] <- ln[, k] + a * (theta - b[k])
  }
  # softmax (numerically stable)
  ln <- ln - apply(ln, 1, max)
  ex <- exp(ln)
  pr <- ex / rowSums(ex)
  # expected score
  as.numeric(pr %*% 0:(K - 1L))
}

# Numerical derivative of log ERF via central differences
d_log_erf <- function(theta, a, b, h = 0.025) {
  ep <- gpcm_erf(theta + h, a, b)
  em <- gpcm_erf(theta - h, a, b)
  ifelse(ep > 1e-12 & em > 1e-12, (log(ep) - log(em)) / (2 * h), NA_real_)
}

theta_s1 <- seq(-5, 4, length.out = 300)

# Three threshold configurations, all with the same a
cfgs_s1 <- list(
  list(b = c(-0.5,  0.5), label = "b = (-0.5, 0.5)"),
  list(b = c( 0.0,  1.0), label = "b = (0.0, 1.0)"),
  list(b = c( 0.5,  1.5), label = "b = (0.5, 1.5)")
)

cols_s1 <- c("#378ADD", "#1D9E75", "#D85A30")

## ── Panel A: log E[X|theta], a = 1 ──────────────────────────────────────────

df_s1_A <- lapply(seq_along(cfgs_s1), function(i) {
  cfg <- cfgs_s1[[i]]
  erf_vals <- gpcm_erf(theta_s1, a = 1, b = cfg$b)
  log_erf  <- ifelse(erf_vals > 1e-12, log(erf_vals), NA_real_)
  # asymptote: a*theta - a*S_1 = theta - b[1]
  asym <- theta_s1 - cfg$b[1]
  data.frame(
    theta    = theta_s1,
    log_erf  = log_erf,
    asym     = asym,
    config   = cfg$label,
    col      = cols_s1[i]
  )
}) |> bind_rows()

p_s1_A <- ggplot(df_s1_A, aes(x = theta, group = config, colour = config)) +
  geom_line(aes(y = log_erf), linewidth = 0.9) +
  geom_line(aes(y = asym), linewidth = 0.5, linetype = "dashed") +
  scale_colour_manual(values = setNames(cols_s1, sapply(cfgs_s1, `[[`, "label"))) +
  coord_cartesian(ylim = c(-8, 1.5)) +
  labs(x = expression(theta), y = expression(log~E*"["*X*"|"*theta*"]"),
       title = "A") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 8),
        panel.grid.minor = element_blank(),
        plot.title      = element_text(face = "bold"))

## ── Panel B: d/dtheta log E[X|theta], a = 1 ─────────────────────────────────

df_s1_B <- lapply(seq_along(cfgs_s1), function(i) {
  cfg <- cfgs_s1[[i]]
  data.frame(
    theta  = theta_s1,
    deriv  = d_log_erf(theta_s1, a = 1, b = cfg$b),
    config = cfg$label,
    col    = cols_s1[i]
  )
}) |> bind_rows()

p_s1_B <- ggplot(df_s1_B, aes(x = theta, y = deriv, group = config, colour = config)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50",
             linewidth = 0.7) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  annotate("text", x = -4.8, y = 1.08, label = "a = 1",
           size = 3, colour = "grey40", hjust = 0) +
  scale_colour_manual(values = setNames(cols_s1, sapply(cfgs_s1, `[[`, "label"))) +
  coord_cartesian(ylim = c(0, 2.2)) +
  labs(x = expression(theta),
       y = expression(frac(d, d*theta)~log~E*"["*X*"|"*theta*"]"),
       title = "B") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

## ── Panel C: d/dtheta log E[X|theta], a = 1 vs a = 2 ────────────────────────

cfgs_s1_C <- list(
  list(a = 1, b = c(0.0, 1.0), label = "a = 1, b = (0, 1)",    lty = "solid"),
  list(a = 1, b = c(0.5, 1.5), label = "a = 1, b = (0.5, 1.5)", lty = "solid"),
  list(a = 2, b = c(0.0, 1.0), label = "a = 2, b = (0, 1)",    lty = "dashed"),
  list(a = 2, b = c(0.5, 1.5), label = "a = 2, b = (0.5, 1.5)", lty = "dashed")
)
cols_s1_C <- c("#378ADD", "#1D9E75", "#378ADD", "#1D9E75")

df_s1_C <- lapply(seq_along(cfgs_s1_C), function(i) {
  cfg <- cfgs_s1_C[[i]]
  data.frame(
    theta   = theta_s1,
    deriv   = d_log_erf(theta_s1, a = cfg$a, b = cfg$b),
    config  = cfg$label,
    lty     = cfg$lty,
    col     = cols_s1_C[i]
  )
}) |> bind_rows()

p_s1_C <- ggplot(df_s1_C, aes(x = theta, y = deriv,
                        group = config, colour = config,
                        linetype = lty)) +
  geom_hline(yintercept = c(1, 2), linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  annotate("text", x = -4.8, y = 1.08, label = "a = 1",
           size = 3, colour = "grey40", hjust = 0) +
  annotate("text", x = -4.8, y = 2.08, label = "a = 2",
           size = 3, colour = "grey40", hjust = 0) +
  scale_colour_manual(
    values = setNames(cols_s1_C, sapply(cfgs_s1_C, `[[`, "label"))
  ) +
  scale_linetype_identity() +
  coord_cartesian(ylim = c(0, 2.8)) +
  labs(x = expression(theta),
       y = expression(frac(d, d*theta)~log~E*"["*X*"|"*theta*"]"),
       title = "C") +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 7.5),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold")) +
  guides(colour = guide_legend(nrow = 2))

fig_s1 <- p_s1_A + p_s1_B + p_s1_C + plot_layout(ncol = 3)

ggsave(file.path(OUT, "step1_isoerf.pdf"), fig_s1, width = 10, height = 3.8)
ggsave(file.path(OUT, "step1_isoerf.png"), fig_s1, width = 10, height = 3.8, dpi = 300)
message("  saved step1_isoerf.pdf / step1_isoerf.png")

message("\nDone. All figures saved to ", OUT)
