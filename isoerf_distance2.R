library(ggplot2)
library(patchwork)

logis <- function(x) 1 / (1 + exp(-x))

# ---- ERF and CRF for each model (K=3) ----

erf_grm <- function(theta, a, b) {
  rowSums(outer(theta, b, function(t, bk) logis(a * (t - bk))))
}
crf_grm <- function(theta, a, b) {
  p0 <- logis(a * (theta - b[1]))
  p1 <- logis(a * (theta - b[2]))
  data.frame(theta = theta, k0 = 1 - p0, k1 = p0 - p1, k2 = p1)
}

erf_gpcm <- function(theta, a, b) {
  sapply(theta, function(t) {
    u1 <- exp(a * (t - b[1]))
    u2 <- exp(a * (t - b[1]) + a * (t - b[2]))
    (u1 + 2 * u2) / (1 + u1 + u2)
  })
}
crf_gpcm <- function(theta, a, b) {
  m <- t(sapply(theta, function(t) {
    u1 <- exp(a * (t - b[1]))
    u2 <- exp(a * (t - b[1]) + a * (t - b[2]))
    d  <- 1 + u1 + u2
    c(k0 = 1/d, k1 = u1/d, k2 = u2/d)
  }))
  as.data.frame(cbind(theta = theta, m))
}

erf_binom <- function(theta, a, b) {
  2 * logis(a * (theta - b[1]))   # b[2] not used
}
crf_binom <- function(theta, a, b) {
  p <- logis(a * (theta - b[1]))
  data.frame(theta = theta, k0 = (1-p)^2, k1 = 2*p*(1-p), k2 = p^2)
}

# ---- Per-model references ----
refs <- list(
  GRM      = list(a = 1, b = c(0, 1),  erf_fn = erf_grm,  crf_fn = crf_grm),
  GPCM     = list(a = 1, b = c(0, 1),  erf_fn = erf_gpcm, crf_fn = crf_gpcm),
  Binomial = list(a = 1, b = c(0.5, NA), erf_fn = erf_binom, crf_fn = crf_binom)
)

erf_dist <- function(erf_ref_fn, erf_fn, a_prime, b_prime) {
  integrate(function(theta) {
    (erf_ref_fn(theta) - erf_fn(theta, a_prime, b_prime))^2 * dnorm(theta)
  }, lower = -5, upper = 5, subdivisions = 200)$value
}

# ---- Grids ----
b1_vals <- seq(-2, 2, length.out = 60)
b2_vals <- seq(-1, 3, length.out = 60)

grids <- list(
  GRM     = subset(expand.grid(b1 = b1_vals, b2 = b2_vals), b1 < b2),
  GPCM    = expand.grid(b1 = b1_vals, b2 = b2_vals),
  Binomial = expand.grid(b1 = b1_vals, b2 = b2_vals)
)

erfs <- list(GRM = erf_grm, GPCM = erf_gpcm, Binomial = erf_binom)
crfs <- list(GRM = crf_grm, GPCM = crf_gpcm, Binomial = crf_binom)

# ---- Plot helpers ----
make_dist_plot <- function(g, erf_ref_fn, b_ref_row, erf_fn, ap,
                           model_name, show_legend = FALSE, highlight = NULL,
                           scale_max = NA, show_isoerf_contour = FALSE) {
  if (!"dist" %in% names(g)) {
    g$dist <- mapply(function(b1, b2) erf_dist(erf_ref_fn, erf_fn, ap, c(b1, b2)),
                     g$b1, g$b2)
    cat(sprintf("%-8s a'=%.1f  min=%.6f  at b1=%.3f b2=%.3f\n",
                model_name, ap, min(g$dist),
                g$b1[which.min(g$dist)], g$b2[which.min(g$dist)]))
  }
  p <- ggplot(g, aes(x = b1, y = b2, fill = sqrt(dist))) +
    geom_tile() +
    { if (!is.na(b_ref_row[2]))
        annotate("point", x = b_ref_row[1], y = b_ref_row[2],
                 shape = 4, size = 3, color = "black")
      else list() } +
    scale_fill_gradient(low = "white", high = "red",
                        limits = c(0, scale_max),
                        name = expression(sqrt(dist)),
                        guide = if (show_legend) "colourbar" else "none") +
    labs(title = sprintf("%s  a'=%.1f", model_name, ap),
         x = expression(b[1]*"'"), y = expression(b[2]*"'")) +
    theme_minimal()
  if (show_isoerf_contour) {
    contour_level <- quantile(sqrt(g$dist), probs = 0.03)
    p <- p + geom_contour(data = g, aes(x = b1, y = b2, z = sqrt(dist)),
                          inherit.aes = FALSE, breaks = contour_level,
                          color = "steelblue", linewidth = 0.8, linetype = "dashed")
  }
  if (!is.null(highlight)) {
    p <- p +
      annotate("point", x = highlight[1], y = highlight[2],
               shape = 21, size = 4, color = "black", fill = "white", stroke = 1.2) +
      annotate("text", x = highlight[1] + 0.15, y = highlight[2] - 0.15,
               label = sprintf("(%.2f, %.2f)", highlight[1], highlight[2]),
               size = 2.8, hjust = 0)
  }
  list(plot = p, grid = g)
}

make_erf_crf <- function(ref, erf_fn, crf_fn, a_cand, b_cand) {
  theta_seq <- seq(-4, 4, length.out = 300)
  erf_ref_fn <- function(theta) ref$erf_fn(theta, ref$a, ref$b)
  rmise <- sqrt(integrate(function(theta)
    (erf_ref_fn(theta) - erf_fn(theta, a_cand, b_cand))^2 * dnorm(theta),
    lower = -5, upper = 5, subdivisions = 200)$value)

  erf_df <- data.frame(
    theta = rep(theta_seq, 2),
    erf   = c(erf_ref_fn(theta_seq), erf_fn(theta_seq, a_cand, b_cand)),
    model = rep(c("Reference", "Candidate"), each = length(theta_seq))
  )
  p_erf <- ggplot(erf_df, aes(x = theta, y = erf, color = model, linetype = model)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("black", "red"), guide = "none") +
    scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
    labs(x = expression(theta), y = expression(E*"["*X*"|"*theta*"]"),
         title = "ERF",
         subtitle = {
           fmt_b <- function(b) paste(sprintf("%.2f", b[!is.na(b)]), collapse = ", ")
           sprintf("-- Ref (a=%.1f, b=(%s))    -- Cand (a'=%.1f, b'=(%s))    RMISE=%.4f",
                   ref$a, fmt_b(ref$b), a_cand, fmt_b(b_cand), rmise)
         }) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 7))

  r  <- ref$crf_fn(theta_seq, ref$a, ref$b)
  cc <- crf_fn(theta_seq, a_cand, b_cand)
  crf_df <- rbind(
    data.frame(theta = theta_seq, prob = r$k0,  k = "k=0", model = "Reference"),
    data.frame(theta = theta_seq, prob = r$k1,  k = "k=1", model = "Reference"),
    data.frame(theta = theta_seq, prob = r$k2,  k = "k=2", model = "Reference"),
    data.frame(theta = theta_seq, prob = cc$k0, k = "k=0", model = "Candidate"),
    data.frame(theta = theta_seq, prob = cc$k1, k = "k=1", model = "Candidate"),
    data.frame(theta = theta_seq, prob = cc$k2, k = "k=2", model = "Candidate")
  )
  p_crf <- ggplot(crf_df, aes(x = theta, y = prob, color = k,
                               linetype = model, linewidth = model)) +
    geom_line() +
    scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
    scale_linewidth_manual(values = c(0.6, 0.9), guide = "none") +
    scale_color_brewer(palette = "Set1", name = "Category") +
    labs(x = expression(theta), y = expression("P(X=k|"*theta*")"),
         title = "CRF", subtitle = "-- Reference    - - Candidate") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 7))

  list(erf = p_erf, crf = p_crf)
}

make_dist_plot_1d <- function(b1_vals, erf_ref_fn, b_ref, erf_fn, ap,
                              model_name, show_legend = FALSE, highlight = NULL,
                              scale_max = NA, precomputed_g = NULL) {
  if (is.null(precomputed_g)) {
    g <- data.frame(b1 = b1_vals)
    g$dist <- sapply(b1_vals, function(b1) erf_dist(erf_ref_fn, erf_fn, ap, c(b1, NA)))
    cat(sprintf("%-8s a'=%.1f  min=%.6f  at b1=%.3f\n",
                model_name, ap, min(g$dist), g$b1[which.min(g$dist)]))
  } else {
    g <- precomputed_g
  }
  # reference a=1; candidate a = ap (0.7 or 1.5); only b1' varies on x-axis
  p <- ggplot(g, aes(x = b1, y = 0, fill = sqrt(dist))) +
    geom_tile(height = 1) +
    annotate("point", x = b_ref, y = 0, shape = 4, size = 4, color = "black", stroke = 1.5) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, scale_max),
                        name = expression(sqrt(dist)),
                        guide = if (show_legend) "colourbar" else "none") +
    scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
    labs(title = sprintf("%s  a'=%.1f", model_name, ap),
         x = expression(b[1]*"'"), y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank())
  if (!is.null(highlight)) {
    p <- p +
      annotate("point", x = highlight, y = 0,
               shape = 21, size = 5, color = "black", fill = "white", stroke = 1.2) +
      annotate("text", x = highlight + 0.1, y = 0.3,
               label = sprintf("(%.2f)", highlight), size = 2.8, hjust = 0)
  }
  list(plot = p, grid = g)
}

# ---- Pass 1: compute all distances, cache grids ----
a_prime <- 1.5
cached <- list()

for (model in c("GRM", "GPCM", "Binomial")) {
  ref <- refs[[model]]
  efn <- erfs[[model]]
  erf_ref_fn <- local({ r <- ref; function(theta) r$erf_fn(theta, r$a, r$b) })

  if (model == "Binomial") {
    g07 <- data.frame(b1 = b1_vals)
    g07$dist <- sapply(b1_vals, function(b1) erf_dist(erf_ref_fn, efn, 0.7, c(b1, NA)))
    g15 <- data.frame(b1 = b1_vals)
    g15$dist <- sapply(b1_vals, function(b1) erf_dist(erf_ref_fn, efn, 1.5, c(b1, NA)))
  } else {
    g <- grids[[model]]
    g07 <- g
    g07$dist <- mapply(function(b1, b2) erf_dist(erf_ref_fn, efn, 0.7, c(b1, b2)), g$b1, g$b2)
    g15 <- g
    g15$dist <- mapply(function(b1, b2) erf_dist(erf_ref_fn, efn, 1.5, c(b1, b2)), g$b1, g$b2)
  }
  cached[[model]] <- list(g07 = g07, g15 = g15, erf_ref_fn = erf_ref_fn)
  cat(sprintf("Computed %s\n", model))
}

scale_max <- max(sqrt(unlist(lapply(cached, function(x) c(x$g07$dist, x$g15$dist)))))
cat(sprintf("Global scale max (sqrt dist): %.4f\n", scale_max))

# ---- Pass 2: build plots with shared color scale ----
rows <- list()

for (model in c("GRM", "GPCM", "Binomial")) {
  ref <- refs[[model]]
  efn <- erfs[[model]]
  cfn <- crfs[[model]]
  cc  <- cached[[model]]

  if (model == "Binomial") {
    r07_plot <- make_dist_plot_1d(b1_vals, cc$erf_ref_fn, ref$b[1], efn, 0.7, model,
                                  show_legend = FALSE, scale_max = scale_max,
                                  precomputed_g = cc$g07)$plot
    cand_b1  <- cc$g15$b1[which.min(cc$g15$dist)]
    r15_plot <- make_dist_plot_1d(b1_vals, cc$erf_ref_fn, ref$b[1], efn, 1.5, model,
                                  show_legend = TRUE, scale_max = scale_max,
                                  precomputed_g = cc$g15, highlight = cand_b1)$plot
    cand <- c(cand_b1, NA)
  } else {
    iso_contour <- (model == "GRM")
    r07_plot <- make_dist_plot(cc$g07, cc$erf_ref_fn, ref$b, efn, 0.7, model,
                               show_legend = FALSE, scale_max = scale_max,
                               show_isoerf_contour = iso_contour)$plot
    cand     <- with(cc$g15, c(b1[which.min(dist)], b2[which.min(dist)]))
    r15_plot <- make_dist_plot(cc$g15, cc$erf_ref_fn, ref$b, efn, 1.5, model,
                               show_legend = TRUE, scale_max = scale_max,
                               highlight = cand, show_isoerf_contour = iso_contour)$plot
  }

  panels <- make_erf_crf(ref, efn, cfn, a_prime, cand)
  rows[[model]] <- (r07_plot + r15_plot) | (panels$erf / panels$crf)
}

p_final <- (rows$GRM / rows$GPCM / rows$Binomial) +
  plot_annotation(
    title    = "Within-model isoERF distance  (GRM / GPCM / Binomial,  K=3,  reference a=1)",
    subtitle = "Left column: candidate a'=0.7    Right column: candidate a'=1.5    x = reference    o = minimum-distance candidate"
  )

ggsave("/home/ben/Dropbox/Apps/Overleaf/binomial_model/src/fig_isoerf_distance2.pdf",
       p_final, width = 13, height = 16)
cat("saved\n")
