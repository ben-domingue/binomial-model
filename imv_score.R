# src/imv_score.R
#
# Polytomous IMVs for the paper, computed with the imv package (>= 0.4.0;
# ben-domingue/imv main). Every omega_c / omega_t in the paper goes through
# imv_by_item(); no script defines its own imv_c / imv_t.
#
# Method (adopted 2026-11-11):
#   * imv::imv_c() / imv::imv_t() with their defaults: weights = "continuity",
#     sigma = 1e-4.
#   * One IMV per item on that item's held-out responses, combined as a mean
#     weighted by the number of held-out responses per item. This is what
#     imv::imv() does for polytomous mirt models (aggregate = "items",
#     item.weights = "responses"); we cannot call imv() itself because most of
#     our models are not mirt fits.
#
# A task in which a model predicts worse than a fair coin is floored at w = 0.5
# by the package (with a warning). We muffle those warnings and count, per
# index, the items with at least one floored task (attributes floor_c, floor_t).

stopifnot(utils::packageVersion("imv") >= "0.4.0")

# x     held-out responses coded 0..K-1
# item  item id for each held-out response
# P0    baseline category probabilities (n x K matrix)
# P1    enhanced category probabilities (n x K matrix)
# Returns c(imv_c = , imv_t = ) with attributes floor_c, floor_t, n_items.
imv_by_item <- function(x, item, P0, P1, K) {
  stopifnot(length(x) == length(item), nrow(P0) == length(x), nrow(P1) == length(x))
  items <- unique(item)
  one <- function(f, i) {
    floored <- FALSE
    v <- withCallingHandlers(
      as.numeric(f(x[i], P0[i, , drop = FALSE], P1[i, , drop = FALSE], K = K)),
      warning = function(w) {
        if (grepl("fair coin", conditionMessage(w))) {
          floored <<- TRUE
          invokeRestart("muffleWarning")
        }
      })
    c(v, floored)
  }
  res <- lapply(items, function(j) {
    i <- item == j
    c(n = sum(i), one(imv::imv_c, i), one(imv::imv_t, i))
  })
  m <- do.call(rbind, res)        # n, c, c_floor, t, t_floor
  wmean <- function(v, n) { ok <- !is.na(v) & n > 0; if (!any(ok)) NA_real_ else sum(n[ok] * v[ok]) / sum(n[ok]) }
  structure(c(imv_c = wmean(m[, 2], m[, 1]), imv_t = wmean(m[, 4], m[, 1])),
            floor_c = sum(m[, 3]), floor_t = sum(m[, 5]), n_items = length(items))
}

# Per-item training-prevalence baseline: row r gets the category proportions of
# item[r] among the training responses (train_x, train_item). Used where the
# paper compares a model to a "prevalence" baseline (imv::imv() baseline = "item").
prevalence_matrix <- function(item, train_x, train_item, K) {
  tab <- t(vapply(split(train_x, factor(train_item, levels = unique(c(item, train_item)))),
                  function(v) tabulate(v + 1L, K) / max(length(v), 1L), numeric(K)))
  tab[as.character(item), , drop = FALSE]
}
