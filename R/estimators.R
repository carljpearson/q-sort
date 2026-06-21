# -----------------------------------------------------------------------------
# estimators.R - Bucket-code estimators and their honest-interval check
#
#   est_naive_mean       the status quo: mean of the 5-4-3-2-1 codes, as shares
#   est_naive_alt        same idea, steeper coding vector (equal-spacing stress-test)
#   est_naive_bootstrap  same point estimate, respondent-level cluster bootstrap CI
#   est_dirichlet        compositional comparator (Dirichlet-multinomial)
#
# Every estimator returns the same tidy shape:
#   item | estimate | conf_low | conf_high   (all on the share scale, sum ~ 1)
# -----------------------------------------------------------------------------

standard_codes <- c(`5` = 5, `4` = 4, `3` = 3, `2` = 2, `1` = 1)
steep_codes    <- c(`5` = 10, `4` = 5, `3` = 3, `2` = 1, `1` = 0)

# Turn bands into per-respondent shares under a given coding vector. Each
# respondent's codes sum to the same constant, so a respondent's share for an
# item is just its code over that constant - a clean per-respondent quantity.
score_shares <- function(qsort, codes = standard_codes) {
  qsort %>%
    mutate(code = codes[as.character(band)]) %>%
    group_by(respondent) %>%
    mutate(share = code / sum(code)) %>%
    ungroup()
}

# Point estimate: average each item's share across respondents.
point_shares <- function(qsort, codes = standard_codes) {
  score_shares(qsort, codes) %>%
    group_by(item) %>%
    summarise(estimate = mean(share), .groups = "drop")
}

# Naive mean with a respondent-level normal-approximation interval. The interval
# is correctly sized; it under-covers because the point estimate is compressed.
est_naive_mean <- function(qsort, codes = standard_codes, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  score_shares(qsort, codes) %>%
    group_by(item) %>%
    summarise(
      estimate = mean(share),
      se = sd(share) / sqrt(n()),
      .groups = "drop"
    ) %>%
    transmute(
      item, estimate,
      conf_low = estimate - z * se,
      conf_high = estimate + z * se
    )
}

# Same machinery, steeper coding vector - isolates the equal-spacing assumption.
est_naive_alt <- function(qsort, codes = steep_codes, level = 0.95) {
  est_naive_mean(qsort, codes = codes, level = level)
}

# Respondent-level (cluster) bootstrap interval around the naive point estimate.
# Resamples whole respondents, never individual item-responses.
est_naive_bootstrap <- function(qsort, codes = standard_codes, level = 0.95, times = 1000) {
  point <- point_shares(qsort, codes)
  ci <- cluster_bootstrap_ci(qsort, function(d) point_shares(d, codes), level, times)
  left_join(point, ci, by = "item")
}

# Compositional comparator: model each respondent's share vector as a Dirichlet.
# Standard codes (no zeros) keep the shares strictly inside the simplex.
est_dirichlet <- function(qsort, level = 0.95, times = 500) {
  point <- point_dirichlet(qsort)
  ci <- cluster_bootstrap_ci(qsort, point_dirichlet, level, times)
  left_join(point, ci, by = "item")
}

point_dirichlet <- function(qsort) {
  wide <- score_shares(qsort) %>%
    select(respondent, item, share) %>%
    pivot_wider(names_from = item, values_from = share, names_sort = TRUE) %>%
    arrange(respondent) %>%
    select(-respondent)

  response <- DirichletReg::DR_data(as.matrix(wide))
  model <- DirichletReg::DirichReg(response ~ 1)

  tibble(
    item = as.integer(colnames(wide)),
    estimate = as.numeric(fitted(model)[1, ])
  )
}

# Generic respondent-level bootstrap: resample whole respondents, recompute the
# point estimate, return percentile intervals per item. `statistic` is a plain
# function of a qsort data frame returning columns item + estimate.
cluster_bootstrap_ci <- function(qsort, statistic, level = 0.95, times = 1000) {
  alpha <- 1 - level
  rsample::group_bootstraps(qsort, group = respondent, times = times)$splits %>%
    map(~ statistic(rsample::analysis(.x))) %>%
    list_rbind() %>%
    group_by(item) %>%
    summarise(
      conf_low = quantile(estimate, alpha / 2, names = FALSE),
      conf_high = quantile(estimate, 1 - alpha / 2, names = FALSE),
      .groups = "drop"
    )
}
