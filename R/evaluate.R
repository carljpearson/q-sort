# -----------------------------------------------------------------------------
# evaluate.R - put estimates on the share scale and score them against truth
#
# Every estimator already returns shares with intervals. Here we:
#   - turn log-worths into shares + delta-method intervals (softmax_shares_ci)
#   - run all estimators on one simulated sample (run_all_estimators)
#   - label items as "extreme" (true top/bottom) or "middle"
#   - score each estimate for bias, squared error, interval coverage and width
# -----------------------------------------------------------------------------

# Delta-method shares + intervals from log-worths. `beta` includes the reference
# item fixed at 0; `vcov_free` is the covariance of the free (non-reference)
# coefficients. Maps both through the softmax.
softmax_shares_ci <- function(beta, vcov_free, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  shares <- exp(beta) / sum(exp(beta))

  # Jacobian of the softmax; drop the reference column (its beta is fixed at 0).
  jacobian <- -outer(shares, shares)
  diag(jacobian) <- shares * (1 - shares)
  jacobian_free <- jacobian[, -1, drop = FALSE]

  variance <- diag(jacobian_free %*% vcov_free %*% t(jacobian_free))
  se <- sqrt(pmax(variance, 0))

  tibble(
    estimate = as.numeric(shares),
    conf_low = shares - z * se,
    conf_high = shares + z * se
  )
}

# Run every estimator on one simulated sample.
run_all_estimators <- function(qsort) {
  bind_rows(
    est_naive_mean(qsort)      %>% mutate(estimator = "naive_mean"),
    est_naive_alt(qsort)       %>% mutate(estimator = "naive_alt_coding"),
    est_naive_bootstrap(qsort) %>% mutate(estimator = "naive_bootstrap"),
    est_pl_ranking(qsort)      %>% mutate(estimator = "pl_ranking"),
    est_pl_trials(qsort)       %>% mutate(estimator = "pl_trials"),
    est_dirichlet(qsort)       %>% mutate(estimator = "dirichlet")
  ) %>%
    relocate(estimator)
}

# Label the true top/bottom items as "extreme", the rest as "middle".
label_strata <- function(truth, prop_extreme = 0.2) {
  n <- nrow(truth)
  k <- max(1, round(n * prop_extreme))
  truth %>%
    arrange(desc(true_share)) %>%
    mutate(
      rank = row_number(),
      stratum = if_else(rank <= k | rank > n - k, "extreme", "middle")
    ) %>%
    select(item, true_share, stratum)
}

# Score every estimate against the planted truth.
evaluate_estimates <- function(estimates, truth, prop_extreme = 0.2) {
  estimates %>%
    left_join(label_strata(truth, prop_extreme), by = "item") %>%
    mutate(
      error = estimate - true_share,
      covered = as.integer(conf_low <= true_share & true_share <= conf_high),
      width = conf_high - conf_low
    )
}
