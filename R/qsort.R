# =============================================================================
# qsort.R - the whole Q-sort estimator study in one readable place
#
# This file replaces the old grab-bag of one-off helpers (dgp.R, encode_*.R,
# estimators.R, evaluate.R) with a small set of functions you can read top to
# bottom:
#
#   make_truth()        plant a known importance profile
#   simulate_qsort()    walk respondents through the four-pass Q-sort
#   fit_*()             the seven estimators, each returning item|estimate|se
#   score_estimates()   compare any estimate against the planted truth
#
# Every estimator returns the same three columns - item, estimate (a share that
# sums to ~1), and se (a standard error on that share scale). Confidence
# intervals are built downstream as estimate +/- z * se, so a single scoring
# step can check coverage at any nominal level. That uniform contract is what
# lets the calibration curve - the point of the study - fall out in a few lines.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(survival)  # clogit() and its strata() special, used by fit_pl_trials
})

# -----------------------------------------------------------------------------
# 1. Truth and data-generating process
# -----------------------------------------------------------------------------
#
# The five Q-sort bands, top to bottom, are the standard 5-4-3-2-1 codes:
#   5 most (1 item) | 4 next-most (2) | 3 middle (the rest) | 2 next-least (2) |
#   1 least (1 item). Needs >= 7 items so the middle band is non-empty.

# True item worths for a list of `n_items`, decaying as a power law. `lumpiness`
# sets how concentrated importance is: 0 is flat, 1 is Zipf-like (1, 1/2, 1/3),
# and >1 lets one or two items run away with it. Shares sum to one.
make_truth <- function(n_items, lumpiness) {
  tibble(item = seq_len(n_items)) %>%
    mutate(
      worth = item^(-lumpiness),
      true_share = worth / sum(worth)
    )
}

# Simulate a whole sample of respondents sorting the same list. Each respondent
# makes four sequential picks - most, next two most, least, next two least -
# drawing without replacement with probability proportional to worth on the
# forward scale (the "most" passes) or 1/worth on the reversed scale (the
# "least" passes). Whatever is never picked is the tied middle. `noise` is a
# temperature on the log-worths: 1 is the nominal model, larger flattens the
# choice probabilities (noisier respondents), smaller sharpens them.
simulate_qsort <- function(truth, n_resp, noise = 1) {
  items   <- truth$item
  forward <- truth$worth^(1 / noise)        # max-type scale, "most" passes
  reverse <- (1 / truth$worth)^(1 / noise)  # min-type scale, "least" passes

  # Draw `k` items out of the still-available `pool`, weighted by `w`. Indexing
  # through `length(pool)` sidesteps sample()'s single-value gotcha.
  draw <- function(pool, w, k) pool[sample.int(length(pool), k, prob = w[pool])]

  one_respondent <- function(r) {
    band <- integer(length(items))
    pool <- items

    most <- draw(pool, forward, 1); band[most] <- 5; pool <- setdiff(pool, most)
    nmost <- draw(pool, forward, 2); band[nmost] <- 4; pool <- setdiff(pool, nmost)
    least <- draw(pool, reverse, 1); band[least] <- 1; pool <- setdiff(pool, least)
    nleast <- draw(pool, reverse, 2); band[nleast] <- 2; pool <- setdiff(pool, nleast)
    band[pool] <- 3  # the tied middle

    tibble(respondent = r, item = items, band = band)
  }

  map(seq_len(n_resp), one_respondent) %>% list_rbind()
}

# -----------------------------------------------------------------------------
# 2. Shared scaffolding the estimators reuse
# -----------------------------------------------------------------------------

# Turn bands into per-respondent shares under a coding vector. Each respondent's
# codes sum to a constant, so a respondent's share for an item is just its code
# over that total - the natural per-respondent compositional quantity.
CODES_STANDARD <- c(`5` = 5, `4` = 4, `3` = 3, `2` = 2, `1` = 1)
CODES_STEEP    <- c(`5` = 10, `4` = 5, `3` = 3, `2` = 1, `1` = 0)

respondent_shares <- function(qsort, codes) {
  qsort %>%
    mutate(code = codes[as.character(band)]) %>%
    group_by(respondent) %>%
    mutate(share = code / sum(code)) %>%
    ungroup()
}

# The per-respondent share matrix (respondents x items) under a coding vector -
# the natural unit for a respondent-level bootstrap, since resampling its rows
# resamples whole people.
share_matrix <- function(qsort, codes) {
  respondent_shares(qsort, codes) %>%
    pivot_wider(id_cols = respondent, names_from = item,
                values_from = share, names_sort = TRUE) %>%
    select(-respondent) %>%
    as.matrix()
}

# Map a coefficient vector to shares via the softmax and propagate its
# covariance to a per-share standard error with the delta method. `beta` and the
# full covariance `vcov` cover the same entries (a fixed reference item simply
# carries a zero row and column). The softmax Jacobian's rows sum to zero, so
# the over-parametrised case is handled without dropping anything by hand.
softmax_with_se <- function(beta, vcov) {
  shares <- exp(beta) / sum(exp(beta))
  jacobian <- diag(shares) - outer(shares, shares)
  variance <- diag(jacobian %*% vcov %*% t(jacobian))
  tibble(estimate = as.numeric(shares), se = sqrt(pmax(variance, 0)))
}

# -----------------------------------------------------------------------------
# 3. The estimators - each returns item | estimate | se on the share scale
# -----------------------------------------------------------------------------

# (a) Naive mean of the bucket codes - the status quo. The SE is the
# respondent-level spread of each item's share; because every respondent
# contributes one share, sd/sqrt(n) is already a cluster-correct SE. It is
# honestly sized - it just sits around a compressed point estimate.
fit_naive <- function(qsort, codes = CODES_STANDARD) {
  respondent_shares(qsort, codes) %>%
    group_by(item) %>%
    summarise(estimate = mean(share), se = sd(share) / sqrt(n()), .groups = "drop")
}

# (b) Same mean with a steeper coding vector - isolates the equal-spacing
# assumption. If the answer moves just because the arbitrary codes changed, the
# even-spacing claim was unearned.
fit_naive_steep <- function(qsort) fit_naive(qsort, codes = CODES_STEEP)

# (c) Naive point estimate, but with a cluster-bootstrap SE - the naive mean's
# fairest shot at an honest interval. We resample whole respondents (rows of the
# share matrix) with replacement and take the spread of the resampled means. If
# it still under-covers, "we sized the interval wrong" is ruled out and the bias
# is the culprit. (In practice it lands close to the analytic SE above, which is
# the point: the trouble is the estimate, not the interval.)
fit_naive_boot <- function(qsort, times = 400) {
  M <- share_matrix(qsort, CODES_STANDARD)
  n <- nrow(M)
  boot_means <- replicate(times, colMeans(M[sample.int(n, n, replace = TRUE), , drop = FALSE]))
  tibble(
    item = as.integer(colnames(M)),
    estimate = colMeans(M),
    se = apply(boot_means, 1, sd)
  )
}

# (d) Dirichlet-multinomial comparator: model each respondent's share vector as
# a draw from one Dirichlet and read off the fitted mean composition. With an
# intercept-only fit the mean shares are exactly the softmax of the per-item
# log-alpha intercepts, so a single fit plus the delta method gives the SE - no
# bootstrap needed. Standard (non-zero) codes keep every share inside the
# simplex. Sits between the naive mean and the choice models.
fit_dirichlet <- function(qsort) {
  wide <- respondent_shares(qsort, CODES_STANDARD) %>%
    pivot_wider(id_cols = respondent, names_from = item,
                values_from = share, names_sort = TRUE) %>%
    select(-respondent)

  # DirichReg uses non-standard evaluation: the response must be a named object
  # in the formula, not an inline expression.
  response <- DirichletReg::DR_data(as.matrix(wide))
  model <- DirichletReg::DirichReg(response ~ 1)

  softmax_with_se(beta = unname(unlist(coef(model))), vcov = vcov(model)) %>%
    mutate(item = as.integer(colnames(wide)), .before = 1)
}

# (e) Plackett-Luce on the single partial ranking the four passes imply
# (most > next most > tied middle > next least > least). PlackettLuce returns
# worths already normalised to the share scale with a covariance matrix, so the
# SE is just its square root. NOTE: the tied middle is a tie of order
# (n_items - 6); PlackettLuce's likelihood cost grows steeply with that order,
# so this estimator is only practical for short lists (see fit_pl_ranking_feasible).
fit_pl_ranking <- function(qsort) {
  ranks <- qsort %>%
    transmute(respondent, item, rank = 6L - band) %>%
    pivot_wider(names_from = item, values_from = rank, names_sort = TRUE) %>%
    arrange(respondent) %>%
    select(-respondent) %>%
    as.matrix()

  model <- PlackettLuce::PlackettLuce(PlackettLuce::as.rankings(ranks))
  worths <- PlackettLuce::itempar(model, vcov = TRUE)
  ord <- order(as.integer(names(worths)))

  tibble(
    item = as.integer(names(worths))[ord],
    estimate = as.numeric(worths)[ord],
    se = sqrt(diag(attr(worths, "vcov")))[ord]
  )
}

# A single switch for whether the single-ranking fit is worth attempting: the
# middle band (n_items - 6 items) becomes a tie of that order, and the fit is
# only tractable while it stays small.
fit_pl_ranking_feasible <- function(n_items, max_tie = 8) (n_items - 6) <= max_tie

# (f) Plackett-Luce as four best/worst trials - the behaviourally faithful
# version. We rebuild each pass as a conditional-logit choice over the items
# that survived to it, two "most" on the forward scale and two "least" on the
# reversed scale, and flip the sign of the item indicators on the "least" tasks
# so one coefficient vector governs both directions (the best-worst convention).
fit_pl_trials <- function(qsort) {
  n_items <- max(qsort$item)

  # One choice task per (respondent, pass): the items still in play, which one
  # was chosen, and the scale direction (+1 for "most", -1 for "least"). Each
  # pass is just a filter on the surviving bands, stacked.
  tasks <- bind_rows(
    qsort %>% filter(band %in% c(5,4,3,2,1)) %>% mutate(task = "most",       chosen = as.integer(band == 5), direction =  1),
    qsort %>% filter(band %in% c(4,3,2,1))   %>% mutate(task = "next_most",  chosen = as.integer(band == 4), direction =  1),
    qsort %>% filter(band %in% c(3,2,1))     %>% mutate(task = "least",      chosen = as.integer(band == 1), direction = -1),
    qsort %>% filter(band %in% c(3,2))       %>% mutate(task = "next_least", chosen = as.integer(band == 2), direction = -1)
  ) %>%
    mutate(stratum = paste(respondent, task, sep = "/")) %>%
    select(stratum, item, direction, chosen)

  # Signed item indicators: each row carries +/-1 in its own item column, 0
  # elsewhere; item 1 is the dropped reference.
  design <- tasks %>%
    mutate(row = row_number(), item_col = paste0("item", item)) %>%
    pivot_wider(id_cols = c(row, stratum, chosen),
                names_from = item_col, values_from = direction, values_fill = 0)

  item_cols <- paste0("item", 2:n_items)
  formula <- as.formula(
    paste("chosen ~", paste(c(item_cols, "strata(stratum)"), collapse = " + "))
  )

  # method = "exact": the two-pick passes are multiple-event strata whose exact
  # conditional likelihood is the unordered "choose 2" probability.
  fit <- survival::clogit(formula, data = design, method = "exact")

  # Item 1 is the reference (coefficient 0, no variance). Pad its row and column
  # into an otherwise-free covariance so one softmax delta handles all items.
  beta <- c(0, coef(fit)[item_cols])
  vcov <- matrix(0, n_items, n_items)
  vcov[2:n_items, 2:n_items] <- stats::vcov(fit)[item_cols, item_cols]

  softmax_with_se(beta, vcov) %>%
    mutate(item = seq_len(n_items), .before = 1)
}

# -----------------------------------------------------------------------------
# 4. Scoring estimates against the planted truth
# -----------------------------------------------------------------------------

# Label the true top and bottom items "extreme" and the rest "middle", so we can
# report recovery separately for the items that decisions hinge on.
label_strata <- function(truth, prop_extreme = 0.2) {
  n <- nrow(truth)
  k <- max(1, round(n * prop_extreme))
  truth %>%
    arrange(desc(true_share)) %>%
    mutate(rank = row_number(),
           stratum = if_else(rank <= k | rank > n - k, "extreme", "middle")) %>%
    select(item, true_share, stratum)
}

# Attach the truth and the per-item error to a set of estimates. Coverage is
# left to the summary step (estimate +/- z * se at whatever nominal level), so
# this stays a plain join + subtraction.
score_estimates <- function(estimates, truth, prop_extreme = 0.2) {
  estimates %>%
    left_join(label_strata(truth, prop_extreme), by = "item") %>%
    mutate(error = estimate - true_share)
}

# -----------------------------------------------------------------------------
# 5. Summaries shared by the run script and the report
# -----------------------------------------------------------------------------

# The two-sided z-multiplier for a nominal confidence level.
z_for <- function(level) qnorm(1 - (1 - level) / 2)

# Spacing recovery: bias, RMSE and mean 95% interval width per estimator and
# stratum, averaged over every scored item and replication handed in.
summarise_recovery <- function(scored) {
  scored %>%
    group_by(estimator, stratum) %>%
    summarise(
      bias    = mean(error),
      rmse    = sqrt(mean(error^2)),
      width95 = mean(2 * z_for(0.95) * se),
      .groups = "drop"
    )
}

# Interval honesty: the fraction of estimate +/- z*se intervals that actually
# cover the truth, at each nominal level. An honest method tracks the diagonal
# (coverage ~ nominal); an overconfident one falls below it.
summarise_coverage <- function(scored, levels = c(0.50, 0.80, 0.90, 0.95)) {
  map(levels, function(L) {
    scored %>%
      group_by(estimator, stratum) %>%
      summarise(
        nominal  = L,
        coverage = mean(abs(error) <= z_for(L) * se),
        .groups  = "drop"
      )
  }) %>%
    list_rbind()
}
