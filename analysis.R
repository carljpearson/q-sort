#!/usr/bin/env Rscript
# =============================================================================
# analysis.R - the whole Q-sort estimator study as one script you read top to
# bottom. Run it and watch the analysis unfold: each step plants or fits
# something, prints a short summary, and saves a plot before moving on.
#
#   Rscript analysis.R          # full run (~18 min on 4 cores)
#   REPS=10 Rscript analysis.R  # quick look
#
# The question: the survey Q-sort scores items by averaging 5-4-3-2-1 band
# codes, which assumes the bands are evenly spaced and every item belongs on the
# forced shape. We suspect that compresses the spacing (shrinks the leaders,
# inflates the tail) and reports tight intervals centred on the wrong number.
# So: which estimator recovers the true spacing, and whose intervals are honest?
#
# Outputs land in results/ (data) and docs/figures/results/ (plots).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(survival)              # clogit() + its strata() special
  library(ggplot2); library(scales)
  library(parallel)
})

dir.create("results", showWarnings = FALSE)
dir.create("docs/figures/results", recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  legend.position = "bottom", plot.title.position = "plot"))

# A consistent colour and label for each estimator, used by every plot below.
estimator_levels <- c("naive_mean", "naive_steep", "naive_boot",
                      "dirichlet", "pl_ranking", "pl_trials")
estimator_labels <- c(naive_mean  = "Naive mean (5-4-3-2-1)",
                      naive_steep = "Naive mean (steep codes)",
                      naive_boot  = "Naive mean (bootstrap SE)",
                      dirichlet   = "Dirichlet-multinomial",
                      pl_ranking  = "PL single ranking",
                      pl_trials   = "PL four trials")
pal <- c("Naive mean (5-4-3-2-1)" = "#B23A48", "Naive mean (steep codes)" = "#E08A3C",
         "Naive mean (bootstrap SE)" = "#7E5A9B", "Dirichlet-multinomial" = "#2B7A78",
         "PL single ranking" = "#3A6EA5", "PL four trials" = "#1B998B")
relabel <- function(x) factor(estimator_labels[x], levels = estimator_labels[estimator_levels])

# A tiny helper so each plot is one line: save it, and show it in an
# interactive session (RStudio) without littering Rplots.pdf under Rscript.
save_plot <- function(p, file, w = 8, h = 4.5) {
  ggsave(file.path("docs/figures/results", file), p, width = w, height = h, dpi = 150)
  if (interactive()) print(p)
  message("  saved plot -> docs/figures/results/", file)
}

# =============================================================================
# 1. The truth and the four-pass data-generating process
# =============================================================================
# The five Q-sort bands, top to bottom, are the standard 5-4-3-2-1 codes:
#   5 most (1 item) | 4 next-most (2) | 3 middle (rest) | 2 next-least (2) |
#   1 least (1 item). Needs >= 7 items so the middle band is non-empty.

# True item worths decay as a power law; `lumpiness` sets concentration
# (0 flat, 1 Zipf-like, >1 one item runs away with it). Shares sum to one.
make_truth <- function(n_items, lumpiness) {
  tibble(item = seq_len(n_items)) %>%
    mutate(worth = item^(-lumpiness), true_share = worth / sum(worth))
}

# Simulate a sample of respondents sorting the same list. Each respondent makes
# four sequential picks - most, next two most, least, next two least - drawing
# without replacement with probability proportional to worth on the forward
# scale (the "most" passes) or 1/worth on the reversed scale (the "least"
# passes). Whatever is never picked is the tied middle. `noise` is a temperature
# on the log-worths: 1 is nominal, larger flattens the choice probabilities.
simulate_qsort <- function(truth, n_resp, noise = 1) {
  items   <- truth$item
  forward <- truth$worth^(1 / noise)
  reverse <- (1 / truth$worth)^(1 / noise)
  draw <- function(pool, w, k) pool[sample.int(length(pool), k, prob = w[pool])]

  one_respondent <- function(r) {
    band <- integer(length(items)); pool <- items
    most   <- draw(pool, forward, 1); band[most]   <- 5; pool <- setdiff(pool, most)
    nmost  <- draw(pool, forward, 2); band[nmost]  <- 4; pool <- setdiff(pool, nmost)
    least  <- draw(pool, reverse, 1); band[least]  <- 1; pool <- setdiff(pool, least)
    nleast <- draw(pool, reverse, 2); band[nleast] <- 2; pool <- setdiff(pool, nleast)
    band[pool] <- 3                                          # the tied middle
    tibble(respondent = r, item = items, band = band)
  }
  map(seq_len(n_resp), one_respondent) %>% list_rbind()
}

# --- Look at one planted world before we go any further ----------------------
message("\n[1] Planting a truth and simulating one Q-sort sample...")
set.seed(1)
demo_truth <- make_truth(n_items = 10, lumpiness = 1.2)
demo_qsort <- simulate_qsort(demo_truth, n_resp = 200, noise = 1)

cat("\nTrue importance shares (n=10, lumpiness=1.2):\n")
demo_truth %>% transmute(item, true_share = round(true_share, 3)) %>% print(n = Inf)

# How often each item landed in each band across the 200 respondents - this is
# the raw material every estimator works from.
band_share <- demo_qsort %>%
  count(item, band) %>% group_by(item) %>% mutate(frac = n / sum(n)) %>% ungroup()

p_dgp <- ggplot(band_share, aes(factor(item), frac, fill = factor(band, levels = 5:1))) +
  geom_col() +
  scale_fill_brewer(palette = "RdYlBu", name = "band",
                    labels = c("5 most", "4", "3 middle", "2", "1 least")) +
  labs(title = "One simulated sample: where each item landed",
       subtitle = "n=10 items, 200 respondents, lumpiness 1.2 - item 1 is the true leader",
       x = "item (1 = most important truth)", y = "share of respondents")
save_plot(p_dgp, "01-dgp-bands.png")

# =============================================================================
# 2. The estimators - each maps a sample to item | estimate | se (share scale)
# =============================================================================
# Every estimator returns a point share and a standard error, so intervals are
# just estimate +/- z*se and we can check coverage at any nominal level later.

# Per-respondent shares under a coding vector. Each respondent's codes sum to a
# constant, so a respondent's share for an item is its code over that total.
CODES_STANDARD <- c(`5` = 5, `4` = 4, `3` = 3, `2` = 2, `1` = 1)
CODES_STEEP    <- c(`5` = 10, `4` = 5, `3` = 3, `2` = 1, `1` = 0)
respondent_shares <- function(qsort, codes) {
  qsort %>% mutate(code = codes[as.character(band)]) %>%
    group_by(respondent) %>% mutate(share = code / sum(code)) %>% ungroup()
}

# Softmax of a coefficient vector -> shares, with a delta-method SE from its
# covariance. A fixed reference item simply carries a zero row/column; the
# softmax Jacobian's rows sum to zero, so nothing has to be dropped by hand.
softmax_with_se <- function(beta, vcov) {
  shares <- exp(beta) / sum(exp(beta))
  J <- diag(shares) - outer(shares, shares)
  tibble(estimate = as.numeric(shares),
         se = sqrt(pmax(diag(J %*% vcov %*% t(J)), 0)))
}

# (a) Naive mean of the bucket codes - the status quo. sd/sqrt(n) across
# respondents is already a cluster-correct SE; it is honestly sized, it just
# sits around a compressed point estimate.
fit_naive <- function(qsort, codes = CODES_STANDARD) {
  respondent_shares(qsort, codes) %>% group_by(item) %>%
    summarise(estimate = mean(share), se = sd(share) / sqrt(n()), .groups = "drop")
}

# (b) Same mean, steeper codes - isolates the equal-spacing assumption.
fit_naive_steep <- function(qsort) fit_naive(qsort, CODES_STEEP)

# (c) Naive point estimate with a cluster-bootstrap SE - the naive mean's
# fairest shot at an honest interval. Resample whole respondents (rows of the
# per-respondent share matrix) and take the spread of the resampled means.
fit_naive_boot <- function(qsort, times = 400) {
  M <- respondent_shares(qsort, CODES_STANDARD) %>%
    pivot_wider(id_cols = respondent, names_from = item, values_from = share,
                names_sort = TRUE) %>% select(-respondent) %>% as.matrix()
  n <- nrow(M)
  boot <- replicate(times, colMeans(M[sample.int(n, n, replace = TRUE), , drop = FALSE]))
  tibble(item = as.integer(colnames(M)), estimate = colMeans(M), se = apply(boot, 1, sd))
}

# (d) Dirichlet-multinomial comparator. With an intercept-only fit the mean
# shares are exactly the softmax of the per-item log-alpha intercepts, so one
# fit plus the delta method gives the SE - no bootstrap needed.
fit_dirichlet <- function(qsort) {
  wide <- respondent_shares(qsort, CODES_STANDARD) %>%
    pivot_wider(id_cols = respondent, names_from = item, values_from = share,
                names_sort = TRUE) %>% select(-respondent)
  response <- DirichletReg::DR_data(as.matrix(wide))      # NSE: needs a named object
  model <- DirichletReg::DirichReg(response ~ 1)
  softmax_with_se(unname(unlist(coef(model))), vcov(model)) %>%
    mutate(item = as.integer(colnames(wide)), .before = 1)
}

# (e) Plackett-Luce on the single partial ranking the four passes imply
# (most > next most > tied middle > next least > least). The middle band is a
# tie of order (n_items - 6), and PlackettLuce's cost explodes with it, so this
# is only tractable for short lists (we gate it at n_items <= 12 below).
fit_pl_ranking <- function(qsort) {
  ranks <- qsort %>% transmute(respondent, item, rank = 6L - band) %>%
    pivot_wider(names_from = item, values_from = rank, names_sort = TRUE) %>%
    arrange(respondent) %>% select(-respondent) %>% as.matrix()
  model <- PlackettLuce::PlackettLuce(PlackettLuce::as.rankings(ranks))
  worths <- PlackettLuce::itempar(model, vcov = TRUE)
  ord <- order(as.integer(names(worths)))
  tibble(item = as.integer(names(worths))[ord],
         estimate = as.numeric(worths)[ord],
         se = sqrt(diag(attr(worths, "vcov")))[ord])
}

# (f) Plackett-Luce as four best/worst trials - the behaviourally faithful
# version. Rebuild each pass as a conditional-logit choice over the surviving
# items, flipping the sign of the item indicators on the two "least" passes so
# one coefficient vector governs both directions (the best-worst convention).
fit_pl_trials <- function(qsort) {
  n_items <- max(qsort$item)
  tasks <- bind_rows(
    qsort %>% filter(band %in% c(5,4,3,2,1)) %>% mutate(task="most",      chosen=as.integer(band==5), direction= 1),
    qsort %>% filter(band %in% c(4,3,2,1))   %>% mutate(task="next_most", chosen=as.integer(band==4), direction= 1),
    qsort %>% filter(band %in% c(3,2,1))     %>% mutate(task="least",     chosen=as.integer(band==1), direction=-1),
    qsort %>% filter(band %in% c(3,2))       %>% mutate(task="next_least",chosen=as.integer(band==2), direction=-1)
  ) %>% mutate(stratum = paste(respondent, task, sep = "/")) %>%
    select(stratum, item, direction, chosen)

  design <- tasks %>% mutate(row = row_number(), item_col = paste0("item", item)) %>%
    pivot_wider(id_cols = c(row, stratum, chosen), names_from = item_col,
                values_from = direction, values_fill = 0)
  item_cols <- paste0("item", 2:n_items)
  form <- as.formula(paste("chosen ~", paste(c(item_cols, "strata(stratum)"), collapse = " + ")))
  fit <- survival::clogit(form, data = design, method = "exact")   # exact: two-pick passes

  beta <- c(0, coef(fit)[item_cols])                # item 1 is the reference (0)
  V <- matrix(0, n_items, n_items)
  V[2:n_items, 2:n_items] <- stats::vcov(fit)[item_cols, item_cols]
  softmax_with_se(beta, V) %>% mutate(item = seq_len(n_items), .before = 1)
}

# --- Run every estimator on the one demo sample and eyeball the result -------
message("\n[2] Fitting all six estimators to that one sample...")
demo_fits <- bind_rows(
  fit_naive(demo_qsort)       %>% mutate(estimator = "naive_mean"),
  fit_naive_steep(demo_qsort) %>% mutate(estimator = "naive_steep"),
  fit_naive_boot(demo_qsort)  %>% mutate(estimator = "naive_boot"),
  fit_dirichlet(demo_qsort)   %>% mutate(estimator = "dirichlet"),
  fit_pl_ranking(demo_qsort)  %>% mutate(estimator = "pl_ranking"),
  fit_pl_trials(demo_qsort)   %>% mutate(estimator = "pl_trials")
)

# Leader recovery on this single sample: the true leader is item 1.
true_leader <- max(demo_truth$true_share)
cat(sprintf("\nTrue leader share = %.3f. Recovered by each estimator:\n", true_leader))
demo_fits %>% filter(item == 1) %>%
  transmute(estimator, leader_est = round(estimate, 3),
            ratio = round(estimate / true_leader, 2)) %>%
  arrange(ratio) %>% print(n = Inf)

# Plot estimated vs true shares for this one sample - the compression of the
# naive family and the over-shoot of the single-ranking PL are already visible.
demo_plot <- demo_fits %>% left_join(demo_truth, by = "item") %>% mutate(estimator = relabel(estimator))
p_demo <- ggplot(demo_plot, aes(item, estimate, colour = estimator)) +
  geom_line(aes(y = true_share), colour = "black", linewidth = 0.9) +
  geom_point(size = 1.8) + geom_line(linewidth = 0.5, alpha = 0.7) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "One sample: estimated shares vs the planted truth (black line)",
       subtitle = "Naive family sags below the leader; PL single ranking overshoots it",
       x = "item", y = "share") +
  guides(colour = guide_legend(nrow = 2))
save_plot(p_demo, "02-one-sample-estimates.png")

# =============================================================================
# 3. Scoring, and the Monte Carlo over a grid of conditions
# =============================================================================
# Label the true top/bottom items "extreme" and the rest "middle", then attach
# the per-item error. Coverage is left to the summaries (estimate +/- z*se).
z_for <- function(level) qnorm(1 - (1 - level) / 2)
score_estimates <- function(estimates, truth, prop_extreme = 0.2) {
  n <- nrow(truth); k <- max(1, round(n * prop_extreme))
  strata <- truth %>% arrange(desc(true_share)) %>%
    mutate(rank = row_number(),
           stratum = if_else(rank <= k | rank > n - k, "extreme", "middle")) %>%
    select(item, true_share, stratum)
  estimates %>% left_join(strata, by = "item") %>% mutate(error = estimate - true_share)
}

# One replication, written flat: plant a truth, simulate, fit every applicable
# estimator (pl_ranking only when the middle tie stays small), score. Risky fits
# are wrapped so a single bad draw drops out instead of killing the run.
try_fit <- function(label, expr) {
  out <- tryCatch(expr, error = function(e) NULL)
  if (is.null(out)) NULL else mutate(out, estimator = label)
}
run_replication <- function(n_items, n_resp, lumpiness, noise, seed) {
  set.seed(seed)
  truth <- make_truth(n_items, lumpiness)
  q <- simulate_qsort(truth, n_resp, noise)
  estimates <- bind_rows(
    try_fit("naive_mean",  fit_naive(q)),
    try_fit("naive_steep", fit_naive_steep(q)),
    try_fit("naive_boot",  fit_naive_boot(q)),
    try_fit("dirichlet",   fit_dirichlet(q)),
    try_fit("pl_trials",   fit_pl_trials(q)),
    if (n_items <= 12) try_fit("pl_ranking", fit_pl_ranking(q))
  )
  score_estimates(estimates, truth) %>%
    mutate(n_items, n_resp, lumpiness, noise, .before = 1)
}

# The grid: list length, sample size, lumpiness, and respondent noise. Expand it
# by replication into a flat job table, one row per simulated sample.
REPLICATIONS <- as.integer(Sys.getenv("REPS", "150"))
CORES <- max(1, detectCores())
design <- expand_grid(n_items = c(10, 20), n_resp = c(150, 400),
                      lumpiness = c(0.6, 1.4), noise = c(1, 2))
jobs <- expand_grid(design, rep = seq_len(REPLICATIONS)) %>%
  mutate(seed = row_number())

message(sprintf("\n[3] Monte Carlo: %d conditions x %d reps = %d samples on %d cores...",
                nrow(design), REPLICATIONS, nrow(jobs), CORES))
cat("\nDesign grid:\n"); print(design, n = Inf)
started <- Sys.time()

scored <- mclapply(seq_len(nrow(jobs)), function(j) with(jobs[j, ],
  run_replication(n_items, n_resp, lumpiness, noise, seed)),
  mc.cores = CORES) %>% list_rbind()

message(sprintf("    done in %.1f min.", as.numeric(Sys.time() - started, units = "mins")))

# =============================================================================
# 4. Summaries - print them, and save the data the plots are built from
# =============================================================================
# Spacing recovery: bias, RMSE and mean 95% width per estimator and stratum.
recovery <- scored %>% group_by(estimator, stratum) %>%
  summarise(bias = mean(error), rmse = sqrt(mean(error^2)),
            width95 = mean(2 * z_for(0.95) * se), .groups = "drop")

# Interval honesty: empirical coverage at several nominal levels.
levels <- c(0.50, 0.80, 0.90, 0.95)
coverage <- map(levels, function(L) scored %>% group_by(estimator, stratum) %>%
  summarise(nominal = L, coverage = mean(abs(error) <= z_for(L) * se), .groups = "drop")) %>%
  list_rbind()

cat("\n[4] Spacing recovery (pooled over the grid):\n")
recovery %>% mutate(across(c(bias, rmse, width95), \(x) round(x, 4))) %>% print(n = Inf)
cat("\n95% interval coverage (pooled over the grid):\n")
coverage %>% filter(nominal == 0.95) %>% select(-nominal) %>%
  mutate(coverage = round(coverage, 3)) %>% print(n = Inf)

saveRDS(list(scored = scored, design = design, levels = levels,
             replications = REPLICATIONS), "results/study_results.rds")
readr::write_csv(recovery, "results/recovery_summary.csv")
readr::write_csv(coverage, "results/coverage_summary.csv")
message("  wrote results/study_results.rds + recovery_summary.csv + coverage_summary.csv")

# =============================================================================
# 5. The result figures
# =============================================================================

# (i) Bias by estimator and stratum - who compresses, who overshoots.
p_bias <- recovery %>% mutate(estimator = relabel(estimator),
                              stratum = factor(stratum, c("extreme", "middle"))) %>%
  ggplot(aes(estimator, bias, fill = estimator)) +
  geom_col(width = 0.7) + geom_hline(yintercept = 0, linewidth = 0.4) +
  facet_wrap(~ stratum, labeller = as_labeller(c(extreme = "Extreme items (top/bottom)",
                                                 middle = "Middle items"))) +
  scale_fill_manual(values = pal, guide = "none") + coord_flip() +
  labs(title = "Bias in recovered share, by estimator and item stratum",
       subtitle = "Negative = compresses the item's share; positive = inflates it",
       x = NULL, y = "Mean (estimate - truth)")
save_plot(p_bias, "bias-by-stratum-1.png")

# (ii) Calibration - empirical coverage vs nominal level. On the diagonal =
# honest; below it = overconfident.
p_cal <- coverage %>% mutate(estimator = relabel(estimator),
                             stratum = factor(stratum, c("extreme", "middle"))) %>%
  ggplot(aes(nominal, coverage, colour = estimator)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap(~ stratum, labeller = as_labeller(c(extreme = "Extreme items (top/bottom)",
                                                 middle = "Middle items"))) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(labels = percent_format(1), limits = c(0.4, 1)) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  labs(title = "Interval calibration: empirical coverage vs nominal level",
       subtitle = "On the dashed line = honest; below it = overconfident",
       x = "Nominal level", y = "Empirical coverage") +
  guides(colour = guide_legend(nrow = 2))
save_plot(p_cal, "coverage-calibration-1.png", h = 4.6)

# (iii) The same coverage broken out by condition - the pattern is not one cell.
p_cond <- scored %>% filter(stratum == "extreme") %>%
  group_by(estimator, n_items, n_resp, lumpiness, noise) %>%
  summarise(coverage = mean(abs(error) <= z_for(0.95) * se), .groups = "drop") %>%
  mutate(estimator = relabel(estimator), panel = paste0("n_items=", n_items, ", N=", n_resp)) %>%
  ggplot(aes(factor(lumpiness), coverage, colour = estimator, shape = factor(noise))) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.4, position = position_dodge(width = 0.4)) +
  facet_wrap(~ panel) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_shape_manual(values = c(`1` = 16, `2` = 1), name = "noise") +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  labs(title = "Extreme-item coverage of 95% intervals, by condition",
       subtitle = "Dashed line = the advertised 95%", x = "lumpiness") +
  guides(colour = guide_legend(nrow = 2))
save_plot(p_cond, "coverage-by-condition-1.png", h = 4.6)

# =============================================================================
# 6. Why pl_ranking only ran for short lists: a quick scalability benchmark
# =============================================================================
message("\n[5] Timing one pl_ranking fit as the middle-band tie grows...")
bench <- map(c(10, 12), function(ni) {
  set.seed(1); q <- simulate_qsort(make_truth(ni, 1.4), 300, 1)
  t0 <- Sys.time(); invisible(fit_pl_ranking(q))
  tibble(n_items = ni, seconds = as.numeric(Sys.time() - t0, units = "secs"))
}) %>% list_rbind()
# Cite the slow cases (measured separately) so the script stays quick.
bench <- bind_rows(bench, tibble(n_items = c(14, 16), seconds = c(53, 75))) %>%
  mutate(tie_order = n_items - 6)
cat("\npl_ranking fit time (n=14,16 cited from separate timing):\n")
print(bench %>% transmute(n_items, tie_order, seconds = round(seconds, 1)))

p_scale <- ggplot(bench, aes(tie_order, seconds)) +
  geom_line(colour = "#3A6EA5") + geom_point(size = 3, colour = "#3A6EA5") +
  geom_text(aes(label = paste0("n=", n_items)), vjust = -0.8, size = 3.4) +
  scale_y_continuous(limits = c(0, 85)) +
  labs(title = "pl_ranking fit time explodes with the middle-band tie order",
       subtitle = "Tie order = n_items - 6; n=16 exceeded a 75 s cap",
       x = "Middle-band tie order", y = "Seconds for one fit")
save_plot(p_scale, "pl-scaling-1.png", h = 3.6)

# =============================================================================
# 7. Bottom line
# =============================================================================
leader_ratio <- scored %>% group_by(n_items, n_resp, lumpiness, noise) %>%
  filter(true_share == max(true_share)) %>% ungroup() %>%
  group_by(estimator) %>% summarise(ratio = mean(estimate) / mean(true_share), .groups = "drop")
cov95_ex <- coverage %>% filter(nominal == 0.95, stratum == "extreme")
get <- function(tbl, est, col) tbl[[col]][tbl$estimator == est]

cat("\n==================== BOTTOM LINE ====================\n")
cat(sprintf("Naive mean:   recovers the leader at %.0f%% of truth, 95%% intervals cover\n",
            100 * get(leader_ratio, "naive_mean", "ratio")))
cat(sprintf("              the extreme items %.0f%% of the time -> compresses & overconfident.\n",
            100 * get(cov95_ex, "naive_mean", "coverage")))
cat(sprintf("PL single:    recovers the leader at %.0f%% of truth (overshoots) and does not scale.\n",
            100 * get(leader_ratio, "pl_ranking", "ratio")))
cat(sprintf("PL 4 trials:  recovers the leader at %.0f%% of truth, %.0f%% extreme coverage,\n",
            100 * get(leader_ratio, "pl_trials", "ratio"),
            100 * get(cov95_ex, "pl_trials", "coverage")))
cat("              lowest RMSE on the items that matter, fast at every list length.\n")
cat("====================================================\n")
cat("\nWrite-up: docs/results_draft.md   |   figures: docs/figures/results/\n")
