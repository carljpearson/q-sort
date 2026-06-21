#!/usr/bin/env Rscript
# =============================================================================
# tests/smoke_test.R - a fast end-to-end check of the whole pipeline
#
# Simulates one Q-sort sample, runs every estimator, and asserts the basics:
# valid shares, sensible standard errors, the right recovery direction, and
# working scoring. Run it before a full study:  Rscript tests/smoke_test.R
# =============================================================================

source("R/qsort.R")

set.seed(1)

ok <- function(label, condition) {
  cat(if (isTRUE(condition)) "PASS " else "FAIL ", label, "\n", sep = "")
  if (!isTRUE(condition)) stop("smoke test failed: ", label)
}

# --- data-generating process -------------------------------------------------
truth <- make_truth(n_items = 10, lumpiness = 1.2)
qsort <- simulate_qsort(truth, n_resp = 200, noise = 1)

ok("truth shares sum to 1", abs(sum(truth$true_share) - 1) < 1e-9)
ok("every respondent x item present", nrow(qsort) == 200 * 10)

band_sizes <- qsort %>% filter(respondent == 1) %>% count(band) %>% arrange(band)
ok("band sizes are 1,2,4,2,1 (least..most)",
   identical(band_sizes$n, c(1L, 2L, 4L, 2L, 1L)))

# --- estimators all share the item | estimate | se contract ------------------
fits <- list(
  naive       = fit_naive(qsort),
  naive_steep = fit_naive_steep(qsort),
  naive_boot  = fit_naive_boot(qsort, times = 300),
  dirichlet   = fit_dirichlet(qsort),
  pl_ranking  = fit_pl_ranking(qsort),
  pl_trials   = fit_pl_trials(qsort)
)

for (nm in names(fits)) {
  est <- fits[[nm]]
  ok(paste0(nm, ": shares sum to ~1"), abs(sum(est$estimate) - 1) < 1e-6)
  ok(paste0(nm, ": standard errors finite & >= 0"), all(est$se >= 0 & is.finite(est$se)))
  ok(paste0(nm, ": one row per item"), nrow(est) == 10)
}

# Recovery sanity: the top true item ranks top under every estimator.
top_item <- truth %>% slice_max(true_share, n = 1) %>% pull(item)
for (nm in names(fits)) {
  ok(paste0(nm, ": recovers the top item"),
     fits[[nm]]$item[which.max(fits[[nm]]$estimate)] == top_item)
}

# Spacing: the four-trial PL should be less compressed than the naive mean.
cat(sprintf("\nleader share - truth %.3f | naive %.3f | pl_trials %.3f | pl_ranking %.3f\n",
            max(truth$true_share), max(fits$naive$estimate),
            max(fits$pl_trials$estimate), max(fits$pl_ranking$estimate)))
ok("pl_trials less compressed than naive at the top",
   max(fits$pl_trials$estimate) >= max(fits$naive$estimate))

# --- scoring and summaries ---------------------------------------------------
scored <- imap(fits, ~ mutate(.x, estimator = .y)) %>%
  list_rbind() %>%
  score_estimates(truth)

ok("scored has both strata", all(c("extreme", "middle") %in% scored$stratum))

cat("\n--- spacing recovery (single sample) ---\n")
print(summarise_recovery(scored), n = Inf)

cat("\n--- 95% coverage (single sample) ---\n")
print(summarise_coverage(scored, levels = 0.95) %>% select(-nominal), n = Inf)

cat("\nAll smoke tests passed.\n")
