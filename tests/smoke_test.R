#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# tests/smoke_test.R - exercise the parts that run without CRAN packages
#
# Covers the DGP, the bucket-code estimators, the four-trial Plackett-Luce
# (survival::clogit), the cluster bootstrap (rsample) and the evaluation. The
# PlackettLuce-ranking and DirichletReg estimators, and the SimDesign harness,
# need CRAN and are checked separately once that access is available.
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(survival)  # clogit + its strata() special

source("R/dgp.R")
source("R/estimators.R")
source("R/encode_trials.R")
source("R/evaluate.R")

set.seed(1)

ok <- function(label, condition) {
  cat(if (isTRUE(condition)) "PASS " else "FAIL ", label, "\n", sep = "")
  if (!isTRUE(condition)) stop("smoke test failed: ", label)
}

# --- DGP -------------------------------------------------------------------
truth <- make_truth(n_items = 10, lumpiness = 1.2)
qsort <- simulate_qsort(truth, n_resp = 200, noise = 1)

ok("truth shares sum to 1", abs(sum(truth$true_share) - 1) < 1e-9)
ok("every respondent x item present", nrow(qsort) == 200 * 10)

band_counts <- qsort %>% count(respondent, band) %>% count(band, name = "rows")
expected <- tibble(band = 1:5, rows = c(200, 200, 200, 200, 200))
sizes_per_resp <- qsort %>% filter(respondent == 1) %>% count(band) %>% arrange(band)
ok("band sizes are 1,2,4,2,1 (least..most)",
   identical(sizes_per_resp$n, c(1L, 2L, 4L, 2L, 1L)))

# --- estimators ------------------------------------------------------------
naive <- est_naive_mean(qsort)
alt   <- est_naive_alt(qsort)
boot  <- est_naive_bootstrap(qsort, times = 200)
pl    <- est_pl_trials(qsort)

for (nm in c("naive", "alt", "boot", "pl")) {
  est <- get(nm)
  ok(paste0(nm, ": shares sum to ~1"), abs(sum(est$estimate) - 1) < 1e-6)
  ok(paste0(nm, ": intervals ordered"), all(est$conf_low <= est$conf_high))
  ok(paste0(nm, ": one row per item"), nrow(est) == 10)
}

# Recovery sanity: the top true item should rank top under every estimator.
top_item <- truth %>% slice_max(true_share, n = 1) %>% pull(item)
for (nm in c("naive", "alt", "boot", "pl")) {
  est <- get(nm)
  ok(paste0(nm, ": recovers the top item"),
     est$item[which.max(est$estimate)] == top_item)
}

# Spacing: the four-trial PL should be less compressed than the naive mean,
# i.e. give the true leader a larger share.
lead_true  <- max(truth$true_share)
lead_naive <- max(naive$estimate)
lead_pl    <- max(pl$estimate)
cat(sprintf("\nleader share - truth %.3f | naive %.3f | pl_trials %.3f\n",
            lead_true, lead_naive, lead_pl))
ok("pl_trials less compressed than naive at the top", lead_pl >= lead_naive)

# --- evaluation ------------------------------------------------------------
estimates <- bind_rows(
  naive %>% mutate(estimator = "naive_mean"),
  alt   %>% mutate(estimator = "naive_alt_coding"),
  boot  %>% mutate(estimator = "naive_bootstrap"),
  pl    %>% mutate(estimator = "pl_trials")
)
scored <- evaluate_estimates(estimates, truth)

ok("scored has stratum labels", all(c("extreme", "middle") %in% scored$stratum))
ok("covered is 0/1", all(scored$covered %in% c(0L, 1L)))

summary_tbl <- scored %>%
  group_by(estimator, stratum) %>%
  summarise(
    bias = mean(error),
    rmse = sqrt(mean(error^2)),
    coverage = mean(covered),
    width = mean(width),
    .groups = "drop"
  )

cat("\n--- per-estimator scoring (single sample) ---\n")
print(summary_tbl, n = Inf)

cat("\nAll smoke tests passed.\n")
