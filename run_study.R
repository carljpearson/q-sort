#!/usr/bin/env Rscript
# =============================================================================
# run_study.R - the Monte Carlo, written as a script you can read top to bottom
#
# For every cell of a design grid we repeatedly: plant a known truth, simulate a
# Q-sort sample, fit every estimator, and score each estimate against the truth.
# The transparent loop below replaces the old SimDesign generate/analyse/
# summarise callbacks - the same study, without the indirection. Parallelism is
# a plain mclapply over (condition x replication) jobs.
#
#   Rscript run_study.R              # full run (settings below)
#   REPS=5 Rscript run_study.R       # quick check with 5 replications
#
# Output lands in results/: the per-item scored draws (the raw material for the
# report) plus two tidy summary tables.
# =============================================================================

source("R/qsort.R")
library(parallel)

# ---- knobs ------------------------------------------------------------------
# Replications and bootstrap depth can be dialled down for a quick look via the
# REPS / BOOT environment variables; the defaults are the full run.
REPLICATIONS   <- as.integer(Sys.getenv("REPS",  "120"))
BOOT_NAIVE     <- as.integer(Sys.getenv("BOOT",  "400"))
NOMINAL_LEVELS <- c(0.50, 0.80, 0.90, 0.95)
CORES          <- max(1, detectCores())
SEED           <- 1

# ---- the design grid --------------------------------------------------------
# List length, sample size, how lumpy the truth is, and how noisy respondents
# are. pl_ranking is attempted only where its middle-band tie stays tractable.
design <- expand_grid(
  n_items   = c(10, 20),
  n_resp    = c(150, 400),
  lumpiness = c(0.6, 1.4),
  noise     = c(1, 2)
) %>%
  mutate(condition = row_number(), .before = 1)

# ---- one replication --------------------------------------------------------
# Plant a truth, simulate a sample, fit every applicable estimator, score it.
# Each fit is wrapped so a single non-convergent draw drops out instead of
# killing the run; pl_ranking is included only when the list is short enough.
run_one <- function(cond, rep, seed) {
  set.seed(seed)
  truth <- make_truth(cond$n_items, cond$lumpiness)
  qsort <- simulate_qsort(truth, cond$n_resp, cond$noise)

  fitters <- list(
    naive_mean  = function(q) fit_naive(q),
    naive_steep = function(q) fit_naive_steep(q),
    naive_boot  = function(q) fit_naive_boot(q, BOOT_NAIVE),
    dirichlet   = function(q) fit_dirichlet(q),
    pl_trials   = function(q) fit_pl_trials(q)
  )
  if (fit_pl_ranking_feasible(cond$n_items)) {
    fitters$pl_ranking <- function(q) fit_pl_ranking(q)
  }

  estimates <- imap(fitters, function(fit, name) {
    out <- tryCatch(fit(qsort), error = function(e) NULL)
    if (is.null(out)) NULL else mutate(out, estimator = name)
  }) %>%
    list_rbind()

  score_estimates(estimates, truth) %>%
    mutate(condition = cond$condition, rep = rep, .before = 1)
}

# ---- run every (condition x replication) job --------------------------------
jobs <- expand_grid(condition = design$condition, rep = seq_len(REPLICATIONS)) %>%
  mutate(seed = SEED + row_number())

cat(sprintf("Running %d conditions x %d replications = %d jobs on %d cores...\n",
            nrow(design), REPLICATIONS, nrow(jobs), CORES))
started <- Sys.time()

scored <- mclapply(seq_len(nrow(jobs)), function(j) {
  cond <- design[design$condition == jobs$condition[j], ]
  run_one(cond, jobs$rep[j], jobs$seed[j])
}, mc.cores = CORES) %>%
  list_rbind() %>%
  left_join(design, by = "condition")

cat(sprintf("Done in %.1f min.\n", as.numeric(Sys.time() - started, units = "mins")))

# ---- summarise and save -----------------------------------------------------
recovery <- summarise_recovery(scored)
coverage <- summarise_coverage(scored, NOMINAL_LEVELS)

dir.create("results", showWarnings = FALSE)
saveRDS(list(scored = scored, design = design, levels = NOMINAL_LEVELS,
             replications = REPLICATIONS),
        "results/study_results.rds")
readr::write_csv(scored,   "results/scored_draws.csv")
readr::write_csv(recovery, "results/recovery_summary.csv")
readr::write_csv(coverage, "results/coverage_summary.csv")

cat("\n--- spacing recovery (pooled over the grid) ---\n")
print(recovery, n = Inf)
cat("\n--- 95% interval coverage (pooled over the grid) ---\n")
print(filter(coverage, nominal == 0.95) %>% select(-nominal), n = Inf)
cat("\nWrote results/ : study_results.rds, scored_draws.csv, recovery_summary.csv, coverage_summary.csv\n")
