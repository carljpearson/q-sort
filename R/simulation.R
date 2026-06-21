# -----------------------------------------------------------------------------
# simulation.R - SimDesign harness (generate / analyse / summarise)
#
# SimDesign manages the replications, parallelism and error re-runs. Our three
# steps stay thin:
#   generate  -> plant a truth and simulate a Q-sort sample
#   analyse   -> run every estimator, score it, return one row of statistics
#   summarise -> average those statistics across replications
# -----------------------------------------------------------------------------

qsort_generate <- function(condition, fixed_objects = NULL) {
  truth <- make_truth(condition$n_items, condition$lumpiness)
  qsort <- simulate_qsort(truth, condition$n_resp, condition$noise)
  list(truth = truth, qsort = qsort)
}

qsort_analyse <- function(condition, dat, fixed_objects = NULL) {
  scored <- run_all_estimators(dat$qsort) %>%
    evaluate_estimates(dat$truth)

  # Collapse to one statistic per estimator x stratum x metric, then flatten to
  # the named numeric vector SimDesign collects across replications. Keeping
  # `mse` here (not RMSE) lets summarise average it correctly before the sqrt.
  per_replication <- scored %>%
    group_by(estimator, stratum) %>%
    summarise(
      bias = mean(error),
      mse = mean(error^2),
      coverage = mean(covered),
      width = mean(width),
      .groups = "drop"
    ) %>%
    pivot_longer(c(bias, mse, coverage, width), names_to = "metric") %>%
    unite("key", estimator, stratum, metric)

  set_names(per_replication$value, per_replication$key)
}

qsort_summarise <- function(condition, results, fixed_objects = NULL) {
  averaged <- colMeans(results, na.rm = TRUE)

  # mse averages across replications, then becomes RMSE.
  is_mse <- endsWith(names(averaged), "_mse")
  averaged[is_mse] <- sqrt(averaged[is_mse])
  names(averaged)[is_mse] <- sub("_mse$", "_rmse", names(averaged)[is_mse])

  averaged
}

# The grid of conditions: list length, sample size, how lumpy the truth is, and
# how noisy respondents are.
qsort_design <- function() {
  SimDesign::createDesign(
    n_items   = c(10, 20),
    n_resp    = c(150, 400),
    lumpiness = c(0.5, 1.5),
    noise     = 1
  )
}

run_study <- function(replications = 100, parallel = FALSE, ...) {
  SimDesign::runSimulation(
    design = qsort_design(),
    replications = replications,
    generate = qsort_generate,
    analyse = qsort_analyse,
    summarise = qsort_summarise,
    packages = c(
      "dplyr", "tidyr", "purrr", "tibble",
      "PlackettLuce", "survival", "rsample", "DirichletReg"
    ),
    parallel = parallel,
    ...
  )
}
