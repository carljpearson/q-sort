# -----------------------------------------------------------------------------
# encode_trials.R - Plackett-Luce as four separate best/worst trials
#
# Keeps the elicitation's real structure: four conditional choices over the set
# that survived to each pass - two "most" on the forward scale, two "least" on
# the reversed scale. We stack them as conditional-logit choice tasks and flip
# the sign of the item indicators on the two "least" tasks, so one coefficient
# vector governs both directions (the standard best-worst convention).
# -----------------------------------------------------------------------------

# One pass for one respondent: keep the surviving bands, mark the chosen items,
# and record the scale direction (+1 forward / -1 reversed).
trial_task <- function(respondent_bands, task, survivors, chosen, direction) {
  respondent_bands %>%
    filter(band %in% survivors) %>%
    transmute(
      task = task,
      item,
      chosen = as.integer(band %in% chosen),
      direction = direction
    )
}

# All four passes for one respondent.
build_respondent_tasks <- function(respondent_bands) {
  bind_rows(
    trial_task(respondent_bands, "most",       c(5, 4, 3, 2, 1), 5,  1),
    trial_task(respondent_bands, "next_most",  c(4, 3, 2, 1),    4,  1),
    trial_task(respondent_bands, "least",      c(3, 2, 1),       1, -1),
    trial_task(respondent_bands, "next_least", c(3, 2),          2, -1)
  )
}

# Stack every respondent's tasks and give each its own choice stratum.
encode_trials <- function(qsort) {
  qsort %>%
    group_by(respondent) %>%
    group_modify(~ build_respondent_tasks(.x)) %>%
    ungroup() %>%
    mutate(stratum = paste(respondent, task, sep = "/"))
}

est_pl_trials <- function(qsort, level = 0.95) {
  n_items <- max(qsort$item)
  tasks <- encode_trials(qsort)

  # Signed item indicators: each row carries +1 (or -1 on a "least" task) in its
  # own item column, 0 elsewhere. Item 1 is the reference, dropped to identify.
  design <- tasks %>%
    mutate(
      row_id = row_number(),
      item_col = paste0("item", item),
      value = direction
    ) %>%
    pivot_wider(
      id_cols = c(row_id, stratum, chosen),
      names_from = item_col,
      values_from = value,
      values_fill = 0
    )

  item_cols <- paste0("item", 2:n_items)
  formula <- as.formula(
    paste("chosen ~", paste(c(item_cols, "strata(stratum)"), collapse = " + "))
  )

  # method = "exact": the passes that pick two items are multiple-event strata,
  # whose exact conditional likelihood is the unordered "choose 2" probability.
  fit <- survival::clogit(formula, data = design, method = "exact")

  beta <- c(item1 = 0, coef(fit)[item_cols])
  vcov_free <- vcov(fit)[item_cols, item_cols, drop = FALSE]

  softmax_shares_ci(beta, vcov_free, level) %>%
    mutate(item = seq_len(n_items), .before = 1)
}
