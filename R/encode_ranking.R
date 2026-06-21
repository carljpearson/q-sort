# -----------------------------------------------------------------------------
# encode_ranking.R - Plackett-Luce on the single partial ranking
#
# Collapses the four passes into the one partial ranking they imply
# (most > next most > middle > next least > least, with ties inside each band)
# and fits one Plackett-Luce likelihood. Worths come back already normalised to
# sum to one (itempar), i.e. on the share scale, with a variance-covariance
# matrix for delta-method intervals.
# -----------------------------------------------------------------------------

# Bands map to ranks: band 5 (most) is rank 1, band 1 (least) is rank 5. Tied
# bands share a rank, which PlackettLuce handles natively.
encode_ranking <- function(qsort) {
  wide <- qsort %>%
    transmute(respondent, item, rank = 6L - band) %>%
    pivot_wider(names_from = item, values_from = rank, names_sort = TRUE) %>%
    arrange(respondent)

  ranking_matrix <- wide %>%
    select(-respondent) %>%
    as.matrix()

  PlackettLuce::as.rankings(ranking_matrix)
}

est_pl_ranking <- function(qsort, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)

  model <- qsort %>%
    encode_ranking() %>%
    PlackettLuce::PlackettLuce()

  worths <- PlackettLuce::itempar(model, vcov = TRUE)  # sum to 1 (share scale)
  ordering <- order(as.integer(names(worths)))

  estimate <- as.numeric(worths)[ordering]
  se <- sqrt(diag(attr(worths, "vcov")))[ordering]

  tibble(
    item = as.integer(names(worths))[ordering],
    estimate = estimate,
    conf_low = estimate - z * se,
    conf_high = estimate + z * se
  )
}
