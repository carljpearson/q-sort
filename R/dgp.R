# -----------------------------------------------------------------------------
# dgp.R - Data-generating process for the survey Q-sort
#
# Plants a known "truth" (item worths / shares) and simulates respondents going
# through the four-pass Q-sort: pick the most important item, then the next two,
# then switch to the reversed scale to pick the least important and the next two
# down. Whatever is never picked lands in the tied middle.
#
# Bands - and the standard 5-4-3-2-1 codes - line up like this:
#   band 5  most         (1 item,  forward / max-type draw)
#   band 4  next most    (2 items, forward / max-type draw)
#   band 3  middle       (the rest, never picked)
#   band 2  next least   (2 items, reversed / min-type draw)
#   band 1  least        (1 item,  reversed / min-type draw)
#
# Needs at least 7 items so the middle band is non-empty.
# -----------------------------------------------------------------------------

# Pass sizes, in the order respondents actually choose.
qsort_pass_sizes <- c(most = 1, next_most = 2, least = 1, next_least = 2)

# True item worths for a list of `n_items`, decaying as a power law.
# `lumpiness` controls how concentrated importance is:
#   0   -> flat (every item equally important)
#   1   -> Zipf-like (1, 1/2, 1/3, ...)
#   >1  -> one or two items run away with it
make_truth <- function(n_items, lumpiness) {
  tibble(item = seq_len(n_items)) %>%
    mutate(
      worth = item^(-lumpiness),
      true_share = worth / sum(worth)
    )
}

# Draw `size` items from `items` without replacement, with probability
# proportional to `weights`. R's sample() applies the weights sequentially
# (draw, remove, renormalise, draw again) - exactly the Plackett-Luce process.
draw_items <- function(items, weights, size) {
  picked <- sample.int(length(items), size = size, replace = FALSE, prob = weights)
  items[picked]
}

# One respondent's four passes -> a band for every item.
simulate_respondent <- function(worths, noise, sizes = qsort_pass_sizes) {
  items <- seq_along(worths)

  # `noise` is a temperature on the log-worths: 1 is the nominal model, larger
  # flattens the choice probabilities (noisier respondents), smaller sharpens.
  forward <- worths^(1 / noise)        # max-type scale, for the "most" passes
  reverse <- (1 / worths)^(1 / noise)  # min-type scale, for the "least" passes

  band <- integer(length(items))
  remaining <- items

  most <- draw_items(remaining, forward[remaining], sizes["most"])
  band[most] <- 5
  remaining <- setdiff(remaining, most)

  next_most <- draw_items(remaining, forward[remaining], sizes["next_most"])
  band[next_most] <- 4
  remaining <- setdiff(remaining, next_most)

  least <- draw_items(remaining, reverse[remaining], sizes["least"])
  band[least] <- 1
  remaining <- setdiff(remaining, least)

  next_least <- draw_items(remaining, reverse[remaining], sizes["next_least"])
  band[next_least] <- 2
  remaining <- setdiff(remaining, next_least)

  band[remaining] <- 3  # the tied middle

  tibble(item = items, band = band)
}

# A full sample: `n_resp` respondents sorting the same list.
simulate_qsort <- function(truth, n_resp, noise = 1) {
  seq_len(n_resp) %>%
    map(function(r) {
      simulate_respondent(truth$worth, noise) %>%
        mutate(respondent = r, .before = 1)
    }) %>%
    list_rbind()
}
