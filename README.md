# Survey Q-Sort Estimator Recovery Study

A Monte Carlo study of which estimators recover true item importance from survey Q-sort data, and whether their confidence intervals hold up.

## Background

The survey Q-sort is a quick way to rank a list of features or needs by importance. Respondents work through four passes. They pick the most important item, then the next two, then the least important, then the next two down. Anything they never pick lands in the middle. Each item gets a code from 5 down to 1 based on where it ended up, and you average those codes across everyone to get a ranking.

The standard scoring carries an assumption most people never stop to check. Averaging the 5-4-3-2-1 codes treats the bands as evenly spaced, so the distance from a 5 to a 4 matches the distance from a 2 to a 1. It also assumes every item belongs somewhere on the forced shape, including ones a respondent has no feeling about. Real importance rarely behaves that way. One item usually runs away with it, a few matter, and a long tail barely registers. The forced shape caps how far the leader can pull ahead, since it only gets one slot in the top band. Combine those and the average compresses the spacing. Leaders shrink, the tail inflates, and the order tends to come out right while the gaps between items get flattened.

I run a simulation here to measure that compression. I check the bias of each candidate estimator and, more to the point, whether its confidence intervals are honest, meaning a stated 95% interval really does contain the truth 95% of the time.

The quantity I am after is population-level importance expressed as shares that sum to one, with spacing that carries weight. If one item comes out worth twice another, I want that to mean twice and not merely higher.

## The question

Across realistic conditions, which estimator recovers the true spacing accurately, and whose intervals cover the truth at the rate they advertise? That is what I am trying to settle.

## How it works

I build a world where the true item worths are known, simulate respondents going through the four passes, fit each estimator, map every estimate onto the same share scale, and compare it against the planted truth. Then I repeat across a grid of conditions. The full framing is in [docs/methodology.md](docs/methodology.md), including the appendix on why I treat the four passes as separate trials rather than one ranking.

## Estimators

I put seven estimators through the simulation. The naive bucket-code mean is the status quo I am scrutinizing, joined by a version that swaps in different coding vectors so I can test the even-spacing assumption directly. Plackett-Luce shows up in two forms, one that collapses the four passes into a single ranking and one that models them as separate trials with the least passes on a reversed scale. A respondent-level cluster bootstrap gives the naive mean its fairest shot at honest intervals. A Dirichlet-multinomial model sits between the two camps as a comparator, and I am holding a hierarchical Bayes layer in reserve for the case where the sample splits into segments. The methodology note describes each one and my reasoning.

## Evaluation

Every run measures how close each estimate's shares sit to the truth, gaps included rather than order alone, and whether the intervals are honest, checked at several nominal levels and reported next to their width. I split all of it between the items at the extremes and the ones in the muddy middle, then repeat across list length, sample size, how lumpy the true importance is, and how noisy respondents are.

## Tech stack

R, leaning on the tidyverse. The whole study is small enough to read top to
bottom, so the Monte Carlo is a plain parallel loop rather than a simulation
framework. The packages that actually do the work:

| Package | Role |
|---|---|
| `dplyr` / `tidyr` / `purrr` / `tibble` | The tidy data layer: the data-generating process, encoders, and scoring are all tidy-table transforms. |
| [`survival`](https://cran.r-project.org/package=survival) (`clogit`) | Exploded conditional logit, the engine for the four-trial best-worst model. |
| [`PlackettLuce`](https://cran.r-project.org/package=PlackettLuce) | Single-ranking Plackett-Luce on the tied partial ranking. Tractable only for short lists (the middle band is a high-order tie). |
| [`DirichletReg`](https://cran.r-project.org/package=DirichletReg) | Compositional comparator; the fitted mean shares are a softmax of its log-alpha intercepts, so the interval comes from one fit by the delta method. |
| `parallel` (`mclapply`) | Runs the (condition × replication) jobs across cores. |
| `ggplot2` / `scales` / `knitr` / `rmarkdown` | Figures and the rendered report. |

The respondent-level cluster bootstrap (for the naive mean) is a few lines of
base R — resample the rows of the per-respondent share matrix — so it needs no
extra package. A hierarchical Bayes layer (Stan) and explicit best-worst choice
engines (`apollo`, `mlogit`, `support.BWS`, `bwsTools`) remain options for later
if heterogeneity or a richer likelihood turns out to matter.

## How the code is laid out

The Q-sort-specific layer is deliberately small and readable. Three files:

- **`R/qsort.R`** — the whole library in one place, readable top to bottom:
  `make_truth()` plants a known importance profile; `simulate_qsort()` walks
  respondents through the four passes (draw the most important by softmax over
  the worths, remove it, draw the next two, switch to the reversed scale for the
  least and next-least, leave the rest as the tied middle); the `fit_*()`
  estimators each return the same `item | estimate | se` on the share scale; and
  `score_estimates()` compares any estimate against the planted truth. Two
  encoders live inside their estimators — `fit_pl_ranking()` builds the tied
  partial ranking, and `fit_pl_trials()` stacks the four passes as clogit tasks
  with the sign of the item dummies flipped on the two "least" passes (the
  standard best-worst convention).
- **`run_study.R`** — the Monte Carlo as a transparent script: a design grid, a
  per-replication function, and an `mclapply` over every (condition ×
  replication) job. It writes the scored draws and two summary tables to
  `results/`.
- **`report.Rmd`** — reads those results and renders [`docs/report.md`](docs/report.md)
  with the recovery and coverage figures.

A uniform contract makes this tidy: because every estimator returns a point
share and a standard error, intervals are built downstream as `estimate ± z·se`,
and one scoring step checks coverage at any nominal level.

## Reproducing the study

```sh
Rscript tests/smoke_test.R     # fast end-to-end check on one sample
Rscript run_study.R            # the full Monte Carlo -> results/
Rscript -e 'rmarkdown::render("report.Rmd")'   # -> docs/report.md
```

`run_study.R` reads a few knobs from the environment, so a quick look is
`REPS=5 Rscript run_study.R`. It needs `PlackettLuce`, `DirichletReg`, the
tidyverse core, `survival`, and (for the report) `ggplot2` / `rmarkdown`.

## Results

The rendered findings live in [`docs/report.md`](docs/report.md). In short:
the naive mean **compresses** the spacing (it recovers the leader at a fraction
of its true share and its 95% intervals badly under-cover the top and bottom
items), the single-ranking Plackett-Luce **over-corrects** and does not scale,
and the **four-trial Plackett-Luce** recovers the spacing with the lowest error
and the most honest intervals while fitting quickly at every list length.

## Key references

- Stephenson, W. (1953). *The Study of Behavior: Q-Technique and Its Methodology.* University of Chicago Press.
- Chrzan, K., & Golovashkina, N. (2006). An empirical test of six stated importance measures. *International Journal of Market Research*, 48(6), 717–740.
- Chrzan, K., & Peitz, M. (2019). Best-Worst Scaling with many items. *Journal of Choice Modelling*, 30, 61–72.
- Louviere, J. J., Flynn, T. N., & Marley, A. A. J. (2015). *Best-Worst Scaling: Theory, Methods and Applications.* Cambridge University Press.
- Turner, H. L., van Etten, J., Firth, D., & Kosmidis, I. (2020). Modelling rankings in R: the PlackettLuce package. *Computational Statistics*, 35, 1027–1057.
- Chalmers, R. P., & Adkins, M. C. (2020). Writing effective and reliable Monte Carlo simulations with the SimDesign package. *The Quantitative Methods for Psychology*, 16(4), 248–280.
- Aizaki, H., & Fogarty, J. (2023). R packages and tutorial for case 1 best-worst scaling. *Journal of Choice Modelling*, 46, 100394.
- White, M. H. II (2021). bwsTools: An R package for case 1 best-worst scaling. *Journal of Choice Modelling*, 39, 100289.

## License

TBD. Copyright Honeycrisp Research.
