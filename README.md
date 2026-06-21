# Survey Q-Sort Estimator Recovery Study

_co-authored with Claude Opus 4.8_

A Monte Carlo parameter-recovery study testing **which estimators recover the true importance of items from survey Q-sort data — and whether their confidence intervals are honest.**

> Honeycrisp Research · Quantitative UX Methodology
> Status: **in development** (scaffolding for implementation in R)

---

## Background

The **survey Q-sort** (the four-pass "most / least" exercise popularized for prioritization) is a fast, software-light way to rank a list of features, needs, or messages. Its standard scoring — average the `5-4-3-2-1` band codes — quietly assumes the bands are **evenly spaced** and that **every item belongs** on the forced shape. Real importance is lumpy and the forced shape caps how far a leader can pull away, so the naive average **compresses the spacing**: it shrinks the leaders and inflates the tail. The ranking usually survives; the spacing does not.

This repo runs a simulation study to quantify that — the **bias** of competing estimators and, crucially, the **confidence-interval coverage** (does a stated 95% interval actually contain the truth 95% of the time?).

**What we are estimating:** population-level importance as **shares that sum to 1**, where *spacing* matters (how much more important one item is than another), not just rank order.

## The question

> Across realistic conditions, which estimator recovers the true **spacing** accurately, and whose intervals are actually **honest** (cover the truth at their stated rate)?

## Approach

Invent a world where the true item worths are **known** → simulate respondents going through the four-pass exercise → fit each estimator → map every estimate to shares on a common scale → compare to the planted truth → repeat across a grid of conditions.

See [`docs/methodology.md`](docs/methodology.md) for the full framing (context, the problem, the model table, and an appendix on the likelihood split that motivates modeling the four passes as separate trials).

## Estimators compared

| Estimator | Role |
|---|---|
| Naive mean of `5-4-3-2-1` codes | Baseline (status quo) |
| Naive mean, alternative coding vectors | Baseline stress-test (isolates the equal-spacing assumption) |
| Plackett-Luce — single partial ranking | Proposed (convenient default) |
| Plackett-Luce — four-trial factorization | Proposed (behaviorally faithful; "least" passes on the reversed scale) |
| Respondent-level (cluster) bootstrap for the naive mean | Required check |
| Dirichlet-multinomial / compositional | Optional comparator |
| Hierarchical Bayes | Later / only if heterogeneity matters |

Full descriptions and rationale in [`docs/methodology.md`](docs/methodology.md).

## Evaluation

- **Spacing accuracy** — bias and RMSE of estimated vs true shares, including the size of the gaps between items (not just whether the order is right).
- **Interval honesty** — empirical vs nominal coverage at several levels (50/80/90/95), reported alongside interval width.
- **Stratified** by item position (extremes vs the muddy middle), and repeated across list length, sample size, how lumpy the truth is, and respondent noise.

---

## Tech stack

R. Key packages by role:

| Package | Role |
|---|---|
| [`SimDesign`](https://cran.r-project.org/package=SimDesign) | Monte Carlo harness (`generate`/`analyse`/`summarise`); built-in `ECR()` coverage, `bias()`, `RMSE()`, `RE()`; auto re-sim of non-convergent fits; parallel + HPC |
| [`PlackettLuce`](https://cran.r-project.org/package=PlackettLuce) | Single-ranking PL — ties of any order, partial rankings, quasi-SEs, `pltree()` for segmentation, optional MVN prior |
| [`survival`](https://cran.r-project.org/package=survival) (`clogit`) | Exploded / stratified conditional logit — engine for the four-trial best-worst model |
| [`support.BWS`](https://cran.r-project.org/package=support.BWS) | Best-worst data encoding (`bws.dataset` → maxdiff / marginal / sequential forms) and the clogit pattern |
| [`apollo`](https://cran.r-project.org/package=apollo) / [`mlogit`](https://cran.r-project.org/package=mlogit) / [`logitr`](https://cran.r-project.org/package=logitr) | Alternative choice-model engines (explicit best+worst likelihood; random coefficients later) |
| [`rsample`](https://cran.r-project.org/package=rsample) (`group_bootstraps`) / `boot` | Respondent-level (cluster) bootstrap |
| [`DirichletReg`](https://cran.r-project.org/package=DirichletReg) / `brms` | Compositional comparator |
| [`bwsTools`](https://cran.r-project.org/package=bwsTools) | Quick BWS scoring sanity checks (`ae_mnl`, empirical Bayes, etc.) |
| *(optional, later)* `apollo` / `RSGHB` / Stan + [`SBC`](https://cran.r-project.org/package=SBC) | Hierarchical Bayes + Bayesian calibration |

## Off-the-shelf vs custom

The estimation engines and the entire Monte-Carlo / coverage machinery already exist and are mature. **The custom layer is small and Q-sort-specific:**

1. **The data-generating process** (`R/dgp.R`) — no package simulates the four-pass full-list response. Given true worths and a choice-noise scale: draw most via softmax over `v`, remove; draw the top-2 set; switch to the reversed scale and draw least via softmax over `1/v` on survivors; draw the bottom-2 set; the rest are the tied middle.
2. **Two encoders** —
   - `R/encode_ranking.R`: map each respondent's four passes to a `PlackettLuce` rankings object with the tie blocks (top-1, tied-2, tied-4, tied-2, bottom-1).
   - `R/encode_trials.R`: explode the four passes into stacked `clogit` tasks, negating the item dummies on the two "least" tasks so a single coefficient vector governs both directions (the standard best-worst sign convention).
3. **Common-scale glue** (`R/evaluate.R`) — normalize every estimator's output to shares and compute coverage separately for extreme vs middle items (lives in SimDesign's `Analyse`/`Summarise`).

Everything else is wiring these together.

## Proposed repository structure

```
qsort-recovery/
├── README.md
├── docs/
│   └── methodology.md          # context, problem, model table, appendix (the likelihood split)
├── figures/
│   ├── coverage_calibration.png / .svg
│   └── likelihood_split.png / .svg
├── R/
│   ├── dgp.R                   # four-pass response simulator (CUSTOM)
│   ├── encode_ranking.R        # PlackettLuce tied-ranking encoder (CUSTOM)
│   ├── encode_trials.R         # reversed-scale clogit explosion (CUSTOM)
│   ├── estimators.R            # naive mean (+ codings), PL single, PL/clogit trials, Dirichlet
│   ├── evaluate.R              # share mapping + position-stratified coverage
│   └── simulation.R            # SimDesign Generate / Analyse / Summarise + design grid
├── analysis/
│   └── run_study.R             # entry point
├── renv.lock
└── LICENSE
```

## Getting started

```r
# R >= 4.3 recommended
install.packages("renv")
renv::init()        # then add the packages above to the lockfile
renv::restore()     # reproducible install

# run the study
source("analysis/run_study.R")
```

## Roadmap

- [ ] `dgp.R` — four-pass response simulator + choice-noise parameter
- [ ] `encode_ranking.R` / `encode_trials.R` — the two encoders
- [ ] `estimators.R` — all estimators returning normalized shares + intervals
- [ ] `evaluate.R` — share mapping, position-stratified coverage
- [ ] `simulation.R` — SimDesign harness wired to `ECR`/`bias`/`RMSE`
- [ ] design grid (K, n, worth separation/shape, noise, homogeneous vs mixture)
- [ ] figures + write-up
- [ ] *(optional)* hierarchical-Bayes layer + SBC calibration

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

TBD — © Honeycrisp Research.
