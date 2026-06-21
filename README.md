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

R. The packages and what each one does:

| Package | Role |
|---|---|
| [`SimDesign`](https://cran.r-project.org/package=SimDesign) | Monte Carlo harness with generate, analyse and summarise steps. Built-in `ECR()` coverage, `bias()`, `RMSE()`, `RE()`. Re-runs non-convergent fits, parallelizes, handles HPC jobs. |
| [`PlackettLuce`](https://cran.r-project.org/package=PlackettLuce) | Single-ranking Plackett-Luce. Ties of any order, partial rankings, quasi standard errors, `pltree()` for segmentation, optional multivariate normal prior. |
| [`survival`](https://cran.r-project.org/package=survival) (`clogit`) | Exploded conditional logit, the engine for the four-trial best-worst model. |
| [`support.BWS`](https://cran.r-project.org/package=support.BWS) | Best-worst data encoding through `bws.dataset` (maxdiff, marginal and sequential forms) and the clogit pattern. |
| [`apollo`](https://cran.r-project.org/package=apollo) / [`mlogit`](https://cran.r-project.org/package=mlogit) / [`logitr`](https://cran.r-project.org/package=logitr) | Alternative choice-model engines for an explicit best-and-worst likelihood, with random coefficients available later. |
| [`rsample`](https://cran.r-project.org/package=rsample) (`group_bootstraps`) / `boot` | Respondent-level cluster bootstrap. |
| [`DirichletReg`](https://cran.r-project.org/package=DirichletReg) / `brms` | Compositional comparator. |
| [`bwsTools`](https://cran.r-project.org/package=bwsTools) | Quick best-worst scoring checks like `ae_mnl` and empirical Bayes. |
| `apollo` / `RSGHB` / Stan + [`SBC`](https://cran.r-project.org/package=SBC) | Hierarchical Bayes and Bayesian calibration, if needed later. |

## What I am writing versus what I am reusing

Most of this is assembly. The estimation engines and the Monte Carlo and coverage apparatus already exist and are stable. What I have to write is a small layer specific to the Q-sort.

The data-generating process in `R/dgp.R` comes first, since no package simulates the four-pass full-list response. Given true worths and a noise scale, it draws the most important item by softmax over the worths, removes it, draws the next two, switches to the reversed scale to draw the least important among the survivors, takes the next two down, and leaves the rest as the tied middle.

Two encoders come next. `R/encode_ranking.R` maps each respondent's four passes into a PlackettLuce rankings object with the tie blocks. `R/encode_trials.R` explodes the same passes into stacked clogit tasks and flips the sign of the item dummies on the two least passes, so one coefficient vector governs both directions. That sign flip is the standard best-worst convention.

The glue in `R/evaluate.R` finishes it off, putting every estimator's output onto the same share scale and computing coverage separately for the extreme and middle items. That work sits inside SimDesign's analyse and summarise steps. The rest is wiring.

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
