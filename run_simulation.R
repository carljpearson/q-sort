#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# run_simulation.R - entry point
#
# Sources the pieces and runs the Monte Carlo study. Lower `replications` for a
# quick check; raise it (and set parallel = TRUE) for a full run.
#
# Packages from CRAN: SimDesign, PlackettLuce, DirichletReg.
# Packages also on apt: dplyr, tidyr, purrr, tibble, survival, rsample.
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(survival)  # clogit + its strata() special

source("R/dgp.R")
source("R/estimators.R")
source("R/encode_ranking.R")
source("R/encode_trials.R")
source("R/evaluate.R")
source("R/simulation.R")

set.seed(1)

results <- run_study(replications = 100, parallel = FALSE)
print(results)
