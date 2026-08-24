context("Tests for the stanr/stanli backend variable filtering")
skip("not run by default")

# These tests require the 'stanr' package (which stanli is embedded in) and
# fit real models, so they are skipped by default like the read_csv_as_stanfit
# tests. Run manually with testthat::test_file() when working on the stanr
# backend.

test_that("stanr backend excludes raw z_1 parameters and keeps lprior/lp__", {
  skip_if_not_installed("stanr")

  set.seed(6231)
  n <- 60
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    g = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )

  fit_stanr <- brm(
    y ~ x + (1 | g), data = dat, backend = "stanr",
    chains = 1, iter = 40, warmup = 20, seed = 8213, refresh = 0
  )

  vars_stanr <- variables(fit_stanr)

  # raw non-centered parameters should never be exposed by default
  # (save_pars$all defaults to FALSE)
  expect_length(grep("^z_1(\\[|$)", vars_stanr), 0)

  # one group-level effect per factor level should be present
  expect_length(grep("^r_g\\[", vars_stanr), nlevels(dat$g))

  # 'lprior' and 'lp__' should be retained like in the rstan/cmdstanr backends
  expect_true("lprior" %in% vars_stanr)
  expect_true("lp__" %in% vars_stanr)
})

test_that("stanr and cmdstanr backends expose the same set of variables", {
  skip_if_not_installed("stanr")
  skip_if_not_installed("cmdstanr")

  set.seed(6231)
  n <- 60
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    g = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )

  fit_common_args <- list(
    formula = y ~ x + (1 | g), data = dat,
    chains = 1, iter = 40, warmup = 20, seed = 8213, refresh = 0
  )

  fit_stanr <- do.call(brm, c(fit_common_args, list(backend = "stanr")))
  fit_cmdstanr <- do.call(brm, c(fit_common_args, list(backend = "cmdstanr")))

  vars_stanr <- sort(variables(fit_stanr))
  vars_cmdstanr <- sort(variables(fit_cmdstanr))

  expect_identical(vars_stanr, vars_cmdstanr)
})
