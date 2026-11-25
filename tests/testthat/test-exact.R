fixed_tree <- function(n, loc, grid) {
  data.frame(node_id = 1:7, parent_id = c(0, 1, 1, 2, 2, 3, 3),
             left = c(2, 4, 6, 0, 0, 0, 0), right = c(3, 5, 7, 0, 0, 0, 0),
             depth = c(0, 1, 1, 2, 2, 2, 2), n = n,
             loc = c(loc, 0, 0, 0, 0), fraction = c(loc / (grid + 1), NA, NA, NA, NA))
}
run_fixed <- function(t, p) {
  ptdensity:::ptd_fixed_cpp(as.matrix(t[, c("left", "right")]), as.integer(t$n), as.integer(t$loc),
                            as.integer(t$depth), p)
}

test_that("fixed-tree messages agree with complete state enumeration", {
  t <- fixed_tree(c(38, 18, 20, 12, 6, 15, 5), c(2, 4, 6), 7)
  priors <- list(density_prior(max_depth = 2, split_grid = 7),
                 density_prior(max_depth = 2, split_grid = 7, states = 3, beta = 1,
                               precision_range = c(0, 2), precision_points = 3),
                 density_prior(max_depth = 2, split_grid = 7, beta = 0))
  for (p in priors) {
    ref <- ref_fixed(t, p)
    cpp <- run_fixed(t, unclass(p))
    expect_equal(cpp$log_evidence, ref$log_evidence, tolerance = 1e-11)
    expect_equal(exp(cpp$tree$nodes$log_mass), ref$mass, tolerance = 1e-11)
    expect_equal(unname(cpp$tree$states), ref$states, tolerance = 1e-11)
    expect_equal(cpp$tree$nodes$theta_mean, ref$theta, tolerance = 1e-11)
    # Children's expected probabilities add up to their parent's.
    m <- exp(cpp$tree$nodes$log_mass)
    expect_equal(m[2] + m[3], m[1]); expect_equal(m[4] + m[5], m[2]); expect_equal(m[6] + m[7], m[3])
  }
})

test_that("tabulated and direct log-gamma values agree beyond the table", {
  t <- fixed_tree(c(5000, 1800, 3200, 1500, 300, 1100, 2100), c(3, 4, 5), 7)
  p <- unclass(density_prior(max_depth = 2, split_grid = 7))
  tab <- run_fixed(t, p); expect_gt(tab$lgamma_cap, 1000)
  p$lgamma_cap <- 0L; direct <- run_fixed(t, p); expect_equal(direct$lgamma_cap, 0L)
  ref <- ref_fixed(t, density_prior(max_depth = 2, split_grid = 7))
  expect_equal(tab$log_evidence, direct$log_evidence, tolerance = 1e-12)
  expect_equal(tab$log_evidence, ref$log_evidence, tolerance = 1e-10)
  expect_equal(exp(tab$tree$nodes$log_mass), ref$mass, tolerance = 1e-10)
})

test_that("the emission table reproduces the direct computation", {
  set.seed(1); x <- matrix(rbeta(600, 2, 5), ncol = 3)
  p <- unclass(density_prior(max_depth = 8)); ctl <- unclass(density_control(particles = 32))
  set.seed(2); a <- ptdensity:::ptd_smc_cpp(x, p, ctl)
  p$no_emission_cache <- TRUE
  set.seed(2); b <- ptdensity:::ptd_smc_cpp(x, p, ctl)
  expect_equal(a$forest, b$forest, tolerance = 1e-12)
  expect_equal(a$log_evidence, b$log_evidence, tolerance = 1e-12)
})

test_that("SMC targets the fully enumerated adaptive tree posterior", {
  skip_on_cran()
  x <- rbind(c(.1, .1), c(.2, .8), c(.3, .4), c(.4, .9), c(.8, .3), c(.9, .6),
             c(.1, .7), c(.2, .6), c(.4, .2), c(.6, .5), c(.7, .8), c(.15, .12))
  u <- rbind(c(.12, .11), c(.5, .5), c(.85, .45), c(.3, .7), c(.95, .95))
  for (spike in c(TRUE, FALSE)) {
    p <- density_prior(max_depth = 2, split_grid = 3, min_node = 2, eta = .1,
                       midpoint_spike = spike, beta = .5)
    ref <- ref_enumerate(x, p)
    expect_equal(sum(ref$weights), 1)
    keys <- vapply(ref$tables, ref_key, character(1))
    expect_false(anyDuplicated(keys) > 0)
    ref_dens <- ref_predict(ref, u)
    runs <- list(c(0, .5), c(1, .5)); if (spike) runs <- c(runs, list(c(1, 1)))
    for (run in runs) {
      set.seed(713)
      fit <- density_pt(x, support = cbind(c(0, 0), c(1, 1)), prior = p,
                        control = density_control(particles = 12000, ess_threshold = run[1],
                                                  resample_power = run[2]))
      expect_equal(fit$diagnostics$log_evidence, ref$log_evidence, tolerance = .025)
      fk <- forest_keys(fit$forest)
      expect_false(anyNA(match(fk, keys)))
      empirical <- numeric(length(keys))
      agg <- tapply(fit$forest$weight, match(fk, keys), sum)
      empirical[as.integer(names(agg))] <- agg
      expect_lt(sum(abs(empirical - ref$weights)) / 2, .07)
      expect_equal(predict(fit, u), ref_dens, tolerance = .03)
      if (run[1] == 1) expect_gt(fit$diagnostics$resamples, 0)
    }
  }
})

test_that("predictive density given the tree is exact for the enumerated trees", {
  # With one particle set to a single deterministic tree (no data split
  # possible beyond the root), prediction equals the enumerated value.
  x <- cbind(c(.1, .2, .3, .7, .9))
  p <- density_prior(max_depth = 1, split_grid = 1, min_node = 2)
  ref <- ref_enumerate(x, p)
  expect_length(ref$weights, 1)
  fit <- density_pt(x, support = c(0, 1), prior = p, control = density_control(particles = 4))
  u <- cbind(c(.25, .75, 0, 1, 1.2, -.1))
  expect_equal(predict(fit, u), ref_predict(ref, u), tolerance = 1e-12)
  expect_equal(fit$diagnostics$log_evidence, ref$log_evidence, tolerance = 1e-12)
})
