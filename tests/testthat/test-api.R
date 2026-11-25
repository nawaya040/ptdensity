make_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  cbind(a = c(rnorm(n / 2, -1, .4), rnorm(n / 2, 1, .7)), b = rnorm(n), c = runif(n))
}
small <- function(...) density_control(particles = 16, ...)
box <- cbind(c(-4, -4, 0), c(4, 4, 1))

test_that("the posterior mean density integrates to one and matches its marginals", {
  x <- make_data()
  set.seed(2); fit <- density_pt(x, support = box, control = small())
  expect_s3_class(fit, "ptdensity_fit")
  for (v in c("a", "b", "c")) {
    m <- marginal(fit, v, n_grid = 1e5)
    expect_equal(unname(sum(m$density) * diff(fit$support[v, ])) / 1e5, 1, tolerance = 1e-3)
  }
  m2 <- marginal(fit, c("a", "c"), n_grid = 1000)
  expect_equal(sum(m2$density) * 8 / 1e6, 1, tolerance = 1e-2)
  # A two-variable marginal of two-variable data is the joint density.
  fit2 <- density_pt(x[, 1:2], support = box[1:2, ], control = small())
  g <- marginal(fit2, 1:2, n_grid = 40)
  expect_equal(g$density, predict(fit2, as.matrix(g[, 1:2])), tolerance = 1e-12)
  # Marginal of one variable equals the joint density integrated over the other.
  m1 <- marginal(fit2, "a", grid = c(-1, .3, 1.2))
  b <- marginal(fit2, "b", n_grid = 1e4)$b
  joint <- sapply(c(-1, .3, 1.2), function(a) mean(predict(fit2, cbind(a, b))) * 8)
  expect_equal(m1$density, joint, tolerance = 1e-3)
})

test_that("draws follow the estimated distribution and stay in the support", {
  x <- make_data()
  set.seed(3); fit <- density_pt(x, support = box, control = small())
  s <- simulate(fit, 50000, seed = 4)
  expect_equal(dim(s), c(50000L, 3L)); expect_named(s, c("a", "b", "c"))
  expect_true(all(s$a >= -4 & s$a <= 4 & s$c >= 0 & s$c <= 1))
  m <- marginal(fit, "a", n_grid = 1e4)
  expect_equal(mean(s$a < 0), sum(m$density[m$a < 0]) * 8 / 1e4, tolerance = .02)
  expect_identical(simulate(fit, 10, seed = 5), simulate(fit, 10, seed = 5))
})

test_that("prediction handles names, order, support boundaries and outside points", {
  x <- make_data(200)
  set.seed(5); fit <- density_pt(x, support = box, control = small())
  new <- data.frame(c = c(.5, .2), b = c(0, 1), a = c(1, -1))
  expect_equal(predict(fit, new), predict(fit, as.matrix(new[, c("a", "b", "c")])))
  expect_equal(predict(fit, cbind(5, 0, .5)), 0)
  expect_equal(predict(fit, cbind(0, 0, 1.5), log = TRUE), -Inf)
  expect_gt(predict(fit, cbind(-4, -4, 0)), 0)
  expect_gt(predict(fit, cbind(4, 4, 1)), 0)
  expect_true(all(marginal(fit, "c", grid = c(0, 1))$density > 0))
  expect_equal(marginal(fit, "c", grid = c(-.1, 1.1))$density, c(0, 0))
  expect_length(predict(fit), 200)
  expect_equal(log(predict(fit)), predict(fit, log = TRUE))
  expect_error(predict(fit, cbind(1, 2)), "variables of the fit")
  fit0 <- density_pt(x, support = box, control = small(keep_data = FALSE))
  expect_null(fit0$data); expect_error(predict(fit0), "not kept")
})

test_that("results are reproducible and do not depend on row order", {
  x <- make_data(150)
  set.seed(9); f1 <- density_pt(x, control = small())
  set.seed(9); f2 <- density_pt(x, control = small())
  expect_identical(f1$forest, f2$forest)
  u <- x[c(3, 50, 120), ]
  set.seed(9); f3 <- density_pt(x[rev(seq_len(nrow(x))), ], control = density_control(particles = 4000))
  set.seed(10); f4 <- density_pt(x, control = density_control(particles = 4000))
  expect_equal(f3$diagnostics$log_evidence, f4$diagnostics$log_evidence, tolerance = .5)
})

test_that("too few observations give the uniform density on the support", {
  x <- cbind(u = c(.1, .4, .8))
  expect_warning(fit <- density_pt(x, support = c(0, 2), control = small()), "No division")
  expect_equal(predict(fit, cbind(c(.5, 1.9))), c(.5, .5))
  expect_equal(fit$diagnostics$log_evidence, 0)
  expect_equal(nrow(partition(fit, "splits")), 0)
  expect_equal(partition(fit)$prob, 1)
})

test_that("representative tree tables are consistent", {
  x <- make_data(300)
  set.seed(6); fit <- density_pt(x, control = small())
  leaves <- partition(fit)
  expect_equal(sum(leaves$prob), 1)
  expect_equal(sum(leaves$n), 300)
  sp <- partition(fit, "splits")
  expect_true(all(sp$prob_left > 0 & sp$prob_left < 1))
  st <- as.matrix(sp[, grep("^state_", names(sp))])
  expect_equal(unname(rowSums(st)), rep(1, nrow(sp)))
  expect_equal(sp$prob_flat, unname(st[, ncol(st)]))
  # Each threshold divides the region's observations as the engine did.
  for (i in seq_len(nrow(sp))) {
    node <- fit$tree[fit$tree$node_id == sp$node_id[i], ]
    child <- fit$tree[node$left, ]
    expect_equal(child$n, sp$n_left[i])
  }
  expect_output(print(fit), "density estimate")
  expect_output(print(summary(fit)), "Support")
  expect_identical(diagnostics(fit), fit$diagnostics)
  expect_equal(fit$diagnostics$distinct_trees, length(fit$forest$weight))
  expect_equal(sum(fit$forest$weight), 1)
})

test_that("inputs are validated", {
  x <- make_data(50)
  expect_error(density_pt(cbind(x, d = 1)), "Constant variables need 'support'")
  expect_error(density_pt(rbind(x, NA)), "non-finite")
  expect_error(density_pt(x, support = box[1:2, ]), "one row per variable")
  expect_error(density_pt(x, support = cbind(c(0, 0, 0), c(1, 1, 1))), "inside 'support'")
  expect_error(density_pt(x, prior = list()), "density_prior")
  expect_error(density_prior(split_grid = 4), "odd")
  expect_error(density_prior(states = 11), "states")
  expect_error(density_prior(precision_range = c(2, 1)), "precision_range")
  expect_error(density_control(resample_power = 0), "greater than 0")
  expect_error(marginal(density_pt(x, control = small()), c("a", "z")), "one or two")
  fit <- density_pt(x, control = small())
  expect_error(marginal(fit, 1:3), "one or two")
  m <- marginal(fit, "b", grid = c(1, -1, 0))
  expect_equal(m$b, c(1, -1, 0))
  expect_equal(m$density, marginal(fit, "b", grid = c(-1, 0, 1))$density[c(3, 1, 2)])
  m2 <- marginal(fit, c("b", "a"), grid = list(c(1, 0), c(.5, -.5, 0)))
  expect_equal(nrow(m2), 6); expect_named(m2, c("b", "a", "density"))
  ref <- marginal(fit, c("b", "a"), grid = list(c(0, 1), c(-.5, 0, .5)))
  expect_equal(m2$density, ref$density[c(6, 5, 2, 1, 4, 3)])
})

test_that("a constant variable is allowed with a supplied support", {
  x <- cbind(a = rnorm(60), b = 1)
  set.seed(1)
  fit <- density_pt(x, support = cbind(c(-5, 0), c(5, 2)), control = small())
  expect_gt(predict(fit, cbind(0, 1)), 0)
})
