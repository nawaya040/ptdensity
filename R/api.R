#' Computational settings
#'
#' Settings that control how [density_pt()] carries out its computation. They
#' affect accuracy and running time.
#'
#' @param particles Number of candidate trees computed in parallel (the
#'   particles of sequential Monte Carlo). More particles give a more accurate
#'   posterior mean density and take longer. Awaya and Ma (2024) used 1000.
#' @param ess_threshold Fraction of `particles` (between 0 and 1). When the
#'   effective number of particles falls below this fraction, the particles
#'   are resampled.
#' @param resample_power Number in (0, 1] used when resampling. Particles are
#'   drawn with probabilities proportional to their weights raised to this
#'   power. The default 0.5 follows Awaya and Ma (2024); 1 gives standard
#'   resampling.
#' @param keep_data If `TRUE`, the data are kept in the result, so that
#'   `predict(fit)` without `newdata` gives the density at the observations.
#' @return A `density_control` object to pass to [density_pt()].
#' @seealso [density_prior()] for the statistical model settings.
#' @examples
#' density_control()
#' density_control(particles = 1000)
#' @export
density_control <- function(particles = 256L, ess_threshold = 0.5,
                            resample_power = 0.5, keep_data = TRUE) {
  .pd_integer(particles, "particles", 2, 100000)
  .pd_number(ess_threshold, "ess_threshold", 0, 1)
  .pd_number(resample_power, "resample_power", 0, 1)
  if (resample_power == 0) stop("'resample_power' must be greater than 0.", call. = FALSE)
  .pd_flag(keep_data, "keep_data")
  structure(list(particles = as.integer(particles), ess_threshold = ess_threshold,
                 resample_power = resample_power, keep_data = keep_data),
            class = "density_control")
}

#' Model settings
#'
#' Settings of the statistical model used by [density_pt()]: how finely the
#' data space may be divided into regions, where the divisions may be placed,
#' and how the smoothness of the density may change from region to region.
#' The defaults are those of Awaya and Ma (2024, Section 5.1 and
#' Supplementary Section F.1).
#'
#' The density is built by dividing the data space into rectangular regions
#' step by step. At each division of a region, the probability of the left
#' part has a beta distribution centered at the proportion of the region's
#' volume on the left side. Its precision `nu` controls how far the density
#' can depart from flat at that division: a large `nu` keeps the division
#' close to proportional to volume. Each region has a hidden state that sets
#' the range of `nu`. With `states = I`, state `i < I` has `log10(nu)`
#' uniform on the `i`-th of `I - 1` equal parts of `precision_range`, and
#' state `I` has `nu` infinite, which makes the density flat within the
#' region. The state can only stay or increase from a region to its
#' subregions, with probabilities proportional to `exp(-beta * step)`, so the
#' density tends to become smoother in smaller regions.
#'
#' @param max_depth Maximum number of successive divisions (1 to 30). Larger
#'   values allow smaller regions and take longer.
#' @param split_grid Number of candidate division points for a variable within
#'   a region (an odd number from 1 to 127). The candidates are equally spaced
#'   and include the middle of the region.
#' @param min_node Regions with fewer observations than this are not divided
#'   further.
#' @param eta Non-negative number. Larger values favor division points closer
#'   to the middle of a region, more strongly in regions with many
#'   observations. `eta = 0` treats all candidates equally.
#' @param midpoint_spike If `TRUE` (default), once a region is divided at its
#'   middle, all its subregions are also divided at their middles, as in
#'   Awaya and Ma (2024, Appendix A). If `FALSE`, every candidate division
#'   point is allowed in every region.
#' @param states Number of hidden smoothness states (2 to 10), including the
#'   final state in which the density is flat.
#' @param beta Non-negative number. Larger values make the state less likely
#'   to increase from a region to its subregions.
#' @param precision_range Two numbers: the range of `log10(nu)` covered by the
#'   states other than the last.
#' @param precision_points Number of equally spaced values of `log10(nu)`
#'   used within each state to compute the uniform distribution of
#'   `log10(nu)` (1 to 50).
#' @return A `density_prior` object to pass to [density_pt()].
#' @seealso [density_control()] for the computational settings.
#' @examples
#' density_prior()
#' density_prior(max_depth = 10, eta = 0.01)
#' @export
density_prior <- function(max_depth = 15L, split_grid = 31L, min_node = 5L,
                          eta = 0.1, midpoint_spike = TRUE, states = 5L,
                          beta = 0.1, precision_range = c(-1, 4),
                          precision_points = 5L) {
  .pd_integer(max_depth, "max_depth", 1, 30)
  .pd_integer(split_grid, "split_grid", 1, 127)
  if (split_grid %% 2 != 1) stop("'split_grid' must be odd (to include the midpoint).", call. = FALSE)
  .pd_integer(min_node, "min_node", 2, .Machine$integer.max)
  .pd_number(eta, "eta", 0, 100)
  .pd_flag(midpoint_spike, "midpoint_spike")
  .pd_integer(states, "states", 2, 10)
  .pd_number(beta, "beta", 0, 100)
  if (!is.numeric(precision_range) || length(precision_range) != 2L ||
      any(!is.finite(precision_range)) || precision_range[1] >= precision_range[2] ||
      any(abs(precision_range) > 10))
    stop("'precision_range' must be two increasing finite numbers between -10 and 10.", call. = FALSE)
  .pd_integer(precision_points, "precision_points", 1, 50)
  structure(list(max_depth = as.integer(max_depth), split_grid = as.integer(split_grid),
                 min_node = as.integer(min_node), eta = eta,
                 midpoint_spike = midpoint_spike, states = as.integer(states),
                 beta = beta, precision_range = as.numeric(precision_range),
                 precision_points = as.integer(precision_points)),
            class = "density_prior")
}

#' Estimate a multivariate density
#'
#' Estimates the density of one sample of observations on one or more numeric
#' variables with the hidden Markov Pólya tree of Awaya and Ma (2024). The
#' estimate is the posterior mean of the density (also called the predictive
#' density).
#'
#' The method divides the data space into rectangular regions step by step,
#' building a tree of divisions, and gives each region a probability. Many
#' candidate trees are explored by sequential Monte Carlo, and the estimate
#' averages the densities of all of them with their posterior weights, so it
#' does not depend on a single choice of regions. The tree that best fits the
#' data (the representative tree) is kept to describe the regions; see
#' [partition()].
#'
#' @param x Data: a numeric matrix or data frame with one row per observation
#'   and one column per variable, or a numeric vector for a single variable.
#' @param support Optional known range of the variables: a matrix with one row
#'   per variable and two columns (lower and upper limit), or a vector
#'   `c(lower, upper)` for a single variable. Without it, the observed range
#'   of each variable is extended by 1% on each side. The estimated density is
#'   zero outside the support.
#' @param prior Model settings from [density_prior()].
#' @param control Computational settings from [density_control()].
#' @return A `ptdensity_fit` object. Use [predict()][predict.ptdensity_fit],
#'   [simulate()][simulate.ptdensity_fit], [marginal()], [partition()] and
#'   [summary()][summary.ptdensity_fit] to look at the results. Its main
#'   components are:
#'   \describe{
#'     \item{`tree`}{The representative tree, one row per region.
#'       [partition()] gives an easier-to-read table.}
#'     \item{`forest`}{The distinct final trees and their weights, used for
#'       the posterior mean density.}
#'     \item{`support`}{The support, one row per variable.}
#'     \item{`diagnostics`}{Information about the computation; see
#'       [diagnostics()].}
#'     \item{`n`, `variables`}{Number of observations and names of the
#'       variables.}
#'     \item{`prior`, `control`}{The settings used.}
#'     \item{`data`}{The observations; `NULL` when
#'       `density_control(keep_data = FALSE)` was used.}
#'   }
#' @details Missing (`NA`) and infinite values are not allowed; remove those
#'   rows first. A variable that takes the same value in all observations
#'   needs `support`. If there are fewer than `min_node` observations (see
#'   [density_prior()]), no division is possible and the estimate is the
#'   uniform density on the support.
#' @references Awaya, N. and Ma, L. (2024). Hidden Markov Pólya Trees for
#'   High-Dimensional Distributions. *Journal of the American Statistical
#'   Association*, 119(545), 189--201. \doi{10.1080/01621459.2022.2105223}
#'
#'   Ma, L. (2017). Adaptive Shrinkage in Pólya Tree Type Models. *Bayesian
#'   Analysis*, 12(3), 779--805. \doi{10.1214/16-BA1021}
#' @seealso [predict()][predict.ptdensity_fit], [marginal()],
#'   [simulate()][simulate.ptdensity_fit], [partition()],
#'   [density_prior()], [density_control()]
#' @examples
#' fit <- density_pt(faithful, control = density_control(particles = 32))
#' fit
#' predict(fit, data.frame(eruptions = c(2, 4.5), waiting = c(55, 80)))
#' head(marginal(fit, "waiting"))
#' @export
density_pt <- function(x, support = NULL, prior = density_prior(),
                       control = density_control()) {
  start <- proc.time()[[3L]]
  # Revalidate lists so editing a settings object cannot bypass checks.
  if (!inherits(prior, "density_prior") || !inherits(control, "density_control"))
    stop("Use density_prior() and density_control() to specify settings.", call. = FALSE)
  prior <- do.call(density_prior, unclass(prior))
  control <- do.call(density_control, unclass(control))
  input <- .pd_prepare(x, support)
  raw <- ptd_smc_cpp(input$unit, unclass(prior), unclass(control))
  log_width <- sum(log(input$support[, 2] - input$support[, 1]))
  tree <- .pd_tree(raw$tree, input, prior, log_width)
  if (!any(tree$left > 0))
    warning("No division was possible: there are fewer observations than 'min_node' ",
            "(see ?density_prior). The estimate is uniform on the support.", call. = FALSE)
  stopped <- tree$left == 0 & tree$depth == prior$max_depth & tree$n >= prior$min_node
  forest <- raw$forest
  forest$log_width <- log_width
  w <- raw$weights
  out <- list(call = match.call(), tree = tree, forest = forest,
              n = nrow(input$original), variables = input$variables,
              support = input$support, support_source = input$support_source,
              prior = prior, control = control,
              data = if (control$keep_data) input$original else NULL,
              diagnostics = list(trace = raw$trace, weights = w,
                particle_scores = raw$particle_scores, final_ess = 1 / sum(w^2),
                distinct_trees = length(forest$weight), resamples = raw$resamples,
                unique_proposals = raw$proposals,
                peak_particle_payload_bytes = raw$peak_particle_payload_bytes,
                forest_bytes = as.numeric(utils::object.size(forest)),
                log_evidence = raw$log_evidence, precisions = raw$precisions,
                depth_limited_leaves = sum(stopped),
                numerical_stop_leaves = sum(tree$stop_reason == "numerical precision"),
                depth_limited_fraction = sum(tree$n[stopped]) / nrow(input$original),
                elapsed_seconds = proc.time()[[3L]] - start))
  class(out) <- "ptdensity_fit"
  out
}

.pd_tree <- function(raw, input, prior, log_width) {
  tree <- raw$nodes
  states <- raw$states
  colnames(states) <- paste0("state_", seq_len(ncol(states)))
  internal <- tree$left > 0L
  tree$variable <- NA_character_
  tree$variable[internal] <- input$variables[tree$variable_id[internal]]
  tree$threshold <- NA_real_
  tree$threshold[internal] <- vapply(which(internal), function(i) {
    j <- tree$variable_id[i]
    thr <- input$support[j, 1] + tree$cut_unit[i] * (input$support[j, 2] - input$support[j, 1])
    .pd_snap(thr, tree$cut_unit[i], input$original[, j], input$unit[, j])
  }, numeric(1))
  tree$condition <- .pd_conditions(tree, input$original, input$variables)
  tree$prob <- exp(tree$log_mass)
  tree$density <- exp(tree$log_mass - tree$log_volume_unit - log_width)
  tree$prob_flat <- states[, ncol(states)]
  cbind(tree, as.data.frame(states))
}

.pd_number <- function(x, name, lower, upper, open = FALSE) {
  valid <- is.numeric(x) && !is.complex(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (valid) valid <- if (open) x > lower && x < upper else x >= lower && x <= upper
  if (!valid) stop(sprintf("'%s' must be a finite number %s %s and %s.",
                          name, if (open) "strictly between" else "between", lower, upper), call. = FALSE)
}
.pd_integer <- function(x, name, lower, upper) {
  .pd_number(x, name, lower, upper)
  if (x != floor(x)) stop(sprintf("'%s' must be an integer.", name), call. = FALSE)
}
.pd_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x))
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
}
.pd_matrix <- function(x, name, min_rows = 1L) {
  if (is.null(dim(x)) && is.numeric(x)) x <- matrix(x, ncol = 1L)
  if (is.data.frame(x)) {
    if (!all(vapply(x, function(z) is.numeric(z) && !is.complex(z), logical(1))))
      stop(sprintf("'%s' must contain only numeric columns.", name), call. = FALSE)
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || nrow(x) < min_rows || ncol(x) < 1L)
    stop(sprintf("'%s' must be a numeric matrix/data frame with at least %d row%s.",
                 name, min_rows, if (min_rows == 1L) "" else "s"), call. = FALSE)
  if (any(!is.finite(x))) stop(sprintf("'%s' contains missing or non-finite values.", name), call. = FALSE)
  storage.mode(x) <- "double"
  x
}
.pd_prepare <- function(x, support) {
  x <- .pd_matrix(x, "x", 2L)
  if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  if (anyNA(colnames(x)) || any(!nzchar(colnames(x))) || anyDuplicated(colnames(x)))
    stop("Variable names must be nonempty and unique.", call. = FALSE)
  ranges <- t(apply(x, 2, range))
  if (is.null(support)) {
    constant <- colnames(x)[ranges[, 1] == ranges[, 2]]
    if (length(constant))
      stop("Constant variables need 'support': ", paste(constant, collapse = ", "), call. = FALSE)
    widths <- ranges[, 2] - ranges[, 1]
    support <- cbind(ranges[, 1] - .01 * widths, ranges[, 2] + .01 * widths)
    source <- "observed range padded by 1%"
  } else {
    if (is.null(dim(support)) && length(support) == 2L) support <- matrix(support, 1L)
    if (!is.matrix(support) || !is.numeric(support) || is.complex(support) ||
        !identical(dim(support), c(ncol(x), 2L)) || any(!is.finite(support)) ||
        any(support[, 1] >= support[, 2]))
      stop("'support' must have finite lower < upper bounds, one row per variable.", call. = FALSE)
    if (any(ranges[, 1] < support[, 1] | ranges[, 2] > support[, 2]))
      stop("All observations must lie inside 'support'.", call. = FALSE)
    storage.mode(support) <- "double"
    source <- "user supplied"
  }
  dimnames(support) <- list(colnames(x), c("lower", "upper"))
  width <- support[, 2] - support[, 1]
  if (any(!is.finite(width)) || any(width <= 0))
    stop("Data scale is too extreme for safe range scaling; rescale the input.", call. = FALSE)
  unit <- .pd_unit(x, support)
  list(unit = unit, original = x, support = support, variables = colnames(x),
       support_source = source)
}
.pd_unit <- function(x, support) {
  u <- sweep(sweep(x, 2, support[, 1], "-"), 2, support[, 2] - support[, 1], "/")
  # Keep points on the support boundary inside [0, 1] despite rounding.
  u[x >= rep(support[, 1], each = nrow(x)) & u < 0] <- 0
  u[x <= rep(support[, 2], each = nrow(x)) & u > 1] <- 1
  u
}

# The engine splits on unit coordinates (unit <= cut goes left). Move the
# original-scale threshold, only if rounding requires it, so that
# original <= threshold reproduces exactly the same assignment.
.pd_snap <- function(thr, cut, original, unit) {
  left <- unit <= cut
  lmax <- if (any(left)) max(original[left]) else -Inf
  rmin <- if (any(!left)) min(original[!left]) else Inf
  if (thr >= lmax && thr < rmin) return(thr)
  if (is.finite(rmin) && is.finite(lmax)) (lmax + rmin) / 2 else thr
}
# Fewest digits that still divide the observations exactly as thr does.
.pd_short <- function(thr, original) {
  left <- original <= thr
  lmax <- if (any(left)) max(original[left]) else -Inf
  rmin <- if (any(!left)) min(original[!left]) else Inf
  for (k in seq_len(15)) {
    s <- signif(thr, k)
    if (s >= lmax && s < rmin) return(format(s, digits = 15, trim = TRUE))
  }
  format(thr, digits = 15, trim = TRUE)
}
# One clause per bound still in force, in the column order of the data.
.pd_conditions <- function(tree, original, variables) {
  n <- nrow(tree); d <- length(variables)
  lo <- matrix(-Inf, n, d); hi <- matrix(Inf, n, d)
  lo_text <- hi_text <- matrix("", n, d)
  label <- character(n)
  for (i in which(tree$left > 0))
    label[i] <- .pd_short(tree$threshold[i], original[, tree$variable_id[i]])
  for (j in seq_len(n)[-1L]) {
    par <- tree$parent_id[j]; v <- tree$variable_id[par]
    lo[j, ] <- lo[par, ]; hi[j, ] <- hi[par, ]
    lo_text[j, ] <- lo_text[par, ]; hi_text[j, ] <- hi_text[par, ]
    if (tree$left[par] == j) { hi[j, v] <- tree$threshold[par]; hi_text[j, v] <- label[par] }
    else { lo[j, v] <- tree$threshold[par]; lo_text[j, v] <- label[par] }
  }
  vapply(seq_len(n), function(j) {
    clauses <- character(0)
    for (v in seq_len(d)) {
      if (is.finite(lo[j, v])) clauses <- c(clauses, paste0(variables[v], " > ", lo_text[j, v]))
      if (is.finite(hi[j, v])) clauses <- c(clauses, paste0(variables[v], " <= ", hi_text[j, v]))
    }
    if (!length(clauses)) "whole support" else paste(clauses, collapse = " & ")
  }, character(1))
}
