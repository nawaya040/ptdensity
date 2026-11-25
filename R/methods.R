#' Posterior mean density at new points
#'
#' Evaluates the estimated density, the posterior mean density of
#' [density_pt()], at given points.
#'
#' The posterior mean density averages, over the final trees and their
#' weights, the density that each tree gives to the region containing the
#' point. It is a probability density on the support: it integrates to one
#' over the support and is zero outside it.
#'
#' @param object Result of [density_pt()].
#' @param newdata Points at which to evaluate the density: a matrix or data
#'   frame with the variables of the fit as columns (matched by name when
#'   both have column names), or a numeric vector for a single variable.
#'   Without it, the density at the observations is returned (this needs
#'   `density_control(keep_data = TRUE)`, the default).
#' @param log If `TRUE`, the natural logarithm of the density is returned.
#' @param ... Not used.
#' @return A numeric vector with one value per row of `newdata`, in the
#'   original units of the variables (the density with respect to the
#'   product of their scales).
#' @seealso [marginal()] for the density of one or two variables.
#' @examples
#' fit <- density_pt(faithful, control = density_control(particles = 32))
#' predict(fit, data.frame(eruptions = c(2, 4.5), waiting = c(55, 80)))
#' # Average log density of held-out points
#' mean(predict(fit, faithful[1:20, ], log = TRUE))
#' @export
predict.ptdensity_fit <- function(object, newdata = NULL, log = FALSE, ...) {
  .pd_dots(...); .pd_flag(log, "log")
  x <- if (is.null(newdata)) .pd_need_data(object) else .pd_newdata(object, newdata)
  ld <- ptd_predict_cpp(object$forest, .pd_unit(x, object$support)) - object$forest$log_width
  if (log) ld else exp(ld)
}

#' Draw from the estimated distribution
#'
#' Draws observations from the predictive distribution of [density_pt()],
#' the distribution whose density is the posterior mean density.
#'
#' Each draw chooses one of the final trees with probability equal to its
#' weight, chooses a region of that tree with probability equal to the
#' posterior mean of its probability, and places the point uniformly in the
#' region.
#'
#' @param object Result of [density_pt()].
#' @param nsim Number of draws.
#' @param seed Optional seed passed to [set.seed()] before drawing.
#' @param ... Not used.
#' @return A data frame with `nsim` rows and one column per variable, in the
#'   original units.
#' @examples
#' fit <- density_pt(faithful, control = density_control(particles = 32))
#' draws <- simulate(fit, nsim = 500, seed = 1)
#' summary(draws)
#' @importFrom stats simulate
#' @export
simulate.ptdensity_fit <- function(object, nsim = 1, seed = NULL, ...) {
  .pd_dots(...); .pd_integer(nsim, "nsim", 1, .Machine$integer.max)
  if (!is.null(seed)) set.seed(seed)
  s <- object$support
  u <- ptd_simulate_cpp(object$forest, as.integer(nsim), nrow(s))
  x <- sweep(sweep(u, 2, s[, 2] - s[, 1], "*"), 2, s[, 1], "+")
  colnames(x) <- object$variables
  as.data.frame(x)
}

#' Density of one or two variables
#'
#' Computes the marginal density of one or two of the variables under the
#' posterior mean density of [density_pt()], on a grid of values. It is the
#' estimated density of those variables with the others integrated out, and
#' is computed exactly from the trees, without simulation.
#'
#' @param x Result of [density_pt()].
#' @param variables One or two variables, by name or by column number.
#' @param n_grid Number of grid values for each variable when `grid` is not
#'   given. The values are the midpoints of `n_grid` equal intervals of the
#'   support.
#' @param grid Optional grid values in the original units: a numeric vector
#'   for one variable, or a list of two numeric vectors for two variables.
#' @return A data frame with one column per variable and a column `density`.
#'   For two variables it has one row per combination of grid values, with
#'   the first variable varying fastest.
#' @seealso [predict()][predict.ptdensity_fit] for the joint density.
#' @examples
#' fit <- density_pt(faithful, control = density_control(particles = 32))
#' m <- marginal(fit, "eruptions")
#' plot(m$eruptions, m$density, type = "l")
#' m2 <- marginal(fit, c("eruptions", "waiting"), n_grid = 50)
#' z <- matrix(m2$density, 50)
#' image(unique(m2$eruptions), unique(m2$waiting), z)
#' @export
marginal <- function(x, variables, n_grid = 100L, grid = NULL) {
  .pd_fit(x)
  if (is.character(variables)) j <- match(variables, x$variables)
  else if (is.numeric(variables) && all(variables == floor(variables))) j <- as.integer(variables)
  else stop("'variables' must be variable names or column numbers.", call. = FALSE)
  if (!length(j) || length(j) > 2L || anyNA(j) || any(j < 1L | j > length(x$variables)) ||
      anyDuplicated(j))
    stop("'variables' must name one or two different variables of the fit.", call. = FALSE)
  s <- x$support[j, , drop = FALSE]
  if (is.null(grid)) {
    .pd_integer(n_grid, "n_grid", 1, 1e5)
    grid <- lapply(seq_along(j), function(a) {
      w <- (s[a, 2] - s[a, 1]) / n_grid
      s[a, 1] + w * (seq_len(n_grid) - .5)
    })
  } else {
    if (!is.list(grid)) grid <- list(grid)
    if (length(grid) != length(j) ||
        !all(vapply(grid, function(g) is.numeric(g) && length(g) >= 1L && all(is.finite(g)), logical(1))))
      stop("'grid' must give finite numeric values for each variable.", call. = FALSE)
  }
  ord <- lapply(grid, order)
  unit <- lapply(seq_along(j), function(a)
    .pd_unit(matrix(grid[[a]][ord[[a]]]), s[a, , drop = FALSE])[, 1])
  dens <- ptd_marginal_cpp(x$forest, as.integer(j - 1L), unit) / prod(s[, 2] - s[, 1])
  # Back to the order of the grid that was given.
  if (length(j) == 1L) {
    out <- data.frame(grid[[1]], density = dens[order(ord[[1]])])
  } else {
    m <- matrix(dens, length(grid[[1]]))[order(ord[[1]]), order(ord[[2]]), drop = FALSE]
    out <- expand.grid(grid[[1]], grid[[2]], KEEP.OUT.ATTRS = FALSE)
    out$density <- as.vector(m)
  }
  names(out)[seq_along(j)] <- x$variables[j]
  out
}

#' Regions of the representative tree
#'
#' Lists the regions of the representative tree, the sampled tree with the
#' highest posterior score. The tree shows where the data space was divided
#' and how the probability is spread over the resulting regions.
#'
#' @param x Result of [density_pt()].
#' @param which `"leaves"` (default): the final regions, which together cover
#'   the support. `"splits"`: the divisions, one per divided region.
#' @return A data frame with one row per region. For `"leaves"`:
#'   \describe{
#'     \item{`node_id`, `depth`}{Position of the region in the tree; depth 0
#'       is the whole support.}
#'     \item{`condition`}{The region, as conditions on the variables.}
#'     \item{`n`}{Number of observations in the region.}
#'     \item{`prob`}{Posterior mean probability of the region given the
#'       representative tree. The values sum to one over the leaves.}
#'     \item{`density`}{`prob` divided by the volume of the region in the
#'       original units: the posterior mean density in the region given the
#'       representative tree.}
#'   }
#'   For `"splits"`, the columns `node_id`, `depth`, `condition` and `n`
#'   describe the divided region, and:
#'   \describe{
#'     \item{`variable`, `threshold`}{The division: the left part is
#'       `variable <= threshold`.}
#'     \item{`n_left`}{Number of observations in the left part.}
#'     \item{`fraction`}{Proportion of the region's width on the left side.}
#'     \item{`prob_left`}{Posterior mean of the conditional probability of
#'       the left part given the region.}
#'     \item{`prob_flat`}{Posterior probability of the last state, in which
#'       the density is flat over the whole region.}
#'     \item{`state_1`, `state_2`, ...}{Posterior probabilities of the hidden
#'       states of the region (see [density_prior()]). Lower states allow
#'       larger departures from a flat density at this division.}
#'   }
#' @examples
#' fit <- density_pt(faithful, control = density_control(particles = 32))
#' head(partition(fit))
#' head(partition(fit, "splits"))
#' @export
partition <- function(x, which = c("leaves", "splits")) {
  .pd_fit(x); which <- match.arg(which)
  t <- x$tree
  if (which == "leaves") {
    out <- t[t$left == 0L, c("node_id", "depth", "condition", "n", "prob", "density")]
  } else {
    rows <- which(t$left > 0L)
    out <- t[rows, c("node_id", "depth", "condition", "n", "variable", "threshold"), drop = FALSE]
    out$n_left <- t$n[t$left[rows]]
    out$fraction <- t$fraction[rows]
    out$prob_left <- t$theta_mean[rows]
    out$prob_flat <- t$prob_flat[rows]
    states <- grep("^state_", names(t), value = TRUE)
    out <- cbind(out, t[rows, states, drop = FALSE])
  }
  rownames(out) <- NULL
  out
}

#' Computational diagnostics
#'
#' Information about the sequential Monte Carlo computation of a fit.
#'
#' @param x Result of [density_pt()].
#' @param ... Not used.
#' @return A list. `trace` has one row per step of the computation, with the
#'   effective sample size before resampling (`ess`), the number of distinct
#'   trees, the number of distinct initial particles still represented and
#'   whether the particles were resampled. Other components: the final
#'   particle `weights` and `particle_scores`, the final effective sample
#'   size `final_ess`, the number of `distinct_trees` kept for prediction,
#'   the number of `resamples`, the size in bytes of the kept trees
#'   (`forest_bytes`), the estimated log marginal likelihood relative to the
#'   uniform distribution on the support (`log_evidence`), the precision
#'   values used for the states (`precisions`), the regions of the
#'   representative tree that reached `max_depth` and the running time.
#' @examples
#' fit <- density_pt(faithful, control = density_control(particles = 32))
#' d <- diagnostics(fit)
#' d$final_ess
#' plot(d$trace$step, d$trace$ess, type = "l")
#' @export
diagnostics <- function(x, ...) UseMethod("diagnostics")

#' @rdname diagnostics
#' @export
diagnostics.ptdensity_fit <- function(x, ...) { .pd_dots(...); x$diagnostics }

#' @rdname diagnostics
#' @export
diagnostics.default <- function(x, ...) {
  if (is.list(x) && !is.null(x$diagnostics)) return(x$diagnostics)
  stop("No diagnostics available for this object.", call. = FALSE)
}

#' Summarize a density estimate
#'
#' `print(fit)` shows the size of the data and of the computation.
#' `summary(fit)` also shows the support, the model settings and a short
#' report of the computation.
#'
#' @param x,object Result of [density_pt()] (or of `summary()`).
#' @param ... Not used.
#' @return `summary()` returns a `summary.ptdensity_fit` object, which is
#'   printed. `print()` returns its argument invisibly.
#' @examples
#' fit <- density_pt(faithful, control = density_control(particles = 32))
#' fit
#' summary(fit)
#' @name summary.ptdensity_fit
NULL

#' @rdname summary.ptdensity_fit
#' @export
print.ptdensity_fit <- function(x, ...) {
  cat("Bayesian density estimate (ptdensity)\n")
  cat(sprintf("Observations: %d; variables: %d\n", x$n, length(x$variables)))
  cat(sprintf("Particles: %d; distinct final trees: %d; elapsed: %.2f s\n",
              x$control$particles, x$diagnostics$distinct_trees, x$diagnostics$elapsed_seconds))
  cat(sprintf("Regions in the representative tree: %d\n", sum(x$tree$left == 0L)))
  cat("Use predict(), marginal(), simulate() and partition() for details.\n")
  invisible(x)
}

#' @rdname summary.ptdensity_fit
#' @export
summary.ptdensity_fit <- function(object, ...) {
  .pd_dots(...)
  structure(list(n = object$n, variables = object$variables, support = object$support,
                 support_source = object$support_source, prior = object$prior,
                 control = object$control, diagnostics = object$diagnostics,
                 leaves = sum(object$tree$left == 0L),
                 max_depth_reached = max(object$tree$depth)),
            class = "summary.ptdensity_fit")
}

#' @rdname summary.ptdensity_fit
#' @export
print.summary.ptdensity_fit <- function(x, ...) {
  d <- x$diagnostics
  cat("Bayesian density estimate (ptdensity)\n")
  cat(sprintf("Observations: %d; variables: %d\n", x$n, length(x$variables)))
  cat(sprintf("Support (%s):\n", x$support_source))
  print(x$support)
  p <- x$prior
  cat(sprintf("Model: max_depth = %d, split_grid = %d, min_node = %d, eta = %g, midpoint_spike = %s\n",
              p$max_depth, p$split_grid, p$min_node, p$eta, p$midpoint_spike))
  cat(sprintf("       states = %d, beta = %g, log10 precision in (%g, %g]\n",
              p$states, p$beta, p$precision_range[1], p$precision_range[2]))
  cat(sprintf("Particles: %d; final ESS = %.1f; resampling steps = %d; distinct final trees = %d\n",
              x$control$particles, d$final_ess, d$resamples, d$distinct_trees))
  cat(sprintf("Representative tree: %d regions, deepest at depth %d\n",
              x$leaves, x$max_depth_reached))
  if (d$depth_limited_leaves > 0)
    cat(sprintf("Regions that reached max_depth: %d (%.1f%% of observations)\n",
                d$depth_limited_leaves, 100 * d$depth_limited_fraction))
  cat(sprintf("Log marginal likelihood (vs. uniform on the support): %.2f\n", d$log_evidence))
  cat(sprintf("Elapsed: %.2f s\n", d$elapsed_seconds))
  invisible(x)
}

#' @export
print.density_prior <- function(x, ...) {
  cat("ptdensity model settings\n"); utils::str(unclass(x), give.head = FALSE); invisible(x)
}

#' @export
print.density_control <- function(x, ...) {
  cat("ptdensity computational settings\n"); utils::str(unclass(x), give.head = FALSE); invisible(x)
}

.pd_fit <- function(x) {
  if (!inherits(x, "ptdensity_fit")) stop("'x' must be a ptdensity_fit object.", call. = FALSE)
}
.pd_dots <- function(...) {
  if (...length()) stop("Unused arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
}
.pd_need_data <- function(x) {
  if (is.null(x$data))
    stop("The data were not kept (density_control(keep_data = FALSE)); supply 'newdata'.", call. = FALSE)
  x$data
}
.pd_newdata <- function(object, newdata) {
  x <- .pd_matrix(newdata, "newdata")
  v <- object$variables
  if (!is.null(colnames(x)) && all(v %in% colnames(x))) x <- x[, v, drop = FALSE]
  else if (ncol(x) != length(v))
    stop(sprintf("'newdata' must have the %d variables of the fit: %s.", length(v),
                 paste(v, collapse = ", ")), call. = FALSE)
  colnames(x) <- v
  x
}
