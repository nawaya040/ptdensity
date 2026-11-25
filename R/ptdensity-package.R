#' ptdensity: Bayesian density estimation with hidden Markov Pólya trees
#'
#' Estimates multivariate densities with the hidden Markov Pólya tree method
#' of Awaya and Ma (2024), using the adaptive smoothness states of Ma (2017).
#'
#' A typical analysis:
#' 1. `fit <- density_pt(x)` estimates the density of the rows of `x`.
#' 2. `predict(fit, newdata)` evaluates the posterior mean density at new
#'    points.
#' 3. `marginal(fit, c("x1", "x2"))` gives the density of one or two
#'    variables on a grid.
#' 4. `simulate(fit, nsim)` draws from the estimated distribution.
#' 5. `partition(fit)` lists the regions of the representative tree.
#'
#' Settings are changed with [density_prior()] (statistical model) and
#' [density_control()] (computation).
#' @references Awaya, N. and Ma, L. (2024). Hidden Markov Pólya Trees for
#'   High-Dimensional Distributions. *Journal of the American Statistical
#'   Association*, 119(545), 189--201. \doi{10.1080/01621459.2022.2105223}
#'
#'   Ma, L. (2017). Adaptive Shrinkage in Pólya Tree Type Models. *Bayesian
#'   Analysis*, 12(3), 779--805. \doi{10.1214/16-BA1021}
#' @useDynLib ptdensity, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
