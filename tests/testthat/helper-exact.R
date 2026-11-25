# Independent scalar R implementation of the model. Enumerates all
# internal-node state configurations and all small adaptive trees.
ref_lse <- function(v) { m <- max(v); if (is.infinite(m)) m else m + log(sum(exp(v - m))) }
ref_transition <- function(p) {
  I <- p$states; xi <- matrix(0, I, I)
  for (i in 1:I) { j <- i:I; xi[i, j] <- exp(p$beta * (i - j)); xi[i, ] <- xi[i, ] / sum(xi[i, ]) }
  xi
}
ref_precisions <- function(p) {
  I <- p$states; K <- p$precision_points; a <- p$precision_range
  w <- (a[2] - a[1]) / (I - 1)
  t(sapply(1:(I - 1), function(i) 10^(a[1] + (i - 1) * w + (1:K - .5) * w / K)))
}
# Log emission and posterior mean of theta in state s for nl/nr observations.
ref_state <- function(s, nl, nr, f, p) {
  if (s == p$states) return(list(em = 0, theta = f))
  nu <- ref_precisions(p)[s, , drop = TRUE]
  terms <- lbeta(nu * f + nl, nu * (1 - f) + nr) - lbeta(nu * f, nu * (1 - f))
  w <- exp(terms - ref_lse(terms))
  list(em = ref_lse(terms) - log(length(nu)) - nl * log(f) - nr * log(1 - f),
       theta = sum(w * (nu * f + nl) / (nu + nl + nr)))
}
# t: node table with left/right (1-based, 0 = leaf), parent_id, fraction, n.
ref_fixed <- function(t, p) {
  xi <- ref_transition(p); I <- p$states
  ids <- which(t$left > 0); N <- nrow(t)
  if (!length(ids)) return(list(log_evidence = 0, mass = 1, states = NULL, theta = NA))
  configs <- as.matrix(expand.grid(rep(list(1:I), length(ids))))
  logterms <- numeric(nrow(configs)); masses <- matrix(0, nrow(configs), N)
  thetas <- matrix(0, nrow(configs), N)
  for (h in seq_len(nrow(configs))) {
    state <- integer(N); state[ids] <- configs[h, ]
    lt <- 0; mass <- numeric(N); mass[1] <- 1
    for (j in ids) {
      pa <- t$parent_id[j]
      lt <- lt + log(if (pa == 0) xi[1, state[j]] else xi[state[pa], state[j]])
      l <- t$left[j]; r <- t$right[j]
      st <- ref_state(state[j], t$n[l], t$n[r], t$fraction[j], p)
      lt <- lt + st$em
      mass[l] <- mass[j] * st$theta; mass[r] <- mass[j] * (1 - st$theta)
      thetas[h, j] <- st$theta
    }
    logterms[h] <- lt; masses[h, ] <- mass
  }
  z <- ref_lse(logterms); w <- exp(logterms - z)
  states <- matrix(NA_real_, N, I)
  for (k in seq_along(ids)) for (s in 1:I) states[ids[k], s] <- sum(w[configs[, k] == s])
  theta <- rep(NA_real_, N); theta[ids] <- colSums(w * thetas)[ids]
  list(log_evidence = z, mass = colSums(w * masses), states = states, theta = theta)
}
# All adaptive trees on the unit cube for data x, with prior weights.
ref_enumerate <- function(x, p) {
  grow <- function(rows, lo, hi, depth, locked = FALSE) {
    base <- list(rows = rows, depth = depth, lo = lo, hi = hi)
    if (depth >= p$max_depth || length(rows) < p$min_node)
      return(list(list(node = base, logprior = 0)))
    fractions <- seq_len(p$split_grid) / (p$split_grid + 1)
    lp <- -p$eta * length(rows) * abs(fractions - .5)
    lp <- lp - ref_lse(lp) - log(ncol(x))
    if (locked) { keep <- fractions == .5; fractions <- fractions[keep]; lp <- -log(ncol(x)) }
    out <- list()
    for (j in seq_len(ncol(x))) for (k in seq_along(fractions)) {
      cut <- lo[j] + fractions[k] * (hi[j] - lo[j]); il <- rows[x[rows, j] <= cut]; ir <- setdiff(rows, il)
      hil <- hi; hil[j] <- cut; lor <- lo; lor[j] <- cut
      lock <- p$midpoint_spike && fractions[k] == .5
      ls <- grow(il, lo, hil, depth + 1, lock); rs <- grow(ir, lor, hi, depth + 1, lock)
      for (l in ls) for (r in rs) out[[length(out) + 1]] <- list(
        node = c(base, list(dim = j, fraction = fractions[k], cut = cut, left = l$node, right = r$node)),
        logprior = lp[k] + l$logprior + r$logprior)
    }
    out
  }
  flatten <- function(node) {
    queue <- list(node); parents <- 0L; out <- list(); i <- 1L
    while (i <= length(queue)) {
      a <- queue[[i]]; l <- r <- 0L
      if (!is.null(a$left)) {
        l <- length(queue) + 1L; r <- l + 1L
        queue <- c(queue, list(a$left, a$right)); parents <- c(parents, i, i)
      }
      out[[i]] <- data.frame(parent_id = parents[i], left = l, right = r, depth = a$depth,
        dim = if (is.null(a$dim)) NA_integer_ else a$dim,
        fraction = if (is.null(a$fraction)) NA_real_ else a$fraction,
        cut = if (is.null(a$cut)) NA_real_ else a$cut, n = length(a$rows),
        lo = I(list(a$lo)), hi = I(list(a$hi)))
      i <- i + 1L
    }
    do.call(rbind, out)
  }
  all <- grow(seq_len(nrow(x)), rep(0, ncol(x)), rep(1, ncol(x)), 0)
  tables <- lapply(all, function(a) flatten(a$node))
  fixed <- lapply(tables, ref_fixed, p = p)
  lp <- vapply(all, `[[`, numeric(1), "logprior") + vapply(fixed, `[[`, numeric(1), "log_evidence")
  z <- ref_lse(lp)
  list(tables = tables, fixed = fixed, weights = exp(lp - z), log_evidence = z)
}
# Predictive density on the unit cube at the rows of u.
ref_predict <- function(ref, u) {
  out <- numeric(nrow(u))
  for (k in seq_along(ref$tables)) {
    t <- ref$tables[[k]]; leaves <- which(t$left == 0)
    for (i in seq_len(nrow(u))) for (j in leaves) {
      lo <- t$lo[[j]]; hi <- t$hi[[j]]
      inside <- all((u[i, ] > lo | (lo == 0 & u[i, ] >= 0)) & u[i, ] <= hi)
      if (inside) out[i] <- out[i] + ref$weights[k] * ref$fixed[[k]]$mass[j] / prod(hi - lo)
    }
  }
  out
}
ref_key <- function(t) paste(ifelse(t$left > 0, paste(t$dim, round(t$cut, 8), t$left, sep = ":"), "-"),
                             collapse = "|")
forest_keys <- function(f) {
  vapply(seq_along(f$weight), function(k) {
    j <- (f$start[k] + 1):f$start[k + 1]; off <- f$start[k]
    paste(ifelse(f$left[j] >= 0, paste(f$dim[j] + 1, round(f$cut[j], 8), f$left[j] - off + 1, sep = ":"), "-"),
          collapse = "|")
  }, character(1))
}
