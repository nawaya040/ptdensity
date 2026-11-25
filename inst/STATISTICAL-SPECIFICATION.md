# Statistical and computational specification (development version)

## Scope

One iid sample of continuous vectors is modeled. The estimate is the
posterior mean density (predictive density) E[q(x) | data] of the hidden
Markov Polya tree of Awaya and Ma (2024) with the adaptive-smoothness states
of their Example 1 (the adaptive Polya tree of Ma, 2017), as used in their
Section 5.1 and Supplementary Section F.1.

Each variable is rescaled linearly to [0,1] (the unit cube of the paper). The
support is supplied by the user or derived from the observed range with 1%
padding on each side; in the second case it is conditioned on the
observations. A variable that is constant in the data needs a supplied
support. The density on the original scale is the unit-scale density divided
by the product of the support widths, and is zero outside the support.

## Tree prior

The tree prior is that of ptcompare. The root is depth zero. Nodes at max_depth
or with fewer than min_node observations are terminal. Other nodes split along
a uniformly selected variable. Relative split fractions are k/(L+1),
k=1,...,L (L = split_grid, odd), with probabilities proportional to
exp(-eta*n_A*abs(k/(L+1)-1/2)): location prior (2) with N_L = L+1. By default
(midpoint_spike = TRUE) the spike-and-slab prior of Appendix A is used: once a
node is divided at its midpoint, its descendants are too. A node is also
terminal, and reported as such, if a side of its unit-scale rectangle is at
most 4*(L+1)*machine_epsilon. Observations with unit <= cut go left.
Terminal nodes carry the uniform distribution.

## State model

There are I states (states). Given V(A)=i < I, theta(A), the conditional
probability of the left child, is Beta(r*nu, (1-r)*nu) with r the split
fraction, and log10(nu) is uniform on (a(i), a(i+1)], where
a(i) = lo + (i-1)*(hi-lo)/(I-1) and (lo, hi) = precision_range (default
(-1, 4)). This uniform prior is approximated by K equally weighted points
(precision_points, default 5), placed at the midpoints of K equal
subintervals of (a(i), a(i+1)]. State I has nu = infinity: theta(A) = r, so
that, with I absorbing, the distribution within A is uniform.

The transition matrix does not depend on depth:
xi(i, i') proportional to exp(beta*(i - i')) for i <= i', and 0 for i > i'
(beta = 0.1, I = 5 by default). Structural zeros are kept as zeros. The root's
state has the distribution xi(1, .), that is, the root's parent is taken to be
in state 1.

Relation to the 2020 SMCMP software: the same model and defaults. The software
placed the K points from a(i) + 1/(2K) to a(i+1) - 1/(2K), which are the
midpoints only when a(i+1) - a(i) = 1; here they are always the midpoints. It
added a small constant to probabilities before taking logarithms; here zeros
stay zeros.

## Likelihood and messages

The internal emission in state i < I is

    (1/K) sum_k B(r nu_k + n_l, (1-r) nu_k + n_r) / B(r nu_k, (1-r) nu_k)
        * r^(-n_l) * (1-r)^(-n_r),

relative to the uniform distribution on the node; in state I it is one.
Terminal emissions are one. Own-state likelihoods U, upward messages M and the
evidence follow the ptcompare recursion with the I-state transition matrix,
all in log space. Downward (outside) messages give each node's exact state
posterior given the whole current tree.

## Predictive density

For a tree T, e_A(i) = E[Q(A) 1{V(A)=i} | data, T] (Proposition 3.2) is
computed from the root down:

    e_root(i) = xi(1,i) U_root(i) / sum_v xi(1,v) U_root(v)
    e_A(i')   = sum_i e_P(i) E[theta_side(P) | V(P)=i, data] xi~(i, i')
    xi~(i,i') = xi(i,i') U_A(i') / M_A(i)

where P is the parent of A and theta_side is theta(P) for a left child and
1 - theta(P) for a right child. Given V(P) = i < I, E[theta(P) | data] mixes
(r nu_k + n_l)/(nu_k + n) over k with weights proportional to the k-th term of
the emission; in state I it is r. E[Q(A) | data, T] = sum_i e_A(i). The
predictive density given T is E[Q(A) | data, T] / vol(A) on the leaf A that
contains x, and the reported density averages it over the final SMC particles
with their normalized weights. Final particles that share a tree are stored
once with their summed weight.

Draws from the predictive distribution choose a stored tree by weight, descend
it choosing each child with probability E[Q(child)]/E[Q(parent)], and draw
uniformly in the leaf; this samples exactly from the weighted mixture of the
conditional predictive measures. One- and two-variable marginal densities of
the predictive measure are computed exactly: the marginal of the uniform
distribution on a leaf is uniform on the projection of the leaf.

## SMC

As in ptcompare: breadth-first growth of the oldest eligible leaf, the
locally optimal proposal over all (variable, location) pairs using the leaf's
state distribution given the rest of the tree, the incremental weight equal to
the proposal's normalizing constant, kappa-tempered systematic resampling
(resample_power, 0.5 by default) when the ESS falls below
ess_threshold * particles, and the sum of log normalizing factors as the
log-evidence estimate. The representative tree maximizes log(tree prior) +
log(tree likelihood) among the final particles.

## Computation

Histograms with cumulative sums give candidate counts. Values of
lgamma(nu_k r_l + m) are tabulated for counts m up to a cap (at most 1024,
reduced so that the table has at most about 2 million entries), and computed
directly above it. lbeta(r nu_k, (1-r) nu_k) and lgamma(nu_k + n_A) are
computed once per location and per node. Emissions of nodes with fewer than
128 observations are tabulated once for all counts and locations (within a
budget of about 2 million entries); they are the same values computed in
advance. Particles share trees and
observation blocks as in ptcompare.

## Verification

tests/testthat/helper-exact.R enumerates all state configurations of small
fixed trees (evidence, state posteriors, E[Q(A)]) and all small adaptive trees
(tree posterior, evidence, predictive density), written independently of the
C++ code.

Reference: Awaya and Ma (2024), JASA 119(545), 189--201,
doi:10.1080/01621459.2022.2105223. The state model is that of Ma (2017),
Bayesian Analysis 12(3), 779--805, doi:10.1214/16-BA1021.
