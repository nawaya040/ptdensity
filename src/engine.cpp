#include <Rcpp.h>
#include <array>
#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <limits>
using namespace Rcpp;
namespace ptd {
const int MAXS = 10;  // maximum number of states; the first I entries are used
const int MAXK = 50;  // maximum number of precision points per state
using SV = std::array<double, MAXS>;
const double neginf = -std::numeric_limits<double>::infinity();
double lse(const double* a, int n) {
  double m = *std::max_element(a, a + n);
  if (m == neginf) return m;
  double s = 0;
  for (int i = 0; i < n; ++i) s += std::exp(a[i] - m);
  return m + std::log(s);
}
double lsev(const std::vector<double>& a) {
  double m = *std::max_element(a.begin(), a.end()), s = 0;
  if (!std::isfinite(m)) Rcpp::stop("Non-finite particle or proposal weights.");
  for (double v : a) s += std::exp(v - m);
  return m + std::log(s);
}

struct Prior {
  int depth, grid, min_node, I, K, cap;
  bool spike;
  double eta, beta;
  SV root;
  std::array<SV, MAXS> trans;  // log transition probabilities, trans[from][to]
  std::vector<double> nu;      // (I-1)*K precisions, index s*K+k
  std::vector<double> base;    // lbeta(nu r_l, nu (1-r_l)), index (s*K+k)*grid+l
  std::vector<double> lg;      // lgamma(nu r_l + m), m = 0..cap
  Prior(const List& p, int n) {
    depth = as<int>(p["max_depth"]); grid = as<int>(p["split_grid"]);
    min_node = as<int>(p["min_node"]); spike = as<bool>(p["midpoint_spike"]);
    eta = as<double>(p["eta"]); beta = as<double>(p["beta"]);
    I = as<int>(p["states"]); K = as<int>(p["precision_points"]);
    NumericVector range = p["precision_range"];
    if (I < 2 || I > MAXS) stop("Invalid number of states.");
    if (K < 1 || K > MAXK) stop("Invalid number of precision points.");
    for (int v = 0; v < MAXS; ++v) trans[v].fill(neginf);
    for (int v = 0; v < I; ++v) {
      double z = neginf;
      for (int w = v; w < I; ++w) {
        trans[v][w] = beta * (v - w);
        double hi = std::max(z, trans[v][w]);
        z = hi + std::log(std::exp(z - hi) + std::exp(trans[v][w] - hi));
      }
      for (int w = v; w < I; ++w) trans[v][w] -= z;
    }
    root = trans[0];
    // log10(nu) | state s ~ Unif(a(s), a(s+1)], approximated by K midpoints.
    double width = (range[1] - range[0]) / (I - 1);
    for (int s = 0; s < I - 1; ++s)
      for (int k = 0; k < K; ++k)
        nu.push_back(std::pow(10.0, range[0] + s * width + (k + .5) * width / K));
    int SK = (I - 1) * K;
    base.resize(SK * grid);
    for (int sk = 0; sk < SK; ++sk)
      for (int l = 0; l < grid; ++l)
        base[sk * grid + l] = R::lbeta(nu[sk] * r(l), nu[sk] * r(grid - 1 - l));
    long budget = 2000000L / ((long)SK * grid);
    cap = (int)std::max(0L, std::min<long>(std::min(n, 1024), budget - 1));
    // Internal: tests lower the cap to compare tabulated and direct values.
    if (p.containsElementNamed("lgamma_cap")) cap = std::min(cap, as<int>(p["lgamma_cap"]));
    lg.resize((size_t)SK * grid * (cap + 1));
    for (int sk = 0; sk < SK; ++sk)
      for (int l = 0; l < grid; ++l)
        for (int m = 0; m <= cap; ++m)
          lg[((size_t)sk * grid + l) * (cap + 1) + m] = std::lgamma(nu[sk] * r(l) + m);
  }
  // Emissions of nodes with n < ncache observations, tabulated once
  // (offset[n] + (nl*grid + l)*(I-1) + s); filled by fill_cache().
  int ncache = 0;
  std::vector<size_t> offset;
  std::vector<double> ecache;
  double r(int l) const { return (l + 1.) / (grid + 1.); }
  // lgamma(nu_sk * r_l + m); r_{grid-1-l} = 1 - r_l gives the other side.
  double lgam(int sk, int l, int m) const {
    return m <= cap ? lg[((size_t)sk * grid + l) * (cap + 1) + m] : std::lgamma(nu[sk] * r(l) + m);
  }
};

// Log emission for every state of a node of n observations split at location l
// with nl on the left. lgn[sk] = lgamma(nu_sk + n). With weights != nullptr the
// posterior weights over the precision points are also returned (state-major).
void emission(int l, int nl, int n, const Prior& p, const std::vector<double>& lgn, SV& em,
              std::vector<double>* weights = nullptr) {
  double r = p.r(l), vol = -nl * std::log(r) - (n - nl) * std::log1p(-r);
  double tmp[MAXK];
  for (int s = 0; s < p.I - 1; ++s) {
    for (int k = 0; k < p.K; ++k) {
      int sk = s * p.K + k;
      tmp[k] = p.lgam(sk, l, nl) + p.lgam(sk, p.grid - 1 - l, n - nl) - lgn[sk] - p.base[sk * p.grid + l];
    }
    double z = lse(tmp, p.K);
    em[s] = z - std::log((double)p.K) + vol;
    if (weights) for (int k = 0; k < p.K; ++k) weights->push_back(std::exp(tmp[k] - z));
  }
  em[p.I - 1] = 0;
}
std::vector<double> lgamma_n(int n, const Prior& p) {
  std::vector<double> out(p.nu.size());
  for (size_t sk = 0; sk < p.nu.size(); ++sk) out[sk] = std::lgamma(p.nu[sk] + n);
  return out;
}
// Tabulates emissions for nodes with fewer than 128 observations, within a
// budget of about 2 million entries. Larger nodes are computed directly.
void fill_cache(Prior& p, int n, bool enabled) {
  long per = (long)p.grid * (p.I - 1), nc = enabled ? std::min(128, n + 1) : 0;
  while (nc > 0 && per * nc * (nc + 1) / 2 > 2000000L) --nc;
  p.ncache = nc; p.offset.assign(nc + 1, 0);
  for (int m = 0; m < nc; ++m) p.offset[m + 1] = p.offset[m] + (size_t)(m + 1) * per;
  p.ecache.resize(p.offset[nc]);
  SV em;
  for (int m = 0; m < nc; ++m) {
    std::vector<double> lgn = lgamma_n(m, p);
    for (int nl = 0; nl <= m; ++nl)
      for (int l = 0; l < p.grid; ++l) {
        emission(l, nl, m, p, lgn, em);
        std::copy(em.begin(), em.begin() + p.I - 1, p.ecache.begin() + p.offset[m] + (nl * p.grid + l) * (p.I - 1));
      }
  }
}
void cached_emission(int l, int nl, int n, const Prior& p, const std::vector<double>& lgn, SV& em) {
  if (n >= p.ncache) { emission(l, nl, n, p, lgn, em); return; }
  const double* src = p.ecache.data() + p.offset[n] + (nl * p.grid + l) * (p.I - 1);
  std::copy(src, src + p.I - 1, em.begin()); em[p.I - 1] = 0;
}
// Posterior mean of theta (left-child probability) given each state.
SV theta_mean(int l, int nl, int n, const Prior& p) {
  SV em, out; out.fill(0); std::vector<double> w;
  emission(l, nl, n, p, lgamma_n(n, p), em, &w);
  for (int s = 0; s < p.I - 1; ++s)
    for (int k = 0; k < p.K; ++k) {
      double nu = p.nu[s * p.K + k];
      out[s] += w[s * p.K + k] * (nu * p.r(l) + nl) / (nu + n);
    }
  out[p.I - 1] = p.r(l);
  return out;
}

struct Node {
  int parent = -1, left = -1, right = -1, depth = 0, dim = -1, n = 0, loc = -1;
  bool numerical_stop = false, midpoint = false;
  double cut = 0, fraction = .5;
  SV emission{}, up{}, message{};
  std::shared_ptr<std::vector<int>> obs;
};
struct Tree { std::vector<Node> nodes; size_t cursor = 0; double logprior = 0; };
SV message(const SV& u, const Prior& p) {
  SV out{}; double tmp[MAXS];
  for (int v = 0; v < p.I; ++v) {
    for (int w = 0; w < p.I; ++w) tmp[w] = p.trans[v][w] + u[w];
    out[v] = lse(tmp, p.I);
  }
  return out;
}
double evidence(const Tree& t, const Prior& p) {
  double tmp[MAXS];
  for (int v = 0; v < p.I; ++v) tmp[v] = p.root[v] + t.nodes[0].up[v];
  return lse(tmp, p.I);
}
void refresh(Tree& t, int id, const Prior& p) {
  while (id >= 0) {
    Node& a = t.nodes[id]; a.up = a.emission;
    if (a.left >= 0) for (int v = 0; v < p.I; ++v)
      a.up[v] += t.nodes[a.left].message[v] + t.nodes[a.right].message[v];
    a.message = message(a.up, p); id = a.parent;
  }
}
// State distribution at a node with that node's subtree removed.
SV outside(const Tree& t, int id, const Prior& p) {
  std::vector<int> path;
  for (int j = id; t.nodes[j].parent >= 0; j = t.nodes[j].parent) path.push_back(j);
  SV out = p.root; double tmp[MAXS], terms[MAXS];
  for (auto it = path.rbegin(); it != path.rend(); ++it) {
    int child = *it; const Node& par = t.nodes[t.nodes[child].parent];
    int sib = par.left == child ? par.right : par.left;
    for (int v = 0; v < p.I; ++v) tmp[v] = out[v] + par.emission[v] + t.nodes[sib].message[v];
    double z = lse(tmp, p.I); for (int v = 0; v < p.I; ++v) tmp[v] -= z;
    for (int w = 0; w < p.I; ++w) {
      for (int v = 0; v < p.I; ++v) terms[v] = tmp[v] + p.trans[v][w];
      out[w] = lse(terms, p.I);
    }
  }
  double z = lse(out.data(), p.I); for (int v = 0; v < p.I; ++v) out[v] -= z;
  return out;
}
// log E[Q(A) | data, T] for every node (Awaya and Ma, 2024, Proposition 3.2).
std::vector<double> log_mass(const Tree& t, const Prior& p) {
  int N = t.nodes.size(), I = p.I;
  std::vector<SV> e(N); std::vector<double> out(N);
  double z = evidence(t, p);
  for (int v = 0; v < I; ++v) e[0][v] = p.root[v] + t.nodes[0].up[v] - z;
  double tmp[MAXS];
  for (int j = 0; j < N; ++j) {
    const Node& a = t.nodes[j];
    out[j] = lse(e[j].data(), I);
    if (a.left < 0) continue;
    SV m = theta_mean(a.loc, t.nodes[a.left].n, a.n, p);
    for (int side = 0; side < 2; ++side) {
      const Node& c = t.nodes[side == 0 ? a.left : a.right];
      SV& ec = e[side == 0 ? a.left : a.right];
      for (int w = 0; w < I; ++w) {
        for (int v = 0; v < I; ++v)
          tmp[v] = e[j][v] + std::log(side == 0 ? m[v] : 1 - m[v]) + p.trans[v][w] + c.up[w] - c.message[v];
        ec[w] = lse(tmp, I);
      }
    }
  }
  return out;
}
void bounds(const Tree& t, int id, int d, std::vector<double>& lo, std::vector<double>& hi) {
  lo.assign(d, 0); hi.assign(d, 1);
  for (int j = id; t.nodes[j].parent >= 0; j = t.nodes[j].parent) {
    const Node& par = t.nodes[t.nodes[j].parent];
    if (par.left == j) hi[par.dim] = std::min(hi[par.dim], par.cut);
    else lo[par.dim] = std::max(lo[par.dim], par.cut);
  }
}
struct Option { int dim, nl, loc; bool midpoint; double cut, logprior; SV em; };
struct Proposal { int id = -1; double increment = 0; std::vector<Option> options; std::vector<double> cdf; };
Proposal propose(Tree& t, const NumericMatrix& x, const Prior& p) {
  Proposal prop;
  std::vector<double> lo, hi;
  while (t.cursor < t.nodes.size()) {
    Node& a = t.nodes[t.cursor];
    if (a.depth < p.depth && a.n >= p.min_node) {
      bounds(t, t.cursor, x.ncol(), lo, hi);
      bool safe = true;
      for (int j = 0; j < x.ncol(); ++j)
        if (hi[j] - lo[j] <= 4 * (p.grid + 1) * std::numeric_limits<double>::epsilon()) safe = false;
      if (safe) break;
      a.numerical_stop = true;
    }
    ++t.cursor;
  }
  if (t.cursor == t.nodes.size()) return prop;
  prop.id = t.cursor;
  const Node& a = t.nodes[prop.id]; int d = x.ncol(), L = p.grid;
  SV pi = outside(t, prop.id, p);
  std::vector<double> lp(L), logq, lgn = lgamma_n(a.n, p);
  for (int k = 0; k < L; ++k) lp[k] = -p.eta * a.n * std::abs((k + 1.) / (L + 1.) - .5);
  double z = lsev(lp); for (double& v : lp) v -= z + std::log((double)d);
  // Awaya and Ma (2024), Appendix A: R(A)=1 iff L(A)=1/2, and R is absorbing.
  const int mid = (L - 1) / 2;
  bool locked = p.spike && a.parent >= 0 && t.nodes[a.parent].midpoint;
  if (locked) for (int k = 0; k < L; ++k) lp[k] = k == mid ? -std::log((double)d) : neginf;
  prop.options.reserve(d * L); logq.reserve(d * L);
  double terms[MAXS];
  for (int j = 0; j < d; ++j) {
    if (j % 8 == 0) checkUserInterrupt();
    std::vector<int> h(L + 1, 0);
    std::vector<double> cuts(L);
    for (int k = 0; k < L; ++k) cuts[k] = lo[j] + (k + 1.) / (L + 1.) * (hi[j] - lo[j]);
    for (int row : *a.obs) {
      double value = x(row, j), scaled = (value - lo[j]) / (hi[j] - lo[j]) * (L + 1.);
      int bin = std::max(0, std::min(L, (int)std::floor(scaled)));
      // Correct floating-point boundary ambiguity using actual cut values.
      while (bin > 0 && value <= cuts[bin - 1]) --bin;
      while (bin < L && value > cuts[bin]) ++bin;
      h[bin]++;
    }
    int nl = 0;
    for (int k = 0; k < L; ++k) {
      nl += h[k];
      if (lp[k] == neginf) continue;
      SV em; cached_emission(k, nl, a.n, p, lgn, em);
      for (int v = 0; v < p.I; ++v) terms[v] = pi[v] + em[v];
      prop.options.push_back({j, nl, k, k == mid, cuts[k], lp[k], em});
      logq.push_back(lp[k] + lse(terms, p.I));
    }
  }
  prop.increment = lsev(logq); double cum = 0;
  for (double v : logq) { cum += std::exp(v - prop.increment); prop.cdf.push_back(cum); }
  prop.cdf.back() = 1; return prop;
}
void split(Tree& t, const Proposal& prop, int choice, const NumericMatrix& x, const Prior& p) {
  const Option& o = prop.options[choice]; int id = prop.id, left = t.nodes.size();
  Node a = t.nodes[id], l, r;
  l.parent = r.parent = id; l.depth = r.depth = a.depth + 1;
  l.n = o.nl; r.n = a.n - o.nl;
  l.obs = std::make_shared<std::vector<int>>(); r.obs = std::make_shared<std::vector<int>>();
  l.obs->reserve(l.n); r.obs->reserve(r.n);
  for (int row : *a.obs) (x(row, o.dim) <= o.cut ? *l.obs : *r.obs).push_back(row);
  if ((int)l.obs->size() != l.n) stop("Internal split count mismatch.");
  a.left = left; a.right = left + 1; a.dim = o.dim; a.cut = o.cut; a.loc = o.loc;
  a.fraction = p.r(o.loc); a.midpoint = o.midpoint; a.emission = o.em; a.obs.reset();
  t.nodes[id] = a;
  t.nodes.push_back(l); t.nodes.push_back(r); ++t.cursor;
  t.logprior += o.logprior; refresh(t, id, p);
}
List table(const Tree& t, const Prior& p) {
  int n = t.nodes.size(), I = p.I;
  IntegerVector id(n), parent(n), left(n), right(n), depth(n), dim(n), count(n);
  NumericVector cut(n), frac(n), tm(n), lv(n);
  CharacterVector reason(n);
  NumericMatrix states(n, I);
  std::vector<double> lm = log_mass(t, p);
  for (int j = 0; j < n; ++j) {
    const Node& a = t.nodes[j]; id[j] = j + 1; parent[j] = a.parent + 1; left[j] = a.left + 1;
    right[j] = a.right + 1; depth[j] = a.depth; count[j] = a.n;
    lv[j] = j == 0 ? 0 : lv[a.parent] + std::log(t.nodes[a.parent].left == j ?
      t.nodes[a.parent].fraction : 1 - t.nodes[a.parent].fraction);
    bool internal = a.left >= 0;
    dim[j] = internal ? a.dim + 1 : NA_INTEGER; cut[j] = internal ? a.cut : NA_REAL;
    frac[j] = internal ? a.fraction : NA_REAL;
    reason[j] = internal ? "split" : (a.numerical_stop ? "numerical precision" :
      (a.depth >= p.depth ? "max_depth" : "min_node"));
    tm[j] = NA_REAL;
    for (int v = 0; v < I; ++v) states(j, v) = NA_REAL;
    if (internal) {
      SV post = outside(t, j, p); for (int v = 0; v < I; ++v) post[v] += a.up[v];
      double z = lse(post.data(), I);
      SV m = theta_mean(a.loc, t.nodes[a.left].n, a.n, p); tm[j] = 0;
      for (int v = 0; v < I; ++v) { states(j, v) = std::exp(post[v] - z); tm[j] += states(j, v) * m[v]; }
    }
  }
  DataFrame df = DataFrame::create(_["node_id"] = id, _["parent_id"] = parent, _["left"] = left,
      _["right"] = right, _["depth"] = depth, _["variable_id"] = dim, _["cut_unit"] = cut,
      _["fraction"] = frac, _["n"] = count, _["log_mass"] = wrap(lm), _["log_volume_unit"] = lv,
      _["theta_mean"] = tm, _["stop_reason"] = reason);
  return List::create(_["nodes"] = df, _["states"] = states);
}
struct Particle { std::shared_ptr<Tree> tree; double logw; int ancestor; };
size_t payload(const std::vector<Particle>& ps) {
  std::unordered_set<const Tree*> ts; std::unordered_set<const std::vector<int>*> blocks;
  size_t bytes = 0;
  for (const auto& q : ps) if (ts.insert(q.tree.get()).second) {
    bytes += sizeof(Tree) + q.tree->nodes.capacity() * sizeof(Node);
    for (const auto& n : q.tree->nodes) if (n.obs && blocks.insert(n.obs.get()).second)
      bytes += sizeof(std::vector<int>) + sizeof(int) * n.obs->capacity();
  }
  return bytes;
}
// Compact storage of the distinct final trees for prediction.
struct Forest {
  std::vector<int> start, left, right, dim;
  std::vector<double> weight, cut, lmass, ldens;
  void add(const Tree& t, double w, const Prior& p) {
    int offset = left.size(); start.push_back(offset); weight.push_back(w);
    std::vector<double> lm = log_mass(t, p), lv(t.nodes.size(), 0);
    for (size_t j = 0; j < t.nodes.size(); ++j) {
      const Node& a = t.nodes[j];
      if (j > 0) { const Node& par = t.nodes[a.parent];
        lv[j] = lv[a.parent] + std::log(par.left == (int)j ? par.fraction : 1 - par.fraction); }
      left.push_back(a.left < 0 ? -1 : a.left + offset); right.push_back(a.right < 0 ? -1 : a.right + offset);
      dim.push_back(a.dim); cut.push_back(a.cut); lmass.push_back(lm[j]); ldens.push_back(lm[j] - lv[j]);
    }
  }
  List wrap_list() const {
    std::vector<int> s = start; s.push_back(left.size());
    return List::create(_["start"] = s, _["left"] = left, _["right"] = right, _["dim"] = dim,
                        _["cut"] = cut, _["log_mass"] = lmass, _["log_density"] = ldens,
                        _["weight"] = weight);
  }
};
} // namespace ptd

// [[Rcpp::export]]
List ptd_smc_cpp(NumericMatrix x, List prior, List control) {
  using namespace ptd;
  Prior p(prior, x.nrow()); int M = as<int>(control["particles"]);
  fill_cache(p, x.nrow(), !prior.containsElementNamed("no_emission_cache"));
  double threshold = as<double>(control["ess_threshold"]);
  double kappa = as<double>(control["resample_power"]);
  auto initial = std::make_shared<Tree>(); Node root;
  root.obs = std::make_shared<std::vector<int>>(); root.obs->reserve(x.nrow());
  for (int i = 0; i < x.nrow(); ++i) root.obs->push_back(i);
  root.n = x.nrow();
  initial->nodes.push_back(root);
  std::vector<Particle> ps;
  for (int i = 0; i < M; ++i) ps.push_back({initial, -std::log((double)M), i});
  initial.reset(); root.obs.reset();
  std::vector<double> ess_trace; std::vector<int> unique_trace, ancestor_trace, resampled_trace;
  double logZ = 0; size_t peak = payload(ps); int rounds = 0, resamples = 0, proposals = 0;
  while (true) {
    checkUserInterrupt();
    std::unordered_map<Tree*, size_t> lookup;
    std::vector<std::vector<int>> groups;
    for (int i = 0; i < M; ++i) {
      Tree* key = ps[i].tree.get(); auto hit = lookup.find(key);
      if (hit == lookup.end()) { lookup[key] = groups.size(); groups.push_back({i}); }
      else groups[hit->second].push_back(i);
    }
    bool changed = false;
    for (const auto& ids : groups) {
      Proposal q = propose(*ps[ids[0]].tree, x, p);
      if (q.id < 0) continue;
      changed = true; ++proposals;
      if (ids.size() == 1) {
        int i = ids[0]; int c = std::lower_bound(q.cdf.begin(), q.cdf.end(), R::runif(0, 1)) - q.cdf.begin();
        split(*ps[i].tree, q, c, x, p); ps[i].logw += q.increment;
      } else {
        auto source = ps[ids[0]].tree;
        std::unordered_map<int, std::shared_ptr<Tree>> children;
        for (int i : ids) {
          int c = std::lower_bound(q.cdf.begin(), q.cdf.end(), R::runif(0, 1)) - q.cdf.begin();
          auto hit = children.find(c);
          if (hit == children.end()) {
            auto child = std::make_shared<Tree>(*source); split(*child, q, c, x, p);
            children[c] = child; ps[i].tree = child;
          } else ps[i].tree = hit->second;
          ps[i].logw += q.increment;
        }
      }
    }
    if (!changed) break;
    ++rounds; std::vector<double> logw; for (const auto& q : ps) logw.push_back(q.logw);
    double z = lsev(logw), sq = 0; logZ += z;
    for (auto& q : ps) { q.logw -= z; sq += std::exp(2 * q.logw); }
    double ess = 1 / sq;
    bool resample = ess < threshold * M;
    if (resample) {
      // Section 3.1 (after Lu et al., 2013): ancestors are drawn with
      // probabilities a proportional to W^kappa and receive weights W/a.
      ++resamples; std::vector<Particle> next; next.reserve(M);
      std::vector<double> la(M);
      for (int i = 0; i < M; ++i) la[i] = kappa * ps[i].logw;
      double za = lsev(la); for (double& v : la) v -= za;
      double u = R::runif(0, 1) / M, cum = std::exp(la[0]); int j = 0;
      for (int i = 0; i < M; ++i) {
        double target = u + (double)i / M;
        while (j < M - 1 && target > cum) cum += std::exp(la[++j]);
        next.push_back({ps[j].tree, ps[j].logw - la[j] - std::log((double)M), ps[j].ancestor});
      }
      ps.swap(next);
      std::vector<double> lw; for (const auto& q : ps) lw.push_back(q.logw);
      double zr = lsev(lw); logZ += zr; for (auto& q : ps) q.logw -= zr;
    }
    std::unordered_set<Tree*> unique; std::unordered_set<int> ancestors;
    for (const auto& q : ps) { unique.insert(q.tree.get()); ancestors.insert(q.ancestor); }
    ess_trace.push_back(ess); unique_trace.push_back(unique.size());
    ancestor_trace.push_back(ancestors.size()); resampled_trace.push_back(resample);
    if ((rounds & (rounds - 1)) == 0) peak = std::max(peak, payload(ps));
  }
  peak = std::max(peak, payload(ps));
  double bestscore = neginf; int best = 0;
  NumericVector weights(M), scores(M);
  std::unordered_map<Tree*, int> index; std::vector<Tree*> distinct; std::vector<double> tw;
  for (int i = 0; i < M; ++i) {
    Tree* t = ps[i].tree.get(); weights[i] = std::exp(ps[i].logw);
    auto hit = index.find(t);
    if (hit == index.end()) { index[t] = distinct.size(); distinct.push_back(t); tw.push_back(weights[i]); }
    else tw[hit->second] += weights[i];
    scores[i] = t->logprior + evidence(*t, p);
    if (scores[i] > bestscore) { bestscore = scores[i]; best = i; }
  }
  Forest f;
  for (size_t k = 0; k < distinct.size(); ++k) { checkUserInterrupt(); f.add(*distinct[k], tw[k], p); }
  DataFrame trace = DataFrame::create(_["step"] = seq_len(rounds), _["ess"] = ess_trace,
      _["unique_trees"] = unique_trace, _["root_ancestors"] = ancestor_trace, _["resampled"] = resampled_trace);
  return List::create(_["tree"] = table(*ps[best].tree, p), _["forest"] = f.wrap_list(),
      _["weights"] = weights, _["particle_scores"] = scores, _["trace"] = trace,
      _["log_evidence"] = logZ, _["resamples"] = resamples, _["proposals"] = proposals,
      _["peak_particle_payload_bytes"] = (double)peak, _["lgamma_cap"] = p.cap,
      _["precisions"] = wrap(p.nu));
}

// Log predictive density (unit scale) at the rows of u; -Inf outside [0,1]^d.
// [[Rcpp::export]]
NumericVector ptd_predict_cpp(List forest, NumericMatrix u) {
  IntegerVector start = forest["start"], left = forest["left"], right = forest["right"], dim = forest["dim"];
  NumericVector cut = forest["cut"], ld = forest["log_density"], w = forest["weight"];
  int n = u.nrow(), d = u.ncol(), T = w.size();
  std::vector<double> m(n, ptd::neginf), s(n, 0);
  std::vector<bool> inside(n, true);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < d; ++j) if (!(u(i, j) >= 0 && u(i, j) <= 1)) inside[i] = false;
  for (int t = 0; t < T; ++t) {
    if (t % 16 == 0) checkUserInterrupt();
    if (w[t] <= 0) continue;
    double lw = std::log(w[t]);
    for (int i = 0; i < n; ++i) {
      if (!inside[i]) continue;
      int j = start[t];
      while (left[j] >= 0) j = u(i, dim[j]) <= cut[j] ? left[j] : right[j];
      double v = lw + ld[j];
      if (v == ptd::neginf) continue;
      if (v > m[i]) { s[i] = s[i] * std::exp(m[i] - v) + 1; m[i] = v; }
      else s[i] += std::exp(v - m[i]);
    }
  }
  NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = inside[i] && s[i] > 0 ? m[i] + std::log(s[i]) : R_NegInf;
  return out;
}

// Exact draws (unit scale) from the weighted mixture of predictive measures.
// [[Rcpp::export]]
NumericMatrix ptd_simulate_cpp(List forest, int nsim, int d) {
  IntegerVector start = forest["start"], left = forest["left"], right = forest["right"], dim = forest["dim"];
  NumericVector cut = forest["cut"], lm = forest["log_mass"], w = forest["weight"];
  std::vector<double> cdf(w.size()); double cum = 0;
  for (int t = 0; t < w.size(); ++t) { cum += w[t]; cdf[t] = cum; }
  NumericMatrix out(nsim, d); std::vector<double> lo(d), hi(d);
  for (int i = 0; i < nsim; ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    int t = std::lower_bound(cdf.begin(), cdf.end(), R::unif_rand() * cum) - cdf.begin();
    t = std::min<int>(t, w.size() - 1);
    std::fill(lo.begin(), lo.end(), 0.); std::fill(hi.begin(), hi.end(), 1.);
    int j = start[t];
    while (left[j] >= 0) {
      if (R::unif_rand() < std::exp(lm[left[j]] - lm[j])) { hi[dim[j]] = cut[j]; j = left[j]; }
      else { lo[dim[j]] = cut[j]; j = right[j]; }
    }
    for (int k = 0; k < d; ++k) out(i, k) = lo[k] + R::unif_rand() * (hi[k] - lo[k]);
  }
  return out;
}

// Exact one- or two-variable marginal density (unit scale) on a grid. Each
// leaf adds its mass over the grid points of its projection; the additions are
// accumulated as differences and summed at the end.
// [[Rcpp::export]]
NumericVector ptd_marginal_cpp(List forest, IntegerVector dims, List grids) {
  IntegerVector start = forest["start"], left = forest["left"], right = forest["right"], dim = forest["dim"];
  NumericVector cut = forest["cut"], lm = forest["log_mass"], w = forest["weight"];
  int q = dims.size();
  std::vector<std::vector<double>> g(q);
  for (int a = 0; a < q; ++a) g[a] = as<std::vector<double>>(grids[a]);
  int n1 = g[0].size(), n2 = q == 2 ? g[1].size() : 1;
  std::vector<double> diff((size_t)(n1 + 1) * (n2 + 1), 0);
  auto range = [](const std::vector<double>& v, double lo, double hi, int& a, int& b) {
    a = (lo <= 0 ? std::lower_bound(v.begin(), v.end(), lo) : std::upper_bound(v.begin(), v.end(), lo)) - v.begin();
    b = std::upper_bound(v.begin(), v.end(), hi) - v.begin();
  };
  struct Item { int j; double lo[2], hi[2]; };
  for (int t = 0; t < w.size(); ++t) {
    if (t % 16 == 0) checkUserInterrupt();
    std::vector<Item> stack{{start[t], {0, 0}, {1, 1}}};
    while (!stack.empty()) {
      Item it = stack.back(); stack.pop_back(); int j = it.j;
      if (left[j] >= 0) {
        Item l = it, r = it; l.j = left[j]; r.j = right[j];
        for (int a = 0; a < q; ++a) if (dim[j] == dims[a]) { l.hi[a] = cut[j]; r.lo[a] = cut[j]; }
        stack.push_back(l); stack.push_back(r); continue;
      }
      double val = w[t] * std::exp(lm[j]);
      int a1, b1, a2 = 0, b2 = 1;
      range(g[0], it.lo[0], it.hi[0], a1, b1); val /= it.hi[0] - it.lo[0];
      if (q == 2) { range(g[1], it.lo[1], it.hi[1], a2, b2); val /= it.hi[1] - it.lo[1]; }
      if (a1 >= b1 || a2 >= b2) continue;
      diff[(size_t)a1 * (n2 + 1) + a2] += val; diff[(size_t)b1 * (n2 + 1) + a2] -= val;
      diff[(size_t)a1 * (n2 + 1) + b2] -= val; diff[(size_t)b1 * (n2 + 1) + b2] += val;
    }
  }
  for (int a = 0; a <= n1; ++a) for (int b = 1; b <= n2; ++b) diff[(size_t)a * (n2 + 1) + b] += diff[(size_t)a * (n2 + 1) + b - 1];
  for (int a = 1; a <= n1; ++a) for (int b = 0; b <= n2; ++b) diff[(size_t)a * (n2 + 1) + b] += diff[(size_t)(a - 1) * (n2 + 1) + b];
  NumericVector out(n1 * n2);
  for (int a = 0; a < n1; ++a) for (int b = 0; b < n2; ++b)
    out[a + n1 * b] = std::max(0., diff[(size_t)a * (n2 + 1) + b]);
  return out;
}

// Internal verification entry: a fixed tree, independently enumerated in R tests.
// [[Rcpp::export]]
List ptd_fixed_cpp(IntegerMatrix children, IntegerVector counts, IntegerVector location,
                   IntegerVector depth, List prior) {
  using namespace ptd;
  int n = children.nrow();
  if (counts.size() != n || children.ncol() != 2 || location.size() != n || depth.size() != n)
    stop("Invalid fixed tree dimensions.");
  Prior p(prior, counts[0]); Tree t; t.nodes.resize(n);
  for (int j = 0; j < n; ++j) {
    Node& a = t.nodes[j]; a.left = children(j, 0) - 1; a.right = children(j, 1) - 1;
    a.depth = depth[j]; a.n = counts[j];
    if (a.depth < 0 || a.depth > p.depth) stop("Invalid fixed tree depth.");
    if (a.left >= 0) {
      if (a.left <= j || a.right <= j || a.left >= n || a.right >= n) stop("Invalid fixed tree child.");
      a.loc = location[j] - 1;
      if (a.loc < 0 || a.loc >= p.grid) stop("Invalid fixed tree location.");
      a.fraction = p.r(a.loc);
      t.nodes[a.left].parent = t.nodes[a.right].parent = j;
      emission(a.loc, counts[a.left], a.n, p, lgamma_n(a.n, p), a.emission);
    }
  }
  for (int j = n - 1; j >= 0; --j) refresh(t, j, p);
  return List::create(_["log_evidence"] = evidence(t, p), _["tree"] = table(t, p), _["lgamma_cap"] = p.cap);
}
