// Prefix index and direct sparse enumeration adapted from sparse_prefix_lbfp.cpp
// by Aurelien Nicosia (2026), MIT. See inst/COPYRIGHTS.
// The stencil is collapsed before traversal and supports ASH and GLBFP.
#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <cmath>
#include <functional>

struct GLBFPVectorHash {
  std::size_t operator()(const std::vector<int>& key) const noexcept {
    std::size_t hash = key.size();
    for (int value : key) hash ^= static_cast<std::size_t>(value) +
      0x9e3779b9U + (hash << 6) + (hash >> 2);
    return hash;
  }
};

struct GLBFPIndex {
  int d;
  std::vector<int> order;
  std::vector<std::unordered_set<std::string> > prefixes;
  std::unordered_map<std::string, int> counts;
  std::unordered_map<std::vector<int>, int, GLBFPVectorHash> direct_counts;
  GLBFPIndex(Rcpp::IntegerMatrix cells, Rcpp::IntegerVector count,
             Rcpp::IntegerVector perm, bool prefix) : d(cells.ncol()) {
    for (int j = 0; j < d; ++j) order.push_back(perm[j] - 1);
    if (prefix) prefixes.resize(d);
    else direct_counts.reserve(cells.nrow() * 2U + 1U);
    for (int r = 0; r < cells.nrow(); ++r) {
      if (!prefix) {
        std::vector<int> coordinates(d);
        for (int j = 0; j < d; ++j) coordinates[j] = cells(r,j);
        direct_counts[coordinates] = count[r];
        continue;
      }
      std::string key;
      for (int depth = 0; depth < d; ++depth) {
        if (depth) key += "|";
        key += std::to_string(cells(r, order[depth]));
        if (prefix) prefixes[depth].insert(key);
      }
      counts[key] = count[r];
    }
  }
};

// [[Rcpp::export]]
SEXP glbfp_index_cpp(Rcpp::IntegerMatrix cells, Rcpp::IntegerVector counts,
                     Rcpp::IntegerVector order, bool prefix) {
  return Rcpp::XPtr<GLBFPIndex>(new GLBFPIndex(cells, counts, order, prefix), true);
}

// A serialized external pointer must be rebuilt from the retained R counts.
// [[Rcpp::export]]
bool glbfp_index_valid_cpp(SEXP pointer) {
  return TYPEOF(pointer) == EXTPTRSXP && R_ExternalPtrAddr(pointer) != NULL;
}

// [[Rcpp::export]]
Rcpp::List glbfp_eval_cpp(SEXP pointer, Rcpp::IntegerMatrix base,
                         Rcpp::NumericMatrix u, Rcpp::IntegerVector m,
                         int n, double volume, bool ash,
                         Rcpp::IntegerMatrix self_keys) {
  Rcpp::XPtr<GLBFPIndex> p(pointer);
  const int nq = base.nrow(), d = p->d;
  Rcpp::NumericVector density(nq), sum_square(nq), self_weight(nq);
  Rcpp::NumericVector visited(nq), nodes(nq);
  // Candidate offsets and storage are reused over all queries.
  std::vector<std::vector<int> > offsets(d);
  std::vector<std::vector<double> > weights(d), lefts(d), rights(d);
  std::vector<int> query_key(d);
  for (int j = 0; j < d; ++j) {
    const int high = ash ? m[j] - 1 : m[j];
    for (int ell = 1 - m[j]; ell <= high; ++ell) {
      offsets[j].push_back(ell);
      lefts[j].push_back(std::max(0.0, 1.0 - std::abs(double(ell)) / m[j]));
      rights[j].push_back(std::max(0.0, 1.0 - std::abs(double(ell - 1)) / m[j]));
      weights[j].push_back(0);
    }
  }
  for (int i = 0; i < nq; ++i) {
    if (i % 256 == 0) Rcpp::checkUserInterrupt();
    for (int j = 0; j < d; ++j) {
      for (std::size_t r = 0; r < offsets[j].size(); ++r) {
        weights[j][r] = ash ? lefts[j][r] :
          (1 - u(i,j)) * lefts[j][r] + u(i,j) * rights[j][r];
      }
    }
    double acc = 0, square = 0, vi = 0, ni = 0;
    std::function<void(int, const std::string&, double)> rec;
    rec = [&](int depth, const std::string& key, double weight) {
      const int j = p->order[depth];
      for (std::size_t r = 0; r < offsets[j].size(); ++r) {
        if (weights[j][r] == 0) continue;
        const int coord = base(i,j) + offsets[j][r];
        const double w = weight * weights[j][r];
        const std::string next = key + (depth ? "|" : "") + std::to_string(coord);
        ++ni;
        if (p->prefixes[depth].find(next) == p->prefixes[depth].end()) continue;
        if (depth == d - 1) {
          auto found = p->counts.find(next);
          if (found != p->counts.end()) {
            ++vi;
            acc += w * found->second;
            square += w * w * found->second;
          }
        } else rec(depth + 1, next, w);
      }
    };
    if (p->prefixes.empty()) {
      std::function<void(int, double)> direct;
      direct = [&](int depth, double weight) {
        const int j = p->order[depth];
        for (std::size_t r = 0; r < offsets[j].size(); ++r) {
          if (weights[j][r] == 0) continue;
          query_key[j] = base(i,j) + offsets[j][r];
          const double w = weight * weights[j][r];
          ++ni;
          if (depth == d - 1) {
            auto found = p->direct_counts.find(query_key);
            if (found != p->direct_counts.end()) {
              ++vi;
              acc += w * found->second;
              square += w * w * found->second;
            }
          } else direct(depth + 1, w);
        }
      };
      direct(0, 1.0);
    } else rec(0, "", 1.0);
    density[i] = acc / (n * volume);
    sum_square[i] = square;
    visited[i] = vi;
    nodes[i] = ni;
    if (self_keys.nrow() == nq) {
      double w = 1;
      for (int j = 0; j < d; ++j) {
        const int ell = self_keys(i,j) - base(i,j);
        const double left = std::max(0.0, 1.0 - std::abs(double(ell)) / m[j]);
        w *= ash ? left : (1 - u(i,j)) * left + u(i,j) *
          std::max(0.0, 1.0 - std::abs(double(ell - 1)) / m[j]);
      }
      self_weight[i] = w;
    }
  }
  return Rcpp::List::create(Rcpp::_["densities"] = density,
    Rcpp::_["sum_square"] = sum_square, Rcpp::_["self_weight"] = self_weight,
    Rcpp::_["visited"] = visited, Rcpp::_["prefix_nodes"] = nodes);
}
