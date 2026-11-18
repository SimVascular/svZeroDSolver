// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file TangentCheck.h
 * @brief Temporary comparison of exact and complex-step tangents
 */
#ifndef SVZERODSOLVER_MODEL_TANGENTCHECK_HPP_
#define SVZERODSOLVER_MODEL_TANGENTCHECK_HPP_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <utility>

/// Largest relative difference per category, printed at exit when the
/// environment variable SVZERO_TANGENT_CHECK is set.
struct TangentCheck {
  bool on = std::getenv("SVZERO_TANGENT_CHECK") != nullptr;  ///< enabled
  std::map<std::string, std::pair<double, long>> worst;      ///< per category
  std::map<std::string, std::pair<double, double>> abs_;  ///< max |a - b|, |b|
  /// Record max|a - b| / max|b| over n entries.
  void add(const std::string& key, const double* a, const double* b, int n) {
    double diff = 0.0, scale = 0.0;
    for (int i = 0; i < n; i++) {
      diff = std::max(diff, std::abs(a[i] - b[i]));
      scale = std::max(scale, std::abs(b[i]));
    }
    auto& w = worst[key];
    w.first = std::max(w.first, scale > 0.0 ? diff / scale : diff);
    w.second++;
    auto& d = abs_[key];
    if (diff > d.first) d = {diff, scale};
  }
  ~TangentCheck() {
    for (const auto& kv : worst)
      std::fprintf(stderr,
                   "tangent check %-28s max rel diff %.3e (%ld evals), max abs "
                   "diff %.3e at |ref| %.3e\n",
                   kv.first.c_str(), kv.second.first, kv.second.second,
                   abs_[kv.first].first, abs_[kv.first].second);
  }
};

/// The process-wide check.
inline TangentCheck& tangent_check() {
  static TangentCheck tc;
  return tc;
}

#endif  // SVZERODSOLVER_MODEL_TANGENTCHECK_HPP_
