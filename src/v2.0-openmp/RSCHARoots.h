#pragma once
#include "RSCHACompatibility.h"
#include "Legendre.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <omp.h>
#include <vector>

namespace RSCHARoots {

struct Value {
  double value = 0.0;
  bool valid = false;
};

inline double predictor(int order, int index, long double theta) {
  return static_cast<double>(
      Legendre::detail::root_predictor(order, index, theta));
}

inline Value scaled(int sign, long double log_abs) {
  if (sign == 0) return {0.0, true};
  if (!std::isfinite(log_abs)) return {};
  return {static_cast<double>(sign) *
              std::exp(static_cast<double>(
                  RSCHACompat::clamp(log_abs, -700.0L, 700.0L))),
          true};
}

inline Value boundary(int order, double degree, double theta) {
  using namespace Legendre::detail;
  EvalResult result;
  const real nu = static_cast<real>(degree);
  const real angle = static_cast<real>(theta);
  if (degree >= 512.0 &&
      (evaluate_pole_residual_fast(nu, order, angle, result) ||
       evaluate_endpoint_p_fast(nu, order, angle, result) ||
       evaluate_near_order_recurrence(nu, order, angle, result) ||
       (degree < 10000.0 && angle <= pi/6.0L &&
        evaluate_order_downward(nu, order, angle, result)) ||
       wkb_evaluate_fast(nu, order, angle, result) ||
       wkb_evaluate_legacy(nu, order, angle, result) ||
       evaluate_canonical_legacy(nu, order, angle, result)) &&
      result.status == status_ok) {
    finalize_metadata(result);
    const Value value = scaled(result.sign, result.log_abs_value);
    if (value.valid)
      return value;
  }
  if (degree >= 1000.0 &&
      evaluate_global_airy_legacy(nu, order, angle, result, true)) {
    finalize_metadata(result);
    const Value value = scaled(result.sign, result.log_abs_value);
    if (value.valid)
      return value;
  }
  double p = 0.0, d1 = 0.0, d2 = 0.0;
  if (CRTable::legendre(order, degree, theta, p, d1, d2) &&
      std::isfinite(p) && p != 0.0)
    return {p, true};
  const auto direct = Legendre::EvaluateLegacySchmidt(order, degree, theta);
  return direct.status == status_ok
             ? scaled(direct.sign, direct.log_abs_value)
             : Value{};
}

inline Value direct_boundary(int order, double degree, double theta) {
  const auto value = Legendre::EvaluateLegacySchmidt(order, degree, theta);
  return value.status == Legendre::detail::status_ok
             ? scaled(value.sign, value.log_abs_value)
             : Value{};
}

inline double direct_root(int order, int index, double theta) {
  const double gap = static_cast<double>(Legendre::detail::pi) / theta;
  const double estimate = predictor(order, index, theta);
  double lo = std::max(static_cast<double>(order), estimate - 0.48 * gap);
  double hi = estimate + 0.48 * gap;
  Value flo = direct_boundary(order, lo, theta);
  Value fhi = direct_boundary(order, hi, theta);
  for (int expand = 0; expand < 8 &&
       (!flo.valid || !fhi.valid ||
        std::signbit(flo.value) == std::signbit(fhi.value));
       ++expand) {
    lo = std::max(static_cast<double>(order), lo - 0.25 * gap);
    hi += 0.25 * gap;
    flo = direct_boundary(order, lo, theta);
    fhi = direct_boundary(order, hi, theta);
  }
  if (!flo.valid || !fhi.valid ||
      std::signbit(flo.value) == std::signbit(fhi.value))
    return std::numeric_limits<double>::quiet_NaN();
  for (int iteration = 0; iteration < 48; ++iteration) {
    if (hi - lo <= 2.0e-13 * std::max(1.0, 0.5 * (lo + hi)))
      break;
    const double mid = lo + 0.5 * (hi - lo);
    const Value value = direct_boundary(order, mid, theta);
    if (!value.valid)
      return std::numeric_limits<double>::quiet_NaN();
    if (value.value == 0.0) return mid;
    if (std::signbit(value.value) == std::signbit(flo.value)) {
      lo = mid;
      flo = value;
    } else {
      hi = mid;
      fhi = value;
    }
  }
  return lo + 0.5 * (hi - lo);
}

inline double anchor_root(int order, int index, double theta) {
  const double gap = static_cast<double>(Legendre::detail::pi) / theta;
  const double estimate = predictor(order, index, theta);
  double lo = std::max(static_cast<double>(order), estimate - 0.24 * gap);
  double hi = estimate + 0.24 * gap;
  Value flo = boundary(order, lo, theta);
  Value fhi = boundary(order, hi, theta);
  for (int expand = 0; expand < 8 &&
       (!flo.valid || !fhi.valid ||
        std::signbit(flo.value) == std::signbit(fhi.value));
       ++expand) {
    lo = std::max(static_cast<double>(order), lo - 0.12 * gap);
    hi += 0.12 * gap;
    flo = boundary(order, lo, theta);
    fhi = boundary(order, hi, theta);
  }
  if (!flo.valid || !fhi.valid ||
      std::signbit(flo.value) == std::signbit(fhi.value))
    return direct_root(order, index, theta);
  int previous_side = 0;
  for (int iteration = 0; iteration < 32; ++iteration) {
    if (hi - lo <= 2.0e-13 * std::max(1.0, 0.5 * (lo + hi)))
      break;
    double candidate =
        (lo * fhi.value - hi * flo.value) / (fhi.value - flo.value);
    if (!std::isfinite(candidate) || candidate <= lo || candidate >= hi ||
        iteration % 10 == 9)
      candidate = lo + 0.5 * (hi - lo);
    const Value value = boundary(order, candidate, theta);
    if (!value.valid)
      return direct_root(order, index, theta);
    if (value.value == 0.0) return candidate;
    const int side =
        std::signbit(value.value) == std::signbit(flo.value) ? -1 : 1;
    if (side < 0) {
      lo = candidate;
      flo = value;
      if (previous_side == side)
        fhi.value *= 0.5;
    } else {
      hi = candidate;
      fhi = value;
      if (previous_side == side)
        flo.value *= 0.5;
    }
    previous_side = side;
  }
  return lo + 0.5 * (hi - lo);
}

// Interpolation is a starting value.  Every returned root must have a
// public-evaluator sign bracket at the requested degree tolerance.
inline double certify_root(int order, int index, double theta, double estimate) {
  if (!std::isfinite(estimate)) return direct_root(order, index, theta);
  const double tolerance = 2.0e-13 * std::max(1.0, std::abs(estimate));
  const double max_radius = 0.24 * static_cast<double>(Legendre::detail::pi) / theta;
  double radius = std::max(tolerance, 8.0 * std::numeric_limits<double>::epsilon() * std::abs(estimate));
  double lo = 0.0, hi = 0.0;
  Value flo, fhi;
  for (int attempt = 0; attempt < 48; ++attempt) {
    lo = std::max(static_cast<double>(order), estimate - radius);
    hi = estimate + radius;
    flo = direct_boundary(order, lo, theta);
    fhi = direct_boundary(order, hi, theta);
    if (flo.valid && fhi.valid) {
      if (flo.value == 0.0) return lo;
      if (fhi.value == 0.0) return hi;
      if (std::signbit(flo.value) != std::signbit(fhi.value)) break;
    }
    if (radius >= max_radius) return direct_root(order, index, theta);
    radius = std::min(max_radius, radius * 4.0);
  }
  if (!flo.valid || !fhi.valid || std::signbit(flo.value) == std::signbit(fhi.value))
    return std::numeric_limits<double>::quiet_NaN();
  for (int iteration = 0; iteration < 64 && hi - lo > 2.0 * tolerance; ++iteration) {
    const double mid = lo + 0.5 * (hi - lo);
    const Value value = direct_boundary(order, mid, theta);
    if (!value.valid) return std::numeric_limits<double>::quiet_NaN();
    if (value.value == 0.0) return mid;
    if (std::signbit(value.value) == std::signbit(flo.value)) {
      lo = mid; flo = value;
    } else hi = mid;
  }
  return lo + 0.5 * (hi - lo);
}

inline void compute(int maximum, double theta_degrees,
                    std::vector<double> &roots, int first_order = 0,
                    int order_step = 1) {
  constexpr int initial = 32;
  constexpr int stride = 16;
  constexpr int interpolation_points = 8;
  const double theta = theta_degrees *
                       static_cast<double>(Legendre::detail::pi) / 180.0;
#pragma omp parallel for schedule(dynamic, 1)
  for (int order = first_order; order <= maximum; order += order_step) {
    const int count = maximum - order + 1;
    std::vector<int> anchors;
    for (int i = 0; i < std::min(initial, count); ++i)
      anchors.push_back(i);
    for (int i = initial; i < count; i += stride)
      anchors.push_back(i);
    if (anchors.back() != count - 1)
      anchors.push_back(count - 1);
    std::vector<double> corrections(anchors.size());
    for (size_t i = 0; i < anchors.size(); ++i) {
      const int index = anchors[i] + 1;
      corrections[i] = anchor_root(order, index, theta) -
                       predictor(order, index, theta);
    }
    for (int i = 0; i < count; ++i) {
      const auto upper = std::lower_bound(anchors.begin(), anchors.end(), i);
      const size_t center = upper == anchors.end()
                                ? anchors.size() - 1
                                : static_cast<size_t>(upper - anchors.begin());
      size_t first = center >= 4 ? center - 4 : 0;
      if (first + interpolation_points > anchors.size())
        first = anchors.size() > interpolation_points
                    ? anchors.size() - interpolation_points
                    : 0;
      const size_t used =
          std::min<size_t>(interpolation_points, anchors.size());
      const double x = 1.0 / static_cast<double>(i + 1);
      double correction = 0.0;
      for (size_t a = 0; a < used; ++a) {
        const double xa =
            1.0 / static_cast<double>(anchors[first + a] + 1);
        double basis = 1.0;
        for (size_t b = 0; b < used; ++b)
          if (a != b) {
            const double xb =
                1.0 / static_cast<double>(anchors[first + b] + 1);
            basis *= (x - xb) / (xa - xb);
          }
        correction += basis * corrections[first + a];
      }
      const int k = order + i;
      double root = certify_root(order, i + 1, theta,
          predictor(order, i + 1, theta) + correction);
      if (i > 0 && !(root > roots[(k - 1) * k / 2 + order]))
        root = direct_root(order, i + 1, theta);
      if (i > 0 && !(root > roots[(k - 1) * k / 2 + order]))
        root = std::numeric_limits<double>::quiet_NaN();
      roots[k * (k + 1) / 2 + order] = root;
    }
  }
}

}
