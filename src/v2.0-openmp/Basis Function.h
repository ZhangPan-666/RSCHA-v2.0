#pragma once
#include "Legendre.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <unordered_map>

namespace DegreeTrig {
constexpr double radians_per_degree =
    3.141592653589793238462643383279502884 / 180.0;
inline double sind(double degree) {
  return std::sin(degree * radians_per_degree);
}
inline double cosd(double degree) {
  return std::cos(degree * radians_per_degree);
}
}

namespace Mehler {

struct BasisRecord {
  double value = std::numeric_limits<double>::quiet_NaN();
  double dtheta = std::numeric_limits<double>::quiet_NaN();
  double d2theta = std::numeric_limits<double>::quiet_NaN();
  int status = -1;
  long double log_abs_value = -std::numeric_limits<long double>::infinity();
  long double log_abs_dtheta = -std::numeric_limits<long double>::infinity();
  long double log_abs_d2theta = -std::numeric_limits<long double>::infinity();
  int sign = 0;
  int dtheta_sign = 0;
  int d2theta_sign = 0;
  long double dtheta_over_value = std::numeric_limits<long double>::quiet_NaN();
  long double d2theta_over_value =
      std::numeric_limits<long double>::quiet_NaN();
  unsigned range_flags = 0;
};

namespace ConicalDetail {

using Real = long double;
constexpr Real PI = 3.141592653589793238462643383279502884L;
inline Real deg2rad(Real x) { return x * PI / 180.0L; }
inline Real sqr(Real x) { return x * x; }
struct PolarData {
  Real s = 0.0L;
  Real log_raw = 0.0L;
  Real ratio = 0.0L;
  Real ratio2 = 0.0L;
};

inline PolarData polar_data(int m, Real lambda, Real theta) {
  PolarData out;
  const Real z = sqr(std::sin(theta / 2.0L));
  Real term = 1.0L;
  Real term_lambda = 0.0L;
  Real f = 1.0L, fz = 0.0L, fzz = 0.0L;
  Real flambda = 0.0L, fzlambda = 0.0L;
  Real term_z = 0.0L, term_zz = 0.0L, term_zlambda = 0.0L;
  for (int k = 1; k <= 64; ++k) {
    const Real kr = static_cast<Real>(k);
    const Real den = kr * (static_cast<Real>(m) + kr);
    const Real factor =
        (static_cast<Real>(m + k) * static_cast<Real>(m + k - 1) + lambda) /
        den;
    const Real previous = term, previous_lambda = term_lambda;
    const Real previous_z = term_z, previous_zz = term_zz;
    term_zlambda = (term_zlambda*z + previous_lambda)*factor +
                   (previous_z*z + previous)/den;
    term_zz = (previous_zz*z + 2*previous_z)*factor;
    term_z = (previous_z*z + previous)*factor;
    term_lambda = previous_lambda*factor*z + previous*z/den;
    term = previous*factor*z;
    f += term;
    flambda += term_lambda;
    fz += term_z;
    fzlambda += term_zlambda;
    fzz += term_zz;
  }
  const Real st = std::sin(theta);
  const Real ct = std::cos(theta);
  const Real logarithmic_fz = fz / f;
  const Real ftheta_over_f = 0.5L * st * logarithmic_fz;
  out.s = 0.5L * st * st * (fzlambda * f - fz * flambda) / sqr(f);
  out.log_raw =
      static_cast<Real>(m) * std::log(st) + std::log(f);
  out.ratio = static_cast<Real>(m) * ct / st + ftheta_over_f;
  out.ratio2 = (m > 1 ? static_cast<Real>(m * (m - 1)) * sqr(ct / st) : 0.0L) -
               static_cast<Real>(m) +
               2.0L * static_cast<Real>(m) * ct / st * ftheta_over_f +
               (0.25L * st * st * fzz + 0.5L * ct * fz) / f;
  return out;
}

inline PolarData positive_data(int m, Real lambda, Real theta) {
  PolarData out;
  const Real z = sqr(std::sin(theta / 2.0L));
  Real term = 1.0L, term_lambda = 0.0L;
  Real f = 1.0L, fz = 0.0L, fzz = 0.0L;
  Real flambda = 0.0L, fzlambda = 0.0L;
  Real term_z = 0.0L, term_zz = 0.0L, term_zlambda = 0.0L;
  bool converged = false;
  for (int k = 1; k <= 256; ++k) {
    const Real kr = k;
    const Real den = kr * (m + kr);
    const Real factor = (kr * (kr - 1.0L) + lambda) / den;
    const Real previous = term, previous_lambda = term_lambda;
    const Real previous_z = term_z, previous_zz = term_zz;
    term_zlambda = (term_zlambda*z + previous_lambda)*factor +
                   (previous_z*z + previous)/den;
    term_zz = (previous_zz*z + 2*previous_z)*factor;
    term_z = (previous_z*z + previous)*factor;
    term_lambda = previous_lambda*factor*z + previous*z/den;
    term = previous*factor*z;
    f += term;
    flambda += term_lambda;
    fz += term_z;
    fzlambda += term_zlambda;
    fzz += term_zz;
    const Real bound = z * (1.0L + lambda / ((kr + 1.0L) * (m + kr + 1.0L)));
    if (bound < 1.0L &&
        term * sqr(kr + 1.0L) < 1.0e-18L * f * std::pow(1.0L-bound,3) &&
        term_lambda * (kr + 1.0L) <= 1.0e-18L * flambda * sqr(1.0L-bound)) {
      converged = true;
      break;
    }
  }
  if (!converged) {
    out.log_raw = std::numeric_limits<Real>::quiet_NaN();
    return out;
  }
  const Real st = std::sin(theta), ct = std::cos(theta);
  const Real hz = fz / f;
  out.s = 0.5L * st * st * (fzlambda * f - fz * flambda) / sqr(f);
  out.log_raw = m * std::log(2.0L * std::tan(theta / 2.0L)) + std::log(f);
  out.ratio = m / st + 0.5L * st * hz;
  out.ratio2 = (m > 1 ? Real(m)*(m-1)/(st*st) : 0.0L) +
               m/(1.0L+ct) + (m+0.5L*ct)*hz + 0.25L*st*st*fzz/f;
  return out;
}

struct GammaData {
  std::complex<Real> logarithm;
  std::complex<Real> digamma;
};

inline GammaData gamma_data(std::complex<Real> z) {
  GammaData out{};
  while (std::abs(z) < 16.0L) {
    out.logarithm -= std::log(z);
    out.digamma -= 1.0L / z;
    z += 1.0L;
  }
  const auto inv = 1.0L / z;
  const auto inv2 = inv * inv;
  auto power = inv;
  out.logarithm += (z-0.5L)*std::log(z)-z+0.5L*std::log(2.0L*PI);
  out.digamma += std::log(z)-0.5L*inv;
  constexpr Real bernoulli[] = {1.0L/6.0L,-1.0L/30.0L,1.0L/42.0L,
      -1.0L/30.0L,5.0L/66.0L,-691.0L/2730.0L,7.0L/6.0L,-3617.0L/510.0L};
  for (int k=1;k<=8;++k) {
    out.logarithm += bernoulli[k-1]*power/(Real(2*k)*(2*k-1));
    out.digamma -= bernoulli[k-1]*power*inv/Real(2*k);
    power *= inv2;
  }
  return out;
}

inline PolarData equator_data(int m, Real tau, Real theta) {
  const auto a=gamma_data({0.5L*m+0.75L,0.5L*tau});
  const auto b=gamma_data({0.5L*m+0.25L,0.5L*tau});
  PolarData out;
  const Real r=2.0L*std::exp(2.0L*(a.logarithm.real()-b.logarithm.real()));
  const Real delta=theta-1.5707963267948965579989817342720925807952880859375L-
                   6.1232339957367660359e-17L;
  out.s=r*(b.digamma.imag()-a.digamma.imag())/(2.0L*tau);
  out.s+=(1.0L-2.0L*r*out.s)*delta;
  out.log_raw=std::lgamma(Real(m)+1.0L)+0.5L*std::log(PI)-
              2.0L*a.logarithm.real()+r*delta;
  const Real lambda=tau*tau+0.25L;
  out.ratio=r+(Real(m)*m+lambda-r*r)*delta;
  const Real st=std::sin(theta),ct=std::cos(theta);
  out.ratio2=lambda+Real(m)*m/(st*st)-(ct/st)*out.ratio;
  return out;
}

inline Real safe_exp(Real x) {
  if (std::isnan(x))
    return x;
  const Real hi = std::log(static_cast<Real>(std::numeric_limits<double>::max()));
  const Real lo = std::log(static_cast<Real>(std::numeric_limits<double>::denorm_min()));
  if (x > hi)
    return std::numeric_limits<Real>::infinity();
  if (x < lo)
    return 0.0L;
  if (x >= std::log(static_cast<Real>(std::numeric_limits<double>::min())))
    return std::exp(x);
  const std::uint64_t bits = static_cast<std::uint64_t>(std::floor(std::exp(x-lo)+0.5L));
  static_assert(std::numeric_limits<double>::is_iec559 && sizeof(double) == sizeof(bits),
                "Subnormal materialization requires IEEE binary64.");
  double value;
  std::memcpy(&value, &bits, sizeof(value));
  return static_cast<Real>(value);
}

inline void finalize_record(BasisRecord &out) noexcept {
  if (out.status != 0)
    return;
  if (out.value != 0.0) {
    out.sign = out.value > 0.0 ? 1 : -1;
    if (!std::isfinite(out.log_abs_value))
      out.log_abs_value = std::log(std::fabs(static_cast<Real>(out.value)));
    if (!std::isfinite(out.dtheta_over_value))
      out.dtheta_over_value = static_cast<Real>(out.dtheta) / out.value;
    if (!std::isfinite(out.d2theta_over_value))
      out.d2theta_over_value = static_cast<Real>(out.d2theta) / out.value;
  }
  if (out.dtheta != 0.0) {
    out.dtheta_sign = out.dtheta > 0.0 ? 1 : -1;
    if (!std::isfinite(out.log_abs_dtheta))
      out.log_abs_dtheta = std::log(std::fabs(static_cast<Real>(out.dtheta)));
  } else if (std::isfinite(out.log_abs_value) &&
             std::isfinite(out.dtheta_over_value) &&
             out.dtheta_over_value != 0.0L) {
    out.dtheta_sign = out.sign * (out.dtheta_over_value > 0.0L ? 1 : -1);
    out.log_abs_dtheta =
        out.log_abs_value + std::log(std::fabs(out.dtheta_over_value));
  }
  if (out.d2theta != 0.0) {
    out.d2theta_sign = out.d2theta > 0.0 ? 1 : -1;
    if (!std::isfinite(out.log_abs_d2theta))
      out.log_abs_d2theta = std::log(std::fabs(static_cast<Real>(out.d2theta)));
  } else if (std::isfinite(out.log_abs_value) &&
             std::isfinite(out.d2theta_over_value) &&
             out.d2theta_over_value != 0.0L) {
    out.d2theta_sign = out.sign * (out.d2theta_over_value > 0.0L ? 1 : -1);
    out.log_abs_d2theta =
        out.log_abs_value + std::log(std::fabs(out.d2theta_over_value));
  }
  Legendre::detail::classify_log_range(
      out.log_abs_value, Legendre::detail::range_value_underflow,
      Legendre::detail::range_value_overflow, out.range_flags);
  Legendre::detail::classify_log_range(
      out.log_abs_dtheta, Legendre::detail::range_dtheta_underflow,
      Legendre::detail::range_dtheta_overflow, out.range_flags);
  Legendre::detail::classify_log_range(
      out.log_abs_d2theta, Legendre::detail::range_d2theta_underflow,
      Legendre::detail::range_d2theta_overflow, out.range_flags);
}

inline BasisRecord evaluate_table_impl(int m, double tau_d,
                                       double theta_deg_d,
                                       double theta0_deg_d) {
  BasisRecord out;
#ifndef RSCHA_CR_TABLE_BUILD
  if (m < 0 || m > 10000 || !std::isfinite(tau_d) || std::fabs(tau_d) > 10000.0 ||
      !(theta_deg_d >= 0.0 && theta_deg_d <= theta0_deg_d) ||
      !(theta0_deg_d > 0.0 && theta0_deg_d <= 90.0))
    return out;
  const double theta = static_cast<double>(deg2rad(theta_deg_d));
  const double theta0 = static_cast<double>(deg2rad(theta0_deg_d));
  double log_raw0 = 0.0, sensitivity0 = 0.0, ratio0 = 0.0;
  const Real tau = std::fabs(static_cast<Real>(tau_d));
  const Real lambda = tau * tau + 0.25L;
  const bool positive = lambda <= 16.0L * (m + 1.0L);
  if (positive) {
    const PolarData endpoint0 = positive_data(m, lambda, theta0);
    log_raw0 = static_cast<double>(endpoint0.log_raw);
    sensitivity0 = static_cast<double>(endpoint0.s);
  } else if (theta0_deg_d == 90.0) {
    const PolarData endpoint0 = equator_data(m, tau, theta0);
    if (!CRTable::mehler(m, tau_d, theta0, log_raw0, sensitivity0, ratio0))
      return out;
    // Table-to-table ratios must cancel the common atlas amplitude offset.
    if (theta_deg_d <= 0.1) log_raw0 = static_cast<double>(endpoint0.log_raw);
    sensitivity0 = static_cast<double>(endpoint0.s);
  } else if (theta0_deg_d <= 0.1) {
    const PolarData endpoint0 = polar_data(m, lambda, theta0);
    log_raw0 = static_cast<double>(endpoint0.log_raw);
    sensitivity0 = static_cast<double>(endpoint0.s);
  } else if (!CRTable::mehler(m, tau_d, theta0, log_raw0, sensitivity0,
                              ratio0)) {
    return out;
  }
  if (!(sensitivity0 > 0.0) || !std::isfinite(log_raw0))
    return out;
  if (theta_deg_d == 0.0 || theta == 0.0) {
    // Every unassigned axis limit is exactly zero, not an invalid field.
    out.value = out.dtheta = out.d2theta = 0.0;
    const Real log_axis = -static_cast<Real>(log_raw0) -
                          0.5L * std::log(static_cast<Real>(sensitivity0));
    const Real amplitude = safe_exp(log_axis);
    if (m == 0) {
      out.value = static_cast<double>(amplitude);
      out.log_abs_d2theta = log_axis + std::log(0.5L * lambda);
      out.d2theta = static_cast<double>(safe_exp(out.log_abs_d2theta));
      out.d2theta_sign = 1;
      out.log_abs_value = log_axis;
      out.sign = 1;
      out.dtheta_over_value = 0.0L;
      out.d2theta_over_value = 0.5L * lambda;
    } else if (m == 1) {
      out.dtheta = static_cast<double>(amplitude);
      out.log_abs_dtheta = log_axis;
      out.dtheta_sign = 1;
    } else if (m == 2) {
      out.log_abs_d2theta = std::log(2.0L) + log_axis;
      out.d2theta = static_cast<double>(safe_exp(out.log_abs_d2theta));
      out.d2theta_sign = 1;
    }
    out.status = 0;
    finalize_record(out);
    return out;
  }
  double log_raw = 0.0, sensitivity = 0.0, ratio = 0.0;
  Real polar_ratio2 = 0.0L;
  if (positive) {
    const PolarData endpoint = positive_data(m, lambda, theta);
    log_raw = static_cast<double>(endpoint.log_raw);
    ratio = static_cast<double>(endpoint.ratio);
    polar_ratio2 = endpoint.ratio2;
  } else if (theta_deg_d == 90.0) {
    const PolarData endpoint = equator_data(m, tau, theta);
    log_raw = log_raw0;
    ratio = static_cast<double>(endpoint.ratio);
    polar_ratio2 = endpoint.ratio2;
  } else if (theta_deg_d <= 0.1) {
    const PolarData endpoint = polar_data(m, lambda, theta);
    log_raw = static_cast<double>(endpoint.log_raw);
    sensitivity = static_cast<double>(endpoint.s);
    ratio = static_cast<double>(endpoint.ratio);
    polar_ratio2 = endpoint.ratio2;
  } else if (!CRTable::mehler(m, tau_d, theta, log_raw, sensitivity, ratio)) {
    return out;
  }
  const Real logf = static_cast<Real>(log_raw - log_raw0) -
                    0.5L * std::log(static_cast<Real>(sensitivity0));
  const Real f = safe_exp(logf);
  const Real st = std::sin(static_cast<Real>(theta));
  const Real ct = std::cos(static_cast<Real>(theta));
  const Real r = static_cast<Real>(ratio);
  const auto derivative = [logf, f](Real factor) {
    if (factor == 0.0L)
      return Real(0.0L);
    return f >= std::numeric_limits<double>::min() ? factor * f
        : std::copysign(safe_exp(logf + std::log(std::fabs(factor))), factor);
  };
  const Real d = derivative(r);
  const Real ratio2 = positive || theta_deg_d <= 0.1
                          ? polar_ratio2
                          : -(ct / st) * r + lambda +
                                sqr(static_cast<Real>(m)) / (st * st);
  // On Windows, long double is binary64: a nonzero subnormal sin^2(theta)
  // can still overflow m*(m-1)/sin^2(theta). Materialize the polar second
  // derivative from its own logarithm before that ratio is used. In this
  // domain (m, tau <= 1e4) the omitted relative correction is below 1e-292.
  const bool tiny_second = m > 1 && theta < 1.0e-150;
  const Real log_second = tiny_second ? logf - 2.0L*std::log(theta) +
      std::log(Real(m)*(m-1)) : (ratio2 == 0.0L ?
      -std::numeric_limits<Real>::infinity() : logf + std::log(std::fabs(ratio2)));
  const Real dd = tiny_second ? safe_exp(log_second) : derivative(ratio2);
  out.value = static_cast<double>(f);
  out.dtheta = static_cast<double>(d);
  out.d2theta = static_cast<double>(dd);
  out.log_abs_value = logf;
  out.sign = 1;
  out.dtheta_over_value = r;
  out.log_abs_dtheta = r == 0.0L ? -std::numeric_limits<Real>::infinity()
      : logf + std::log(std::fabs(r));
  out.dtheta_sign = (r > 0.0L) - (r < 0.0L);
  out.d2theta_over_value = ratio2;
  out.log_abs_d2theta = log_second;
  out.d2theta_sign = tiny_second ? 1 : ((ratio2 > 0.0L) - (ratio2 < 0.0L));
  out.status = std::isfinite(out.value) && std::isfinite(out.dtheta) &&
                       std::isfinite(out.d2theta)
                   ? 0
                   : -1;
  finalize_record(out);
#else
  (void)m;
  (void)tau_d;
  (void)theta_deg_d;
  (void)theta0_deg_d;
#endif
  return out;
}

}

inline BasisRecord Evaluate(double m, double tau, double theta, double theta0) {
  const int order = std::isfinite(m) && m >= 0.0 && m <= 10000.0 &&
      std::fabs(m - std::round(m)) <=
          8.0 * std::numeric_limits<double>::epsilon() * (1.0 + m)
      ? static_cast<int>(std::round(m)) : -1;
  return ConicalDetail::evaluate_table_impl(order, tau, theta, theta0);
}

inline void Mehler_1rank(double m, double tau, double theta, double theta0,
                         double &value, double &dtheta) {
  const BasisRecord result = Evaluate(m, tau, theta, theta0);
  value = result.status == 0 ? result.value : std::numeric_limits<double>::quiet_NaN();
  dtheta = result.status == 0 ? result.dtheta : std::numeric_limits<double>::quiet_NaN();
}

inline void Mehler_2rank(double m, double tau, double theta, double theta0,
                         double &value, double &dtheta, double &d2theta) {
  const BasisRecord result = Evaluate(m, tau, theta, theta0);
  value = result.status == 0 ? result.value : std::numeric_limits<double>::quiet_NaN();
  dtheta = result.status == 0 ? result.dtheta : std::numeric_limits<double>::quiet_NaN();
  d2theta = result.status == 0 ? result.d2theta : std::numeric_limits<double>::quiet_NaN();
}

}

namespace ConstantDetail {

using Real = long double;
constexpr Real PI = 3.141592653589793238462643383279502884L;
constexpr int STATUS_OK = 0;
constexpr int STATUS_INVALID_ORDER = 1;
constexpr int STATUS_INVALID_ANGLE = 2;
constexpr int STATUS_SERIES_FAILURE = 3;
constexpr int STATUS_NONFINITE = 4;

inline Real deg2rad(Real degree) noexcept { return degree * PI / 180.0L; }

inline Real sqr(Real x) noexcept { return x * x; }

inline Real tan_half_degree(Real degree) noexcept {
  if (degree <= 90.0L)
    return std::tan(deg2rad(0.5L * degree));
  const Real complement = 180.0L - degree;
  return 1.0L / std::tan(deg2rad(0.5L * complement));
}

inline Real sin_degree(Real degree) noexcept {
  if (degree <= 90.0L)
    return std::sin(deg2rad(degree));
  return std::sin(deg2rad(180.0L - degree));
}

inline Real cos_degree(Real degree) noexcept {
  if (degree <= 90.0L)
    return std::cos(deg2rad(degree));
  return -std::cos(deg2rad(180.0L - degree));
}

inline Real sin_half_degree(Real degree) noexcept {
  if (degree <= 90.0L)
    return std::sin(deg2rad(0.5L * degree));
  return std::cos(deg2rad(0.5L * (180.0L - degree)));
}

struct BasisResult {
  double value = std::numeric_limits<double>::quiet_NaN();
  double dtheta = std::numeric_limits<double>::quiet_NaN();
  double d2theta = std::numeric_limits<double>::quiet_NaN();
  int status = STATUS_OK;
  long double log_abs_value = -std::numeric_limits<long double>::infinity();
  long double log_abs_dtheta = -std::numeric_limits<long double>::infinity();
  long double log_abs_d2theta = -std::numeric_limits<long double>::infinity();
  int sign = 0;
  int dtheta_sign = 0;
  int d2theta_sign = 0;
  long double dtheta_over_value = std::numeric_limits<long double>::quiet_NaN();
  long double d2theta_over_value =
      std::numeric_limits<long double>::quiet_NaN();
  unsigned range_flags = 0;
};

struct NormData {
  Real u0 = 0.0L;
  Real log_norm_factor = 0.0L;
  int status = STATUS_OK;
};

struct CacheKey {
  int m = 0;
  std::uint64_t theta0_bits = 0;
  bool operator==(const CacheKey &other) const noexcept {
    return m == other.m && theta0_bits == other.theta0_bits;
  }
};

struct CacheKeyHash {
  std::size_t operator()(const CacheKey &key) const noexcept {
    const std::size_t h1 = std::hash<int>{}(key.m);
    const std::size_t h2 = std::hash<std::uint64_t>{}(key.theta0_bits);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6U) + (h1 >> 2U));
  }
};

inline std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

inline bool integer_order(double m, int &order) noexcept {
  if (!std::isfinite(m) || m < 0.0)
    return false;
  const double rounded = std::nearbyint(m);
  const double scale = std::max(1.0, std::fabs(m));
  if (std::fabs(m - rounded) >
      32.0 * std::numeric_limits<double>::epsilon() * scale)
    return false;
  if (rounded > static_cast<double>(std::numeric_limits<int>::max()))
    return false;
  order = static_cast<int>(rounded);
  return true;
}

inline double from_log(Real log_abs) noexcept {
  return static_cast<double>(Mehler::ConicalDetail::safe_exp(log_abs));
}

inline void finalize_result(BasisResult &out) noexcept {
  if (out.status != STATUS_OK)
    return;
  if (out.value != 0.0) {
    out.sign = out.value > 0.0 ? 1 : -1;
    if (!std::isfinite(out.log_abs_value))
      out.log_abs_value = std::log(std::fabs(static_cast<Real>(out.value)));
  }
  if (out.dtheta != 0.0) {
    out.dtheta_sign = out.dtheta > 0.0 ? 1 : -1;
    if (!std::isfinite(out.log_abs_dtheta))
      out.log_abs_dtheta = std::log(std::fabs(static_cast<Real>(out.dtheta)));
  } else if (std::isfinite(out.log_abs_value) &&
             std::isfinite(out.dtheta_over_value) &&
             out.dtheta_over_value != 0.0L) {
    out.dtheta_sign = out.sign * (out.dtheta_over_value > 0.0L ? 1 : -1);
    out.log_abs_dtheta =
        out.log_abs_value + std::log(std::fabs(out.dtheta_over_value));
  }
  if (out.d2theta != 0.0) {
    out.d2theta_sign = out.d2theta > 0.0 ? 1 : -1;
    if (!std::isfinite(out.log_abs_d2theta))
      out.log_abs_d2theta = std::log(std::fabs(static_cast<Real>(out.d2theta)));
  } else if (std::isfinite(out.log_abs_value) &&
             std::isfinite(out.d2theta_over_value) &&
             out.d2theta_over_value != 0.0L) {
    out.d2theta_sign = out.sign * (out.d2theta_over_value > 0.0L ? 1 : -1);
    out.log_abs_d2theta =
        out.log_abs_value + std::log(std::fabs(out.d2theta_over_value));
  }
  Legendre::detail::classify_log_range(
      out.log_abs_value, Legendre::detail::range_value_underflow,
      Legendre::detail::range_value_overflow, out.range_flags);
  Legendre::detail::classify_log_range(
      out.log_abs_dtheta, Legendre::detail::range_dtheta_underflow,
      Legendre::detail::range_dtheta_overflow, out.range_flags);
  Legendre::detail::classify_log_range(
      out.log_abs_d2theta, Legendre::detail::range_d2theta_underflow,
      Legendre::detail::range_d2theta_overflow, out.range_flags);
}

inline std::shared_ptr<const NormData> axis_norm_data(int m,
                                                      double theta0_degree) {
  auto data = std::make_shared<NormData>();
  data->u0 = tan_half_degree(static_cast<Real>(theta0_degree));
  const Real q = data->u0 * data->u0;
  if (!(q > 0.0L) || !std::isfinite(q)) {
    data->status = STATUS_INVALID_ANGLE;
    return data;
  }
  const Real a = q / (1.0L + q);
  const Real one_minus_a = 1.0L / (1.0L + q);
  Real term = 1.0L / static_cast<Real>(m + 1);
  Real sum = term;
  for (int k = 0; k < 3; ++k) {
    term *= static_cast<Real>(k + 2) * a /
            static_cast<Real>(m + k + 2);
    sum += term;
  }
  sum *= one_minus_a * one_minus_a;
  data->log_norm_factor =
      0.5L * (std::log(2.0L) + std::log(q) + std::log(sum));
  return data;
}

inline std::shared_ptr<const NormData> get_norm_data(int m,
                                                     double theta0_degree) {
  static std::unordered_map<CacheKey, std::shared_ptr<const NormData>,
                            CacheKeyHash>
      cache;
  static std::mutex cache_mutex;

  if (theta0_degree == 0.0)
    theta0_degree = 0.0;
  const CacheKey key{m, double_bits(theta0_degree)};
  {
    const std::lock_guard<std::mutex> lock(cache_mutex);
    const auto it = cache.find(key);
    if (it != cache.end())
      return it->second;
  }

  std::shared_ptr<const NormData> value;
  if (m == 0) {
    auto data = std::make_shared<NormData>();
    const Real theta0 = deg2rad(static_cast<Real>(theta0_degree));
    if (!(theta0 > 0.0L && theta0 < PI)) {
      data->status = STATUS_INVALID_ANGLE;
    } else {
      const Real cap_measure =
          2.0L * sqr(sin_half_degree(static_cast<Real>(theta0_degree)));
      data->log_norm_factor = 0.5L * std::log(cap_measure);
    }
    value = data;
  } else if (theta0_degree < 0.01) {
    value = axis_norm_data(m, theta0_degree);
  } else {
#ifndef RSCHA_CR_TABLE_BUILD
    double tabulated = 0.0;
    if (CRTable::constant_norm(m, static_cast<double>(deg2rad(theta0_degree)),
                               tabulated)) {
      auto data = std::make_shared<NormData>();
      data->log_norm_factor = static_cast<Real>(tabulated);
      data->u0 = tan_half_degree(static_cast<Real>(theta0_degree));
      value = data;
    }
#endif
    if (!value) {
      auto data = std::make_shared<NormData>();
      data->status = STATUS_NONFINITE;
      value = data;
    }
  }
  {
    const std::lock_guard<std::mutex> lock(cache_mutex);
    const auto result = cache.emplace(key, value);
    return result.first->second;
  }
}

inline BasisResult evaluate(int m, double theta_degree,
                            double theta0_degree) noexcept {
  BasisResult out;
  try {
    const auto norm = get_norm_data(m, theta0_degree);
    if (!norm || norm->status != STATUS_OK) {
      out.status = norm ? norm->status : STATUS_NONFINITE;
      return out;
    }

    Real theta_degrees = static_cast<Real>(theta_degree);
    if (!std::isfinite(theta_degrees) ||
        !(theta_degrees >= 0.0L && theta_degrees < 180.0L)) {
      out.status = STATUS_INVALID_ANGLE;
      return out;
    }
    const Real angle_tol =
        32.0L * std::numeric_limits<Real>::epsilon() *
        std::max<Real>(1.0L, std::fabs(static_cast<Real>(theta0_degree)));
    if (theta_degrees > static_cast<Real>(theta0_degree)) {
      if (theta_degrees - static_cast<Real>(theta0_degree) > angle_tol) {
        out.status = STATUS_INVALID_ANGLE;
        return out;
      }
      theta_degrees = static_cast<Real>(theta0_degree);
    }

    if (m == 0) {
      const Real log_value = -norm->log_norm_factor;
      out.value = from_log(log_value);
      out.dtheta = 0.0;
      out.d2theta = 0.0;
      out.log_abs_value = log_value;
      out.sign = 1;
      out.dtheta_over_value = 0.0L;
      out.d2theta_over_value = 0.0L;
      finalize_result(out);
      return out;
    }

    if (theta_degrees == 0.0L || theta_degrees * PI / 180.0L == 0.0L) {
      const Real log_A =
          -static_cast<Real>(m) * (std::log(2.0L) + std::log(norm->u0)) -
          norm->log_norm_factor;
      out.value = 0.0;
      out.dtheta = (m == 1) ? from_log(log_A) : 0.0;
      out.d2theta = (m == 2) ? from_log(std::log(2.0L) + log_A) : 0.0;
      if (m == 1) {
        out.log_abs_dtheta = log_A;
        out.dtheta_sign = 1;
      }
      if (m == 2) {
        out.log_abs_d2theta = std::log(2.0L) + log_A;
        out.d2theta_sign = 1;
      }
      finalize_result(out);
      return out;
    }

    const Real u = tan_half_degree(theta_degrees);
    if (!(u > 0.0L) || !std::isfinite(u)) {
      out.status = STATUS_INVALID_ANGLE;
      return out;
    }
    const Real log_value =
        static_cast<Real>(m) * (std::log(u) - std::log(norm->u0)) -
        norm->log_norm_factor;
    const Real sin_theta = sin_degree(theta_degrees);
    const Real cos_theta = cos_degree(theta_degrees);
    const Real d1_coefficient = static_cast<Real>(m) / sin_theta;

    const Real d2_coefficient = m == 1 ? 1.0L / (1.0L + cos_theta) :
        static_cast<Real>(m) * (m - cos_theta) / (sin_theta * sin_theta);
    const Real log_d2 = log_value + (m == 1 ? -std::log1p(cos_theta) :
        std::log(static_cast<Real>(m)) + std::log(m - cos_theta) - 2.0L*std::log(sin_theta));

    out.value = from_log(log_value);
    out.dtheta = from_log(log_value + std::log(d1_coefficient));
    out.d2theta = from_log(log_d2);
    out.log_abs_value = log_value;
    out.sign = 1;
    out.dtheta_over_value = d1_coefficient;
    out.log_abs_dtheta = log_value + std::log(d1_coefficient);
    out.dtheta_sign = 1;
    out.d2theta_over_value = d2_coefficient;
    out.log_abs_d2theta = log_d2;
    out.d2theta_sign = 1;

    if (!std::isfinite(log_value))
      out.status = STATUS_NONFINITE;
    finalize_result(out);
    return out;
  } catch (...) {
    out.status = STATUS_NONFINITE;
    return out;
  }
}

}

namespace Constant {
inline ConstantDetail::BasisResult Evaluate(double m, double theta, double theta0) {
  int order = 0;
  if (!ConstantDetail::integer_order(m, order)) {
    ConstantDetail::BasisResult result;
    result.status = ConstantDetail::STATUS_INVALID_ORDER;
    return result;
  }
  return ConstantDetail::evaluate(order, theta, theta0);
}

inline void Constant_1rank(double m, double theta, double theta0, double &value,
                           double &dtheta) {
  const auto result = Evaluate(m, theta, theta0);
  value = result.status == ConstantDetail::STATUS_OK ? result.value
      : std::numeric_limits<double>::quiet_NaN();
  dtheta = result.status == ConstantDetail::STATUS_OK ? result.dtheta
      : std::numeric_limits<double>::quiet_NaN();
}

inline void Constant_2rank(double m, double theta, double theta0, double &value,
                           double &dtheta, double &d2theta) {
  const auto result = Evaluate(m, theta, theta0);
  value = result.status == ConstantDetail::STATUS_OK ? result.value
      : std::numeric_limits<double>::quiet_NaN();
  dtheta = result.status == ConstantDetail::STATUS_OK ? result.dtheta
      : std::numeric_limits<double>::quiet_NaN();
  d2theta = result.status == ConstantDetail::STATUS_OK ? result.d2theta
      : std::numeric_limits<double>::quiet_NaN();
}
}

namespace Legendre {

inline void Legendre_2rank(double m, double n, double theta_degree,
                           double &value, double &dtheta, double &d2theta) {
  const auto result = EvaluateLegacySchmidt(
      m, n, theta_degree * DegreeTrig::radians_per_degree);
  if (result.status != 0) {
    value = dtheta = d2theta = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  value = static_cast<double>(result.value);
  dtheta = static_cast<double>(result.dtheta);
  d2theta = static_cast<double>(result.d2theta);
}
}
