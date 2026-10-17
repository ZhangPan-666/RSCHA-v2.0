#pragma once
#include "RSCHACompatibility.h"

// Preserve signed-log and compensated arithmetic in OpenMP builds as well.
#pragma float_control(precise, on, push)

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <complex>
#include <limits>
#include <unordered_map>

#include "CanonicalResidualLegendre.h"
#include "CRTable.h"
#include "CRNearFallback.h"

namespace Legendre {

namespace detail {

using real = long double;
constexpr real pi = 3.141592653589793238462643383279502884L;
constexpr real half_pi = pi / 2.0L;
constexpr real log2_value = 0.693147180559945309417232121458176568L;
constexpr real default_rel_tol = 2.0e-15L;
constexpr real default_abs_tol = 2.0e-18L;
constexpr int status_ok = 0;
constexpr int status_domain = 1;
constexpr int status_series = 2;
constexpr int status_integrator = 3;
constexpr int status_phase = 4;
constexpr int status_range = 6;

inline real sqr(real x) noexcept { return x * x; }

inline bool finite(real x) noexcept {
  return std::isfinite(static_cast<long double>(x));
}

inline real safe_exp(real x) noexcept {
  const real hi = std::log(std::numeric_limits<real>::max());
  const real lo = std::log(std::numeric_limits<real>::min());
  if (x > hi)
    return std::numeric_limits<real>::infinity();
  if (x < lo)
    return 0.0L;
  return std::exp(x);
}

// Positive-argument digamma: recur to the stable Bernoulli region.
inline real digamma_positive(real x) noexcept {
  real value = 0.0L;
  while (x < 12.0L) {
    value -= 1.0L / x;
    x += 1.0L;
  }
  const real inv = 1.0L / x;
  const real inv2 = inv * inv;
  value +=
      std::log(x) - 0.5L * inv -
      inv2 * (1.0L / 12.0L -
              inv2 * (1.0L / 120.0L -
                      inv2 * (1.0L / 252.0L -
                              inv2 * (1.0L / 240.0L -
                                      inv2 * (1.0L / 132.0L -
                                              inv2 * (691.0L / 32760.0L))))));
  return value;
}

inline double to_double(real x) noexcept {
  if (std::isnan(x))
    return std::numeric_limits<double>::quiet_NaN();
  if (x > static_cast<real>(DBL_MAX))
    return std::numeric_limits<double>::infinity();
  if (x < -static_cast<real>(DBL_MAX))
    return -std::numeric_limits<double>::infinity();
  if (std::fabs(x) < static_cast<real>(DBL_MIN))
    return 0.0;
  return static_cast<double>(x);
}

inline bool parse_order(double value, int &order) noexcept {
  if (!std::isfinite(value) || value < 0.0 || value > 1.0e8)
    return false;
  const double rounded = std::round(value);
  if (std::fabs(value - rounded) > 1.0e-12)
    return false;
  order = static_cast<int>(rounded);
  return true;
}

inline bool valid_mode_domain(int order, double degree, double theta) noexcept {
  return std::isfinite(degree) && std::isfinite(theta) &&
         degree >= static_cast<double>(order) && degree <= 1.0e8 &&
         theta > 0.0 && theta < static_cast<double>(pi);
}

struct EvalResult {
  real value = 0.0L;
  real dtheta = 0.0L;
  real d2theta = 0.0L;
  real log_abs_value = -std::numeric_limits<real>::infinity();
  real log_abs_value_lo = 0.0L;
  int sign = 0;
  real residual = 0.0L;
  real phase = 0.0L;
  real phase_derivative = 0.0L;
  real phase_kappa = 0.0L;
  int status = status_ok;
  real log_abs_dtheta = -std::numeric_limits<real>::infinity();
  real log_abs_d2theta = -std::numeric_limits<real>::infinity();
  real log_abs_dtheta_lo = 0.0L;
  real log_abs_d2theta_lo = 0.0L;
  int dtheta_sign = 0;
  int d2theta_sign = 0;
  real dtheta_over_value = std::numeric_limits<real>::quiet_NaN();
  real d2theta_over_value = std::numeric_limits<real>::quiet_NaN();
  unsigned range_flags = 0;
  const char* evaluation_method = "unclassified";
};

constexpr unsigned range_value_underflow = 1U << 0;
constexpr unsigned range_dtheta_underflow = 1U << 1;
constexpr unsigned range_d2theta_underflow = 1U << 2;
constexpr unsigned range_value_overflow = 1U << 3;
constexpr unsigned range_dtheta_overflow = 1U << 4;
constexpr unsigned range_d2theta_overflow = 1U << 5;

inline void classify_log_range(real log_abs, unsigned underflow_flag,
                               unsigned overflow_flag,
                               unsigned &flags) noexcept {
  if (!finite(log_abs))
    return;
  const real log_min =
      std::log(static_cast<real>(std::numeric_limits<double>::min()));
  const real log_max =
      std::log(static_cast<real>(std::numeric_limits<double>::max()));
  if (log_abs < log_min)
    flags |= underflow_flag;
  if (log_abs > log_max)
    flags |= overflow_flag;
}

inline void finalize_metadata(EvalResult &out) noexcept {
  if (out.status != status_ok)
    return;
  if (out.value != 0.0L) {
    out.sign = out.value > 0.0L ? 1 : -1;
    if (!finite(out.log_abs_value))
      out.log_abs_value = std::log(std::fabs(out.value));
    if (!finite(out.dtheta_over_value))
      out.dtheta_over_value = out.dtheta / out.value;
    if (!finite(out.d2theta_over_value))
      out.d2theta_over_value = out.d2theta / out.value;
  }
  if (out.dtheta != 0.0L) {
    out.dtheta_sign = out.dtheta > 0.0L ? 1 : -1;
    out.log_abs_dtheta = std::log(std::fabs(out.dtheta));
    out.log_abs_dtheta_lo=0.0L;
  } else if (finite(out.log_abs_value) && finite(out.dtheta_over_value) &&
             out.dtheta_over_value != 0.0L) {
    out.dtheta_sign = out.sign * (out.dtheta_over_value > 0.0L ? 1 : -1);
    const real increment = std::log(std::fabs(out.dtheta_over_value));
    const auto sum=crli::detail::DoubleDouble{static_cast<double>(out.log_abs_value),
        static_cast<double>(out.log_abs_value_lo)}+
        crli::detail::DoubleDouble{static_cast<double>(increment),0.0};
    out.log_abs_dtheta = sum.hi;
    out.log_abs_dtheta_lo = sum.lo;
  }
  if (out.d2theta != 0.0L) {
    out.d2theta_sign = out.d2theta > 0.0L ? 1 : -1;
    out.log_abs_d2theta = std::log(std::fabs(out.d2theta));
    out.log_abs_d2theta_lo=0.0L;
  } else if (finite(out.log_abs_value) && finite(out.d2theta_over_value) &&
             out.d2theta_over_value != 0.0L) {
    out.d2theta_sign = out.sign * (out.d2theta_over_value > 0.0L ? 1 : -1);
    const real increment = std::log(std::fabs(out.d2theta_over_value));
    const auto sum=crli::detail::DoubleDouble{static_cast<double>(out.log_abs_value),
        static_cast<double>(out.log_abs_value_lo)}+
        crli::detail::DoubleDouble{static_cast<double>(increment),0.0};
    out.log_abs_d2theta = sum.hi;
    out.log_abs_d2theta_lo = sum.lo;
  }
  classify_log_range(out.log_abs_value, range_value_underflow,
                     range_value_overflow, out.range_flags);
  classify_log_range(out.log_abs_dtheta, range_dtheta_underflow,
                     range_dtheta_overflow, out.range_flags);
  classify_log_range(out.log_abs_d2theta, range_d2theta_underflow,
                     range_d2theta_overflow, out.range_flags);
}

struct PolarData {
  real log_abs_L = 0.0L;
  real log_abs_L_kappa = 0.0L;
  int sign = 1;
  real r_t = 0.0L;
  real r_t_kappa = 0.0L;
  real d2theta_over_value = std::numeric_limits<real>::quiet_NaN();
  real L = 0.0L;
  real Ltheta = 0.0L;
  bool converged = false;
};

struct PhaseState {

  real a = 0.0L;
  real q = 0.0L;
  real alpha = 0.0L;
  real a_kappa = 0.0L;
  real q_kappa = 0.0L;
  real alpha_kappa = 0.0L;
};

struct MatchData {
  real theta = 0.0L;
  real delta = 0.0L;
  real delta_kappa = 0.0L;
  real log_C = 0.0L;
  real phi_match = 0.0L;
  int status = status_ok;
};

template <std::size_t N> struct TaylorJet {
  std::array<real, N + 1> c{};
};

template <std::size_t N>
inline TaylorJet<N> jet_add(const TaylorJet<N> &x,
                            const TaylorJet<N> &y) noexcept {
  TaylorJet<N> z;
  for (std::size_t i = 0; i <= N; ++i)
    z.c[i] = x.c[i] + y.c[i];
  return z;
}

template <std::size_t N>
inline TaylorJet<N> jet_scale(const TaylorJet<N> &x, real a) noexcept {
  TaylorJet<N> z;
  for (std::size_t i = 0; i <= N; ++i)
    z.c[i] = a * x.c[i];
  return z;
}

template <std::size_t N>
inline TaylorJet<N> jet_mul(const TaylorJet<N> &x,
                            const TaylorJet<N> &y) noexcept {
  TaylorJet<N> z;
  for (std::size_t i = 0; i <= N; ++i)
    for (std::size_t j = 0; j <= i; ++j)
      z.c[i] += x.c[j] * y.c[i - j];
  return z;
}

template <std::size_t N>
inline TaylorJet<N> jet_inv(const TaylorJet<N> &x) noexcept {
  TaylorJet<N> z;
  z.c[0] = 1.0L / x.c[0];
  for (std::size_t i = 1; i <= N; ++i) {
    real sum = 0.0L;
    for (std::size_t j = 1; j <= i; ++j)
      sum += x.c[j] * z.c[i - j];
    z.c[i] = -sum / x.c[0];
  }
  return z;
}

template <std::size_t N>
inline TaylorJet<N> jet_sqrt(const TaylorJet<N> &x) noexcept {
  TaylorJet<N> z;
  z.c[0] = std::sqrt(x.c[0]);
  for (std::size_t i = 1; i <= N; ++i) {
    real sum = 0.0L;
    for (std::size_t j = 1; j < i; ++j)
      sum += z.c[j] * z.c[i - j];
    z.c[i] = (x.c[i] - sum) / (2.0L * z.c[0]);
  }
  return z;
}

template <std::size_t N>
inline TaylorJet<N> jet_derivative(const TaylorJet<N> &x) noexcept {
  TaylorJet<N> z;
  for (std::size_t i = 0; i < N; ++i)
    z.c[i] = static_cast<real>(i + 1) * x.c[i + 1];
  return z;
}

template <std::size_t N> inline TaylorJet<N> sine_jet(real theta) noexcept {
  TaylorJet<N> z;
  const real s = std::sin(theta), c = std::cos(theta);
  real factorial = 1.0L;
  for (std::size_t i = 0; i <= N; ++i) {
    if (i > 0)
      factorial *= static_cast<real>(i);
    switch (i & 3U) {
    case 0:
      z.c[i] = s / factorial;
      break;
    case 1:
      z.c[i] = c / factorial;
      break;
    case 2:
      z.c[i] = -s / factorial;
      break;
    default:
      z.c[i] = -c / factorial;
      break;
    }
  }
  return z;
}

struct WkbDensity {
  real p = 0.0L;
  real ptheta = 0.0L;
  real correction = 0.0L;
  real relative_defect = std::numeric_limits<real>::infinity();
  bool valid = false;
};

template <std::size_t N>
inline TaylorJet<N> jet_prefix(const TaylorJet<18> &source) noexcept {
  TaylorJet<N> result;
  for (std::size_t i=0; i<=N; ++i) result.c[i]=source.c[i];
  return result;
}

template <std::size_t N>
inline bool wkb_density_update(const TaylorJet<18> &q, TaylorJet<18> &p,
                               real &previous_change, int iter,
                               bool &valid) noexcept {
  const auto pinv = jet_inv(jet_prefix<N>(p));
  TaylorJet<N> first, second;
  for (std::size_t i=0; i<=N; ++i) {
    first.c[i]=static_cast<real>(i+1)*p.c[i+1];
    second.c[i]=static_cast<real>(i+1)*
        (static_cast<real>(i+2)*p.c[i+2]);
  }
  const auto ratio1=jet_mul(first,pinv);
  const auto ratio2=jet_mul(second,pinv);
  auto corrected=jet_add(jet_prefix<N>(q),jet_scale(jet_mul(ratio1,ratio1),0.75L));
  corrected=jet_add(corrected,jet_scale(ratio2,-0.5L));
  if (!(corrected.c[0]>0.0L) || !finite(corrected.c[0])) {
    valid=false; return false;
  }
  const auto next=jet_sqrt(corrected);
  const real change=std::fabs(next.c[0]-p.c[0]);
  if (iter>=3 && change>previous_change) return false;
  for (std::size_t i=0; i<=N; ++i) p.c[i]=next.c[i];
  previous_change=change;
  return true;
}

// Table-free Kummer iteration, accepted only when its defect is at roundoff
// level.
inline WkbDensity wkb_density(real kappa, int m, real theta,
                              bool coalescent = false) noexcept {
  constexpr std::size_t order = 18;
  const auto st = sine_jet<order>(theta);
  const auto st2 = jet_mul(st, st);
  const auto csc2 = jet_inv(st2);
  TaylorJet<order> q;
  q.c[0] = kappa * kappa;
  const real eta2 = static_cast<real>(m) * m - 0.25L;
  q = jet_add(q, jet_scale(csc2, -eta2));
  if(coalescent) {
    const real cot=std::cos(theta)/std::sin(theta);
    q.c[0]=(kappa-m)*(kappa+m)+0.25L-eta2*cot*cot;
  }
  WkbDensity out;
  if (!(q.c[0] > 0.0L))
    return out;

  TaylorJet<order> p = jet_sqrt(q);
  real previous_change = std::numeric_limits<real>::infinity();
  bool valid=true;
  bool advance=wkb_density_update<16>(q,p,previous_change,0,valid);
  if(advance) advance=wkb_density_update<14>(q,p,previous_change,1,valid);
  if(advance) advance=wkb_density_update<12>(q,p,previous_change,2,valid);
  if(advance) advance=wkb_density_update<10>(q,p,previous_change,3,valid);
  if(advance) advance=wkb_density_update<8>(q,p,previous_change,4,valid);
  if(advance) advance=wkb_density_update<6>(q,p,previous_change,5,valid);
  if(advance) wkb_density_update<4>(q,p,previous_change,6,valid);
  if(!valid) return out;

  const auto pinv = jet_inv(p);
  const auto ratio1 = jet_mul(jet_derivative(p), pinv);
  const auto ratio2 = jet_mul(jet_derivative(jet_derivative(p)), pinv);
  const real source = 0.75L * ratio1.c[0] * ratio1.c[0] - 0.5L * ratio2.c[0];
  const real corrected_q = q.c[0] + source;
  const real defect = p.c[0] * p.c[0] - corrected_q;
  out.p = p.c[0];
  out.ptheta = p.c[1];
  // Form the small phase correction without subtracting O(kappa) values.
  out.correction = source / (std::sqrt(corrected_q) + std::sqrt(q.c[0]));
  out.relative_defect = std::fabs(defect) / q.c[0];
  out.valid = finite(out.p) && finite(out.ptheta) && out.p > 0.0L &&
              out.relative_defect <= 3.0e-13L;
  return out;
}

// Reduce the phase modulo 2*pi before large O(kappa) terms lose precision.
#pragma float_control(precise, on, push)
inline real reduced_leading_phase(real kappa, int m, real theta,
                                  real correction) noexcept {
  constexpr real half_hi = 1.57079632679489655800L;
  constexpr real half_lo = 6.12323399573676603587e-17L;
  constexpr real two_pi_hi = 6.28318530717958623200L;
  constexpr real two_pi_lo = 2.44929359829470641435e-16L;
  const real delta = half_hi - theta;
  const real b = delta - half_hi;
  const real delta_error = (half_hi - (delta - b)) - (theta + b);
  const real product = kappa * delta;
  const real product_error =
      std::fma(kappa, delta, -product) + kappa * (delta_error + half_lo);

  const real eta2 = static_cast<real>(m) * m - 0.25L;
  const real s = std::sin(theta), u = std::cos(theta);
  real small = 0.0L;
  if (eta2 >= 0.0L) {
    const real eta = std::sqrt(eta2), rho = eta / kappa;
    const real rad = std::sqrt(std::max(0.0L, s * s - rho * rho));
    const real s_minus_rad = (rho * rho) / (s + rad);
    const real dphi = std::atan2(u * s_minus_rad, rad * s + u * u);
    small = kappa * dphi - eta * std::atan2(rho * u, rad);
  } else {
    const real eta = std::sqrt(-eta2), h = eta / kappa;
    const real rad = std::sqrt(s * s + h * h);
    const real s_minus_rad = -(h * h) / (s + rad);
    const real dphi = std::atan2(u * s_minus_rad, rad * s + u * u);
    small = kappa * dphi + eta * std::atanh(h * u / rad);
  }
  const real phase_hi = -product;
  const real phase_lo = -(product_error + small + correction);
  const real quotient = std::nearbyint(phase_hi / two_pi_hi);
  return std::fma(-quotient, two_pi_hi, phase_hi) +
         (phase_lo - quotient * two_pi_lo);
}
#pragma float_control(pop)

template <class F>
inline bool gauss_phase_integral(F &&f, real a, real b, real &value,
                                 real tolerance, int depth = 0) {
  static constexpr real x8[4] = {0.183434642495649805L, 0.525532409916328986L,
                                 0.796666477413626740L, 0.960289856497536232L};
  static constexpr real w8[4] = {0.362683783378361983L, 0.313706645877887287L,
                                 0.222381034453374471L, 0.101228536290376259L};
  static constexpr real x16[8] = {
      0.0950125098376374402L, 0.281603550779258913L, 0.458016777657227386L,
      0.617876244402643748L,  0.755404408355003034L, 0.865631202387831744L,
      0.944575023073232576L,  0.989400934991649933L};
  static constexpr real w16[8] = {0.189450610455068496L, 0.182603415044923589L,
                                  0.169156519395002538L, 0.149595988816576733L,
                                  0.124628971255533872L, 0.095158511682492785L,
                                  0.062253523938647893L, 0.027152459411754095L};
  const real mid = 0.5L * (a + b), half = 0.5L * (b - a);
  real i8 = 0.0L, i16 = 0.0L;
  for (int i = 0; i < 4; ++i) {
    const real fm = f(mid - half * x8[i]), fp = f(mid + half * x8[i]);
    if (!finite(fm) || !finite(fp))
      return false;
    i8 += w8[i] * (fm + fp);
  }
  for (int i = 0; i < 8; ++i) {
    const real fm = f(mid - half * x16[i]), fp = f(mid + half * x16[i]);
    if (!finite(fm) || !finite(fp))
      return false;
    i16 += w16[i] * (fm + fp);
  }
  i8 *= half;
  i16 *= half;
  if (std::fabs(i16 - i8) <= tolerance) {
    value = i16;
    return true;
  }
  if(depth >= 12) return false;
  real left = 0.0L, right = 0.0L;
  if (!gauss_phase_integral(f, a, mid, left, 0.5L * tolerance, depth + 1))
    return false;
  if (!gauss_phase_integral(f, mid, b, right, 0.5L * tolerance, depth + 1))
    return false;
  value = left + right;
  return true;
}

inline real sinpi_reduced(real x) noexcept {
  const real reduced = std::remainder(x, 2.0L);
  const real integer = std::nearbyint(reduced);
  const real delta = reduced - integer;
  const real sign = integer == 0.0L ? 1.0L : -1.0L;
  return sign * std::sin(pi * delta);
}

inline real log_gamma_ratio_half(real x) noexcept {
  if (!(x > 0.0L))
    return std::numeric_limits<real>::quiet_NaN();
  if (x < 32.0L)
    return std::lgammal(x) - std::lgammal(x + 0.5L);
  const real inv = 1.0L / x, inv2 = inv * inv;
  const real correction =
      inv * (1.0L / 8.0L +
             inv2 * (-1.0L / 192.0L +
                     inv2 * (1.0L / 640.0L +
                             inv2 * (-17.0L / 14336.0L +
                                     inv2 * (31.0L / 18432.0L -
                                             inv2 * 691.0L / 180224.0L)))));
  return -0.5L * std::log(x) + correction;
}

inline void equator_legacy_values(real nu, int m, real &value,
                                  real &derivative) noexcept {
  const real mr = static_cast<real>(m);
  const real A = 0.5L * (nu - mr + 1.0L);
  const real B = 0.5L * (nu + mr + 1.0L);
  const real longitude = (m == 0) ? 0.0L : std::log(2.0L);
  const real phase0 = sinpi_reduced(0.5L * (nu - mr + 1.0L));
  if (phase0 == 0.0L)
    value = 0.0L;
  else {
    const real logabs =
        -0.5L * std::log(pi) + 0.5L * longitude + std::log(std::fabs(phase0)) +
        0.5L * (log_gamma_ratio_half(A) + log_gamma_ratio_half(B));
    const int sign = (phase0 > 0.0L ? 1 : -1);
    value = static_cast<real>(sign) * safe_exp(logabs);
  }
  const real phase1 = sinpi_reduced(0.5L * (nu - mr));
  if (phase1 == 0.0L)
    derivative = 0.0L;
  else {
    const real logabs =
        std::log(nu + mr) - 0.5L * std::log(pi) + 0.5L * longitude +
        std::log(std::fabs(phase1)) +
        0.5L * (-log_gamma_ratio_half(A) + log_gamma_ratio_half(B - 0.5L) -
                std::log(B - 0.5L));
    const int rawsign = (phase1 > 0.0L ? 1 : -1);
    derivative = -static_cast<real>(rawsign) * safe_exp(logabs);
  }
}

// Exact degree derivatives of the same equatorial data used by the value
// evaluator.  Keeping these formulas together prevents the sensitivity IVP
// from being initialized on a different normalization branch.
inline bool equator_legacy_sensitivities(real nu, int m, real &value,
                                         real &derivative, real &value_nu,
                                         real &derivative_nu) noexcept {
  const real mr = static_cast<real>(m);
  const real A = 0.5L * (nu - mr + 1.0L);
  const real B = 0.5L * (nu + mr + 1.0L);
  if (!(A > 0.0L && B > 0.5L))
    return false;
  const real longitude = (m == 0) ? 0.0L : std::log(2.0L);
  const real parity = 1.0L;
  const auto ratio_psi = [](real x) {
    return digamma_positive(x) - digamma_positive(x + 0.5L);
  };

  const real x0 = 0.5L * (nu - mr + 1.0L);
  const real phase0 = sinpi_reduced(x0);
  const real phase0_nu = 0.5L * pi * sinpi_reduced(std::remainder(x0, 2.0L) + 0.5L);
  const real base0 = -0.5L * std::log(pi) + 0.5L * longitude +
                     0.5L * (log_gamma_ratio_half(A) +
                             log_gamma_ratio_half(B));
  const real base0_nu = 0.25L * (ratio_psi(A) + ratio_psi(B));
  const real amplitude0 = safe_exp(base0);
  value = parity * phase0 * amplitude0;
  value_nu = parity * amplitude0 * (phase0_nu + phase0 * base0_nu);

  const real x1 = 0.5L * (nu - mr);
  const real phase1 = sinpi_reduced(x1);
  const real phase1_nu = 0.5L * pi * sinpi_reduced(std::remainder(x1, 2.0L) + 0.5L);
  const real base1 = std::log(nu + mr) - 0.5L * std::log(pi) +
                     0.5L * longitude +
                     0.5L * (-log_gamma_ratio_half(A) +
                             log_gamma_ratio_half(B - 0.5L) -
                             std::log(B - 0.5L));
  const real base1_nu = 1.0L / (nu + mr) +
      0.25L * (-ratio_psi(A) + ratio_psi(B - 0.5L)) -
      0.25L / (B - 0.5L);
  const real amplitude1 = safe_exp(base1);
  derivative = -parity * phase1 * amplitude1;
  derivative_nu = -parity * amplitude1 *
                  (phase1_nu + phase1 * base1_nu);
  return finite(value) && finite(derivative) && finite(value_nu) &&
         finite(derivative_nu);
}

inline bool wkb_evaluate_legacy(real nu, int m, real theta, EvalResult &out) {
  const real kappa = nu + 0.5L;
  if (kappa < 32.0L || !(theta > 0.0L && theta <= half_pi))
    return false;
  const real eta2 = static_cast<real>(m) * m - 0.25L;
  const real qtheta = kappa * kappa - eta2 / sqr(std::sin(theta));
  if (!(qtheta > 0.0L) || qtheta < 0.04L * kappa * kappa)
    return false;
  const WkbDensity pt = wkb_density(kappa, m, theta);
  const WkbDensity p0 = wkb_density(kappa, m, half_pi);
  if (!pt.valid || !p0.valid)
    return false;

  const auto correction = [=](real x) -> real {
    const WkbDensity d = wkb_density(kappa, m, x);
    if (!d.valid)
      return std::numeric_limits<real>::quiet_NaN();
    // Integrate the stable correction already formed by rationalization.
    // Subtracting two O(kappa) densities injects roundoff into the quadrature.
    return d.correction;
  };
  real corr = 0.0L;
  if (!gauss_phase_integral(correction, theta, half_pi, corr, 5.0e-13L))
    return false;
  const real phase = reduced_leading_phase(kappa, m, theta, corr);

  const real mr = static_cast<real>(m);
  real y0 = 0.0L, dy0 = 0.0L;
  equator_legacy_values(nu, m, y0, dy0);
  const real rootp0 = std::sqrt(p0.p);
  const real c1 = y0 * rootp0;
  const real c2 = dy0 / rootp0;
  const real sn = std::sin(phase), cs = std::cos(phase);
  const real shape = c1 * cs + c2 * sn;
  const real shape_perp = -c1 * sn + c2 * cs;
  const real roott = std::sqrt(pt.p);
  const real U = shape / roott;
  const real a = -0.5L * pt.ptheta / pt.p;
  const real Ud = (pt.p * shape_perp + a * shape) / roott;
  const real st = std::sin(theta), ct = std::cos(theta), cot = ct / st;
  EvalResult candidate;
  candidate.value = U / std::sqrt(st);
  candidate.dtheta = (Ud - 0.5L * cot * U) / std::sqrt(st);
  const real coeff = nu * (nu + 1.0L) - mr * mr / (st * st);
  candidate.d2theta = -cot * candidate.dtheta - coeff * candidate.value;
  if (!finite(candidate.value) || !finite(candidate.dtheta) ||
      !finite(candidate.d2theta))
    return false;
  candidate.log_abs_value = candidate.value == 0.0L
                                ? -std::numeric_limits<real>::infinity()
                                : std::log(std::fabs(candidate.value));
  candidate.sign = (candidate.value > 0.0L) - (candidate.value < 0.0L);
  candidate.phase = phase;
  candidate.phase_derivative = pt.p;
  candidate.residual = std::max(pt.relative_defect, p0.relative_defect);
  candidate.status = status_ok;
  out = candidate;
  return true;
}

inline bool wkb_evaluate_fast(real nu, int m, real theta,
                              EvalResult &out) noexcept {
  const real kappa = nu + 0.5L;
  if (kappa < 512.0L || kappa * theta < 1000.0L ||
      !(theta > 0.0L && theta <= half_pi))
    return false;
  const real eta2 = static_cast<real>(m) * m - 0.25L;
  const real st = std::sin(theta);
  const real qtheta = kappa * kappa - eta2 / (st * st);
  if (!(qtheta > 0.04L * kappa * kappa))
    return false;
  const WkbDensity pt = wkb_density(kappa, m, theta);
  const WkbDensity p0 = wkb_density(kappa, m, half_pi);
  if (!pt.valid || !p0.valid)
    return false;
  static constexpr real nodes[8] = {
      0.0950125098376374402L, 0.281603550779258913L,
      0.458016777657227386L, 0.617876244402643748L,
      0.755404408355003034L, 0.865631202387831744L,
      0.944575023073232576L, 0.989400934991649933L};
  static constexpr real weights[8] = {
      0.189450610455068496L, 0.182603415044923589L,
      0.169156519395002538L, 0.149595988816576733L,
      0.124628971255533872L, 0.095158511682492785L,
      0.062253523938647893L, 0.027152459411754095L};
  const real middle = 0.5L * (theta + half_pi);
  const real half = 0.5L * (half_pi - theta);
  real correction = 0.0L;
  for (int i = 0; i < 8; ++i)
    for (real side : {-1.0L, 1.0L}) {
      const real point = middle + side * half * nodes[i];
      const WkbDensity density = wkb_density(kappa, m, point);
      if (!density.valid)
        return false;
      correction += weights[i] * density.correction;
    }
  correction *= half;
  static constexpr real coarse_nodes[4] = {
      0.183434642495649805L, 0.525532409916328986L,
      0.796666477413626740L, 0.960289856497536232L};
  static constexpr real coarse_weights[4] = {
      0.362683783378361983L, 0.313706645877887287L,
      0.222381034453374471L, 0.101228536290376259L};
  real coarse_correction = 0.0L;
  for (int i = 0; i < 4; ++i)
    for (real side : {-1.0L, 1.0L}) {
      const real point = middle + side * half * coarse_nodes[i];
      const WkbDensity density = wkb_density(kappa, m, point);
      if (!density.valid)
        return false;
      coarse_correction +=
          coarse_weights[i] * density.correction;
    }
  coarse_correction *= half;
  if (std::fabs(correction - coarse_correction) >
          2.0e-12L * (1.0L + std::fabs(correction)) &&
      !(kappa >= 1000000.0L && theta >= 88.0L * pi / 180.0L))
    return false;
  const real phase = reduced_leading_phase(kappa, m, theta, correction);
  real y0 = 0.0L, dy0 = 0.0L;
  equator_legacy_values(nu, m, y0, dy0);
  const real rootp0 = std::sqrt(p0.p);
  const real c1 = y0 * rootp0;
  const real c2 = dy0 / rootp0;
  const real sn = std::sin(phase), cs = std::cos(phase);
  const real shape = c1 * cs + c2 * sn;
  const real shape_perp = -c1 * sn + c2 * cs;
  const real roott = std::sqrt(pt.p);
  const real U = shape / roott;
  const real amplitude_rate = -0.5L * pt.ptheta / pt.p;
  const real Ud = (pt.p * shape_perp + amplitude_rate * shape) / roott;
  const real ct = std::cos(theta), cot = ct / st;
  out.value = U / std::sqrt(st);
  out.dtheta = (Ud - 0.5L * cot * U) / std::sqrt(st);
  const real mr = static_cast<real>(m);
  const real coefficient = nu * (nu + 1.0L) - mr * mr / (st * st);
  out.d2theta = -cot * out.dtheta - coefficient * out.value;
  if (!finite(out.value) || !finite(out.dtheta) || !finite(out.d2theta))
    return false;
  out.log_abs_value = out.value == 0.0L
                          ? -std::numeric_limits<real>::infinity()
                          : std::log(std::fabs(out.value));
  out.sign = (out.value > 0.0L) - (out.value < 0.0L);
  out.phase = phase;
  out.phase_derivative = pt.p;
  out.residual = std::max(pt.relative_defect, p0.relative_defect);
  out.status = status_ok;
  return true;
}

template <std::size_t N, class RHS>
bool dopri5_integrate(real x0, real x1, std::array<real, N> &y, RHS &&rhs,
                      real rel_tol = default_rel_tol,
                      real abs_tol = default_abs_tol, real max_step = 0.05L,
                      int max_steps = 200000) {
  if (x0 == x1)
    return true;
  const real direction = (x1 > x0) ? 1.0L : -1.0L;
  real h = direction *
           std::min(max_step, std::max(std::fabs(x1 - x0) / 16.0L, 1.0e-8L));
  real x = x0;

  std::array<real, N> k1{}, k2{}, k3{}, k4{}, k5{}, k6{}, k7{};
  std::array<real, N> yt{}, y5{}, y4{};

  for (int step = 0; step < max_steps; ++step) {
    if ((direction > 0.0L && x + h > x1) || (direction < 0.0L && x + h < x1))
      h = x1 - x;

    rhs(x, y, k1);
    for (std::size_t i = 0; i < N; ++i)
      yt[i] = y[i] + h * (1.0L / 5.0L) * k1[i];
    rhs(x + h * (1.0L / 5.0L), yt, k2);

    for (std::size_t i = 0; i < N; ++i)
      yt[i] = y[i] + h * ((3.0L / 40.0L) * k1[i] + (9.0L / 40.0L) * k2[i]);
    rhs(x + h * (3.0L / 10.0L), yt, k3);

    for (std::size_t i = 0; i < N; ++i)
      yt[i] = y[i] + h * ((44.0L / 45.0L) * k1[i] - (56.0L / 15.0L) * k2[i] +
                          (32.0L / 9.0L) * k3[i]);
    rhs(x + h * (4.0L / 5.0L), yt, k4);

    for (std::size_t i = 0; i < N; ++i)
      yt[i] = y[i] +
              h * ((19372.0L / 6561.0L) * k1[i] - (25360.0L / 2187.0L) * k2[i] +
                   (64448.0L / 6561.0L) * k3[i] - (212.0L / 729.0L) * k4[i]);
    rhs(x + h * (8.0L / 9.0L), yt, k5);

    for (std::size_t i = 0; i < N; ++i)
      yt[i] =
          y[i] + h * ((9017.0L / 3168.0L) * k1[i] - (355.0L / 33.0L) * k2[i] +
                      (46732.0L / 5247.0L) * k3[i] + (49.0L / 176.0L) * k4[i] -
                      (5103.0L / 18656.0L) * k5[i]);
    rhs(x + h, yt, k6);

    for (std::size_t i = 0; i < N; ++i)
      y5[i] =
          y[i] + h * ((35.0L / 384.0L) * k1[i] + (500.0L / 1113.0L) * k3[i] +
                      (125.0L / 192.0L) * k4[i] - (2187.0L / 6784.0L) * k5[i] +
                      (11.0L / 84.0L) * k6[i]);
    rhs(x + h, y5, k7);

    for (std::size_t i = 0; i < N; ++i)
      y4[i] = y[i] +
              h * ((5179.0L / 57600.0L) * k1[i] + (7571.0L / 16695.0L) * k3[i] +
                   (393.0L / 640.0L) * k4[i] - (92097.0L / 339200.0L) * k5[i] +
                   (187.0L / 2100.0L) * k6[i] + (1.0L / 40.0L) * k7[i]);

    real err_norm = 0.0L;
    for (std::size_t i = 0; i < N; ++i) {
      const real scale =
          abs_tol + rel_tol * std::max(std::fabs(y[i]), std::fabs(y5[i]));
      const real e = std::fabs(y5[i] - y4[i]) / scale;
      err_norm = std::max(err_norm, e);
    }

    if (!finite(err_norm))
      return false;

    if (err_norm <= 1.0L) {
      x += h;
      y = y5;
      if (x == x1)
        return true;
    }

    real factor = 0.9L * std::pow(std::max(err_norm, 1.0e-16L), -0.2L);
    factor = RSCHACompat::clamp(factor, 0.1L, 5.0L);
    h *= factor;
    if (std::fabs(h) > max_step)
      h = direction * max_step;
    if (std::fabs(h) < 32.0L * std::numeric_limits<real>::epsilon() *
                           std::max(1.0L, std::fabs(x)))
      return false;
  }
  return false;
}

inline PolarData polar_series(real kappa, int m, real theta,
                              real tolerance = 5.0e-19L,
                              int max_terms = 20000) {
  PolarData out;
  if (!(theta > 0.0L && theta < pi) || m < 0 || !(kappa > 0.0L))
    return out;

  const real s = std::sin(theta);
  const real c = std::cos(theta);
  const real z = sqr(std::sin(theta / 2.0L));
  const real dz = 0.5L * s;

  real sum = 1.0L, sum_k = 0.0L;
  real sum_z = 0.0L, sum_z_k = 0.0L, sum_zz = 0.0L;
  real comp_sum = 0.0L, comp_sum_k = 0.0L, comp_sum_z = 0.0L,
       comp_sum_z_k = 0.0L, comp_sum_zz = 0.0L;
  real abs_sum = 1.0L;
  real term = 1.0L, term_k = 0.0L;
  real term_z = 0.0L, term_zz = 0.0L, term_z_k = 0.0L;
  bool converged = false;

  for (int j = 1; j <= max_terms && !converged; ++j) {
    const real jj = static_cast<real>(j);
    const real den = jj * (static_cast<real>(m) + jj);
    const real center = static_cast<real>(m) + jj - 0.5L;
    const real A = (center - kappa) * (center + kappa);
    const real ratio = (A / den) * z;
    const real dratio = (-2.0L * kappa / den) * z;

    const real previous = term;
    const real previous_k = term_k;
    const real previous_z = term_z, previous_zz = term_zz;
    const real previous_z_k = term_z_k;
    term = previous * ratio;
    term_k = previous_k * ratio + previous * dratio;
    term_z = (previous_z * z + previous) * (A / den);
    term_zz = (previous_zz * z + 2.0L * previous_z) * (A / den);
    term_z_k = (previous_z_k * z + previous_k) * (A / den) +
               (previous_z * z + previous) * (-2.0L * kappa / den);

    auto compensated_add = [](real value, real &total,
                              real &compensation) noexcept {
      const real adjusted = value - compensation;
      const real updated = total + adjusted;
      compensation = (updated - total) - adjusted;
      total = updated;
    };
    compensated_add(term, sum, comp_sum);
    compensated_add(term_k, sum_k, comp_sum_k);
    abs_sum += std::fabs(term);
    compensated_add(term_z, sum_z, comp_sum_z);
    compensated_add(term_z_k, sum_z_k, comp_sum_z_k);
    compensated_add(term_zz, sum_zz, comp_sum_zz);

    const real next_center = static_cast<real>(m) + jj + 0.5L;
    const real next_den = (jj + 1.0L) * (static_cast<real>(m) + jj + 1.0L);
    const real next_ratio_abs =
        std::fabs((next_center * next_center - kappa * kappa) / next_den * z);
    const real tail = (next_ratio_abs < 0.95L)
                          ? std::fabs(term) * next_ratio_abs /
                                std::max(1.0e-30L, 1.0L - next_ratio_abs)
                          : std::numeric_limits<real>::infinity();

    const real cancellation =
        abs_sum / std::max(std::fabs(sum), std::numeric_limits<real>::min());
    if (j > 8 && tail <= tolerance * std::max(std::fabs(sum), 1.0L) &&
        std::fabs(term_k) <=
            50.0L * tolerance * std::max(std::fabs(sum_k), 1.0L) &&
        cancellation < 1.0e7L) {
      converged = true;
    }
    if (!finite(sum) || !finite(sum_k) || !finite(term))
      break;
  }

  if (!converged || sum == 0.0L || !finite(sum))
    return out;

  const real log_abs_L =
      static_cast<real>(m) * std::log(s) + std::log(std::fabs(sum));
  const int sign = (sum > 0.0L) ? 1 : -1;
  const real ratio_z = sum_z / sum;
  const real ratio_z_k = (sum_z_k * sum - sum_z * sum_k) / (sum * sum);
  const real Ltheta_over_L = static_cast<real>(m) * c / s + dz * ratio_z;
  const real r_t = s * Ltheta_over_L;
  const real h_t = s * dz * ratio_z_k;

  out.log_abs_L = log_abs_L;
  out.log_abs_L_kappa = sum_k / sum;
  out.sign = sign;
  out.r_t = r_t;
  out.r_t_kappa = h_t;
  out.d2theta_over_value =
      (m > 1 ? static_cast<real>(m) * (m - 1.0L) * sqr(c / s) : 0.0L) -
      static_cast<real>(m) + (static_cast<real>(m) + 0.5L) * c * ratio_z +
      0.25L * s * s * sum_zz / sum;
  const real mag = safe_exp(log_abs_L);
  out.L = static_cast<real>(sign) * mag;
  out.Ltheta = out.L * Ltheta_over_L;
  out.converged = true;
  return out;
}

inline real choose_series_start(real kappa, int m, real target) {
  const real scale = std::max(kappa, static_cast<real>(m + 1));
  real theta = std::min<real>(1.0e-2L, 0.20L / scale);
  theta = std::max<real>(theta, 2.0e-8L);
  theta = std::min(theta, target * 0.25L);
  theta = std::max<real>(theta, std::min(target * 0.5L, 2.0e-10L));
  return theta;
}

inline real choose_match_angle(real kappa, int m, real start) {
  real match = start;
  const real eta2 = static_cast<real>(m) * m - 0.25L;
  if (m >= 1 && eta2 > 0.0L) {
    const real x = std::sqrt(eta2) / kappa;
    const real turn = (x < 1.0L) ? std::asin(x) : half_pi;

    match = std::max(start, turn + 0.08L / kappa);
  } else {
    match = std::max(start, 0.60L / kappa);
  }

  match = std::min(match, 0.98L * half_pi);
  return std::max(match, start);
}

inline bool propagate_logarithmic(real kappa, int m, real theta0, real theta1,
                                  PolarData &state) {
  if (theta1 == theta0)
    return true;
  const real t0 = std::log(std::tan(theta0 / 2.0L));
  const real t1 = std::log(std::tan(theta1 / 2.0L));
  std::array<real, 4> y{state.r_t, state.r_t_kappa, state.log_abs_L,
                        state.log_abs_L_kappa};
  auto rhs = [=](real t, const std::array<real, 4> &v, std::array<real, 4> &f) {
    const real sech = 1.0L / std::cosh(t);
    const real sech2 = sech * sech;
    const real r = v[0];
    const real h = v[1];
    f[0] = static_cast<real>(m) * m - (kappa * kappa - 0.25L) * sech2 - r * r;
    f[1] = -2.0L * kappa * sech2 - 2.0L * r * h;
    f[2] = r;
    f[3] = h;
  };
  const bool ok = dopri5_integrate(t0, t1, y, rhs, 2.0e-15L, 2.0e-18L, 0.12L);
  if (!ok)
    return false;
  state.r_t = y[0];
  state.r_t_kappa = y[1];
  state.d2theta_over_value = std::numeric_limits<real>::quiet_NaN();
  state.log_abs_L = y[2];
  state.log_abs_L_kappa = y[3];
  state.sign = 1;
  state.L = safe_exp(state.log_abs_L);
  state.Ltheta = state.L * state.r_t / std::sin(theta1);
  return true;
}

inline PhaseState phase_anchor(real kappa, int m) {
  // Exact equatorial phase data; a local WKB seed selects the wrong branch.
  const real nu = kappa - 0.5L;
  const real mr = static_cast<real>(m);
  const real a1 = 0.5L * (nu - mr + 2.0L);
  const real a2 = 0.5L * (nu + mr + 2.0L);
  const real b1 = 0.5L * (nu - mr + 1.0L);
  const real b2 = 0.5L * (nu + mr + 1.0L);
  const real log_p0 = log2_value + std::lgamma(a1) + std::lgamma(a2) -
                      std::lgamma(b1) - std::lgamma(b2);
  PhaseState s;
  s.a = 0.0L;
  s.q = log_p0;
  s.alpha = 0.0L;
  s.a_kappa = 0.0L;
  s.q_kappa = 0.5L * (digamma_positive(a1) + digamma_positive(a2) -
                      digamma_positive(b1) - digamma_positive(b2));
  s.alpha_kappa = 0.0L;
  return s;
}

inline bool integrate_phase(real kappa, int m, real from, real to,
                            PhaseState &s) {
  if (from == to)
    return true;
  std::array<real, 6> y{s.a, s.q, s.alpha, s.a_kappa, s.q_kappa, s.alpha_kappa};
  const real eta2 = static_cast<real>(m) * m - 0.25L;
  auto rhs = [=](real theta, const std::array<real, 6> &v,
                 std::array<real, 6> &f) {
    const real st = std::sin(theta);
    const real csc2 = 1.0L / (st * st);
    const real Q = kappa * kappa - eta2 * csc2;
    const real p2 = safe_exp(2.0L * v[1]);
    const real p = safe_exp(v[1]);
    f[0] = -Q - v[0] * v[0] + p2;
    f[1] = -2.0L * v[0];
    f[2] = p;
    f[3] = -2.0L * kappa - 2.0L * v[0] * v[3] + 2.0L * p2 * v[4];
    f[4] = -2.0L * v[3];
    f[5] = p * v[4];
  };
  const bool ok =
      dopri5_integrate(from, to, y, rhs, 3.0e-15L, 3.0e-18L, 0.035L);
  if (!ok)
    return false;
  s.a = y[0];
  s.q = y[1];
  s.alpha = y[2];
  s.a_kappa = y[3];
  s.q_kappa = y[4];
  s.alpha_kappa = y[5];
  return finite(s.a) && finite(s.q) && finite(s.alpha);
}

inline MatchData build_match(real theta_match, const PolarData &polar,
                             const PhaseState &phase) {
  MatchData out;
  out.theta = theta_match;
  const real st = std::sin(theta_match);
  const real ct = std::cos(theta_match);
  const real g = polar.r_t / st + 0.5L * ct / st;
  const real gk = polar.r_t_kappa / st;
  const real p = safe_exp(phase.q);
  if (!(p > 0.0L) || !finite(p)) {
    out.status = status_phase;
    return out;
  }
  const real x = g - phase.a;
  real phi = std::atan2(p, x);
  if (phi <= 0.0L)
    phi += pi;
  const real pk = p * phase.q_kappa;
  const real den = x * x + p * p;
  const real phi_k = (x * pk - p * (gk - phase.a_kappa)) / den;

  const real logU = polar.log_abs_L + 0.5L * std::log(st);
  const real sin_phi = std::sin(phi);
  if (!(sin_phi > 0.0L) || !finite(phi_k)) {
    out.status = status_phase;
    return out;
  }

  out.phi_match = phi;
  out.delta = phi - phase.alpha;
  out.delta_kappa = phi_k - phase.alpha_kappa;
  out.log_C = logU + 0.5L * phase.q - std::log(sin_phi);
  out.status = status_ok;
  return out;
}

inline crli::detail::DoubleDouble pole_log_gamma(crli::detail::DoubleDouble x) {
  using namespace crli::detail;
  if(x.hi<16.0)return {std::lgamma(x.hi+x.lo),0.0};
  const double inverse=1.0/x.hi,z=inverse*inverse;
  const double correction=inverse*(1.0/12.0+z*(-1.0/360.0+z*(1.0/1260.0+
      z*(-1.0/1680.0+z*(1.0/1188.0+z*(-691.0/360360.0+z/156.0))))));
  return (x-DoubleDouble{0.5,0.0})*dd_log_positive(x)-
      x+DoubleDouble{0.91893853320467274178,0.0}+DoubleDouble{correction,0.0};
}

inline crli::detail::DoubleDouble pole_log_scale(real nu, int m) {
  using namespace crli::detail;
  const DoubleDouble degree{static_cast<double>(nu),0.0}, order{static_cast<double>(m),0.0};
  const DoubleDouble logarithm_two = dd_log_positive(DoubleDouble{2.0,0.0});
  return (pole_log_gamma(degree+order+DoubleDouble{1.0,0.0}) -
          pole_log_gamma(degree-order+DoubleDouble{1.0,0.0}))*0.5 -
      pole_log_gamma(order+DoubleDouble{1.0,0.0}) - logarithm_two*static_cast<double>(m) +
      (m == 0 ? DoubleDouble{} : logarithm_two*0.5);
}

inline real legacy_log_scale(real nu, int m) {
  const auto value = pole_log_scale(nu, m);
  return value.hi + value.lo;
}

inline bool evaluate_convergent_pole_schmidt(real nu,int m,real theta,EvalResult& out) {
  if(nu<90000.0L || m<1000 || !(theta>0.0L && theta<0.5L))return false;
  const real sh=std::sin(theta/2),z=sh*sh,kappa=nu+0.5L;
  if(kappa*kappa*z/(m+1.0L)>1.0L)return false;
  real term=1.0L,sum=1.0L,first=0.0L,second=0.0L;int quiet=0;
  for(int k=1;k<=256;k++) {
    term*=((k-1.0L-nu)*(k+nu)/((k+m)*static_cast<real>(k)))*z;
    sum+=term;first+=k*term;second+=k*(k-1.0L)*term;
    if(std::fabs(term)*k*k<=2.0e-19L*std::fabs(sum))++quiet;else quiet=0;
    if(quiet>=3)break;
    if(k==256 || !finite(sum))return false;
  }
  if(!(sum>0.0L))return false;
  using namespace crli::detail;
  const DoubleDouble degree{static_cast<double>(nu),0.0},order{static_cast<double>(m),0.0},one{1.0,0.0};
  const DoubleDouble logarithm=(pole_log_gamma(degree+order+one)-
      pole_log_gamma(degree-order+one))*0.5-
      pole_log_gamma(order+one)+dd_log_positive(DoubleDouble{static_cast<double>(std::tan(theta/2)),0.0})*static_cast<double>(m)+
      DoubleDouble{0.5*static_cast<double>(log2_value)+static_cast<double>(std::log(sum)),0.0};
  const real sine=std::sin(theta),cosine=std::cos(theta),zratio=std::cos(theta/2)/sh;
  const real h1=first/sum,h2=second/sum;
  out.log_abs_value=logarithm.hi;out.log_abs_value_lo=logarithm.lo;out.sign=1;
  out.dtheta_over_value=m/sine+h1*zratio;
  out.d2theta_over_value=m*(m-cosine)/(sine*sine)+2*m/sine*h1*zratio+
      h2*zratio*zratio+h1*cosine/(2*z);
  out.value=safe_exp(out.log_abs_value+out.log_abs_value_lo);
  out.dtheta=out.value*out.dtheta_over_value;out.d2theta=out.value*out.d2theta_over_value;
  out.status=status_ok;
  return finite(out.log_abs_value)&&finite(out.dtheta_over_value)&&finite(out.d2theta_over_value);
}

inline bool evaluate_integer_near_order(real nu, int m, real theta,
                                        EvalResult &out) {
  if (!(theta > 0.0L && theta < pi) || m < 0)
    return false;
  const real rounded = std::round(nu);
  const real gap = rounded - static_cast<real>(m);
  if (std::fabs(nu - rounded) > 1.0e-12L || rounded < m ||
      (gap > 4096.0L && !(gap <= 12000.0L && theta >= 1.4L) &&
       !(gap <= 20000.0L && static_cast<real>(m) / (nu + 0.5L) > 0.9L)))
    return false;
  const int degree = static_cast<int>(rounded);
  const real st = std::sin(theta), ct = std::cos(theta), cot = ct / st;
  real log_base = 0.0L;
  crli::detail::DoubleDouble retained_log_base{0.0,0.0};
  if (m > 0) {
    const real mr = static_cast<real>(m);
    if (m >= 64) {
      const real inverse = 1.0L / mr;
      const real inverse2 = inverse * inverse;
      log_base = 0.5L * log2_value -
                 0.25L * (std::log(pi) + std::log(mr)) -
                 inverse / 16.0L + inverse * inverse2 / 384.0L -
                 inverse * inverse2 * inverse2 / 1280.0L;
    } else {
      log_base = 0.5L * std::lgamma(2.0L * mr + 1.0L) -
                 mr * log2_value - std::lgamma(mr + 1.0L) +
                 0.5L * log2_value;
    }
    // At m=1e8 a single rounded log(sin(theta))*m loses useful digits.
    // Keep both product and sum tails for signed-log output after underflow.
    using namespace crli::detail;
    const DoubleDouble sine_squared=DoubleDouble{1.0,0.0}-
        DoubleDouble{static_cast<double>(ct),0.0}*DoubleDouble{static_cast<double>(ct),0.0};
    const DoubleDouble log_sine=theta>0.5L && theta<pi-0.5L ? dd_log_positive(sine_squared)*0.5
        : dd_log_positive(DoubleDouble{static_cast<double>(st),0.0});
    retained_log_base=DoubleDouble{static_cast<double>(log_base),0.0}+log_sine*static_cast<double>(mr);
    log_base=retained_log_base.hi;
  }

  real previous = 0.0L, value = 1.0L;
  real previous_d = 0.0L, derivative = static_cast<real>(m) * cot;
  real previous_d2 = 0.0L;
  real second_derivative = static_cast<real>(m) *
      ((static_cast<real>(m) - 1.0L) * cot * cot - 1.0L);
  real recurrence_log_scale = 0.0L;
  for (int j = m; j < degree; ++j) {
    const real jr = static_cast<real>(j);
    const real denominator =
        (jr - m + 1.0L) * (jr + m + 1.0L);
    const real a = (2.0L * jr + 1.0L) / std::sqrt(denominator);
    const real b = (j == m)
                       ? 0.0L
                       : std::sqrt((jr + m) * (jr - m) / denominator);
    const real next = a * ct * value - b * previous;
    const real next_d = a * (-st * value + ct * derivative) -
                        b * previous_d;
    const real next_d2 =
        a * (-ct * value - 2.0L * st * derivative +
             ct * second_derivative) - b * previous_d2;
    previous = value;
    value = next;
    previous_d = derivative;
    derivative = next_d;
    previous_d2 = second_derivative;
    second_derivative = next_d2;
    const real magnitude = std::max(
        {std::fabs(previous), std::fabs(value), std::fabs(previous_d),
         std::fabs(derivative), std::fabs(previous_d2),
         std::fabs(second_derivative)});
    if (!finite(magnitude) || magnitude == 0.0L)
      return false;
    if (magnitude > 1.0e100L || magnitude < 1.0e-100L) {
      previous /= magnitude;
      value /= magnitude;
      previous_d /= magnitude;
      derivative /= magnitude;
      previous_d2 /= magnitude;
      second_derivative /= magnitude;
      recurrence_log_scale += std::log(magnitude);
    }
  }

  auto scaled = [&](real coefficient) {
    if (coefficient == 0.0L)
      return 0.0L;
    return std::copysign(
        safe_exp(log_base + recurrence_log_scale +
                 std::log(std::fabs(coefficient))), coefficient);
  };
  out.value = scaled(value);
  out.dtheta = scaled(derivative);
  out.d2theta = scaled(second_derivative);
  out.log_abs_value = value == 0.0L
                          ? -std::numeric_limits<real>::infinity()
                          : log_base + recurrence_log_scale +
                                std::log(std::fabs(value));
  if(value!=0.0L) {
    using namespace crli::detail;
    const DoubleDouble total=retained_log_base+
        DoubleDouble{static_cast<double>(recurrence_log_scale),0.0}+
        DoubleDouble{static_cast<double>(std::log(std::fabs(value))),0.0};
    out.log_abs_value=total.hi;out.log_abs_value_lo=total.lo;
  }
  out.sign = (value > 0.0L) - (value < 0.0L);
  out.dtheta_over_value = value == 0.0L
                              ? std::numeric_limits<real>::quiet_NaN()
                              : derivative / value;
  out.d2theta_over_value = value == 0.0L
                               ? std::numeric_limits<real>::quiet_NaN()
                               : second_derivative / value;
  out.residual = 0.0L;
  out.status = status_ok;
  return finite(out.value) && finite(out.dtheta) && finite(out.d2theta);
}

inline bool evaluate_polar_log(real nu, int m, real theta, real scale_log,
                               EvalResult &out) {
  const real kappa = nu + 0.5L;
  const real start = choose_series_start(kappa, m, half_pi);
  const real match = choose_match_angle(kappa, m, start);
  PolarData polar = polar_series(kappa, m, theta);
  if (polar.converged) {
    const real st = std::sin(theta), ct = std::cos(theta);
    const real total_log = polar.log_abs_L + scale_log;
    out.value = static_cast<real>(polar.sign) * safe_exp(total_log);
    out.dtheta = out.value * polar.r_t / st;
    const real coeff = nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
    out.dtheta_over_value = polar.r_t / st;
    out.d2theta_over_value = finite(polar.d2theta_over_value)
                                 ? polar.d2theta_over_value
                                 : -(ct / st) * out.dtheta_over_value - coeff;
    out.d2theta = out.d2theta_over_value * out.value;
    out.log_abs_value = total_log;
    out.sign = polar.sign;
    out.status = status_ok;
    return true;
  }
  if (!(theta < match))
    return false;
  polar = polar_series(kappa, m, start);
  if (!polar.converged)
    return false;
  bool ok = true;
  if (theta <= start)
    polar = polar_series(kappa, m, theta);
  else
    ok = propagate_logarithmic(kappa, m, start, theta, polar);
  if (!ok || !polar.converged)
    return false;
  const real st = std::sin(theta), ct = std::cos(theta),
             total_log = polar.log_abs_L + scale_log;
  out.value = static_cast<real>(polar.sign) * safe_exp(total_log);
  out.dtheta = out.value * polar.r_t / st;
  const real coeff = nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
  out.dtheta_over_value = polar.r_t / st;
  out.d2theta_over_value = finite(polar.d2theta_over_value)
                               ? polar.d2theta_over_value
                               : -(ct / st) * out.dtheta_over_value - coeff;
  out.d2theta = out.d2theta_over_value * out.value;
  out.log_abs_value = total_log;
  out.sign = polar.sign;
  out.status = status_ok;
  return true;
}

// Stable minimal-branch downward recurrence near and below the turning point.
inline bool evaluate_order_downward(real nu, int target, real theta,
                                    EvalResult &out) {
  if (!(theta > 0.0L && theta <= half_pi) || target < 0)
    return false;
  const int top = static_cast<int>(std::floor(nu));
  if (top < target + 1)
    return false;

  struct Pair {
    real y = 0.0L, d = 0.0L;
  };
  auto accurate_scale = [&](int order) {
    // Compensated log-Gamma normalization is O(1) in the bounded-degree chart.
    if (nu < 10000.0L) return legacy_log_scale(nu, order);
    real sum = 0.0L, compensation = 0.0L;
    const real base = nu - static_cast<real>(order);
    for (int k = 1; k <= order; ++k) {
      const real kr = static_cast<real>(k);
      const real term = 0.5L * (std::log(base + kr) + std::log(nu + kr)) -
                        std::log(2.0L * kr);
      const real adjusted = term - compensation;
      const real updated = sum + adjusted;
      compensation = (updated - sum) - adjusted;
      sum = updated;
    }
    return sum + ((order == 0) ? 0.0L : 0.5L * std::log(2.0L));
  };
  auto top_value = [&](int order, Pair &v, real &logmag) {
    const PolarData p = polar_series(nu + 0.5L, order, theta, 5.0e-19L, 20000);
    if (!p.converged)
      return false;
    logmag = p.log_abs_L + accurate_scale(order);
    const real st = std::sin(theta);
    v.y = static_cast<real>(p.sign);
    v.d = v.y * p.r_t / st;
    return finite(logmag) && finite(v.d);
  };

  Pair upper, lower;
  real log_upper = 0.0L, log_lower = 0.0L;
  if (!top_value(top, upper, log_upper) ||
      !top_value(top - 1, lower, log_lower))
    return false;
  real log_scale = std::max(log_upper, log_lower);
  const real eu = safe_exp(log_upper - log_scale);
  const real el = safe_exp(log_lower - log_scale);
  upper.y *= eu;
  upper.d *= eu;
  lower.y *= el;
  lower.d *= el;

  const real st = std::sin(theta), ct = std::cos(theta);
  const real cot = ct / st, csc2 = 1.0L / (st * st);
  for (int j = top - 2; j >= target; --j) {
    const real A0 =
        (nu - static_cast<real>(j)) * (nu + static_cast<real>(j) + 1.0L);
    const real A1 =
        (nu - static_cast<real>(j) - 1.0L) * (nu + static_cast<real>(j) + 2.0L);
    if (!(A0 > 0.0L && A1 > 0.0L))
      return false;
    const real longitude = (j == 0) ? std::sqrt(2.0L) : 1.0L;
    const real rootA1 = std::sqrt(A1);
    const real B = 2.0L * static_cast<real>(j + 1) * cot / rootA1;
    const real Bd = -2.0L * static_cast<real>(j + 1) * csc2 / rootA1;
    const real C = longitude * std::sqrt(A0 / A1);
    Pair next;
    next.y = (B * lower.y - upper.y) / C;
    next.d = (Bd * lower.y + B * lower.d - upper.d) / C;
    upper = lower;
    lower = next;

    const real magnitude = std::max(std::fabs(lower.y), std::fabs(lower.d));
    if (!finite(magnitude) || magnitude == 0.0L)
      return false;
    if (magnitude > 1.0e100L || magnitude < 1.0e-100L) {
      upper.y /= magnitude;
      upper.d /= magnitude;
      lower.y /= magnitude;
      lower.d /= magnitude;
      log_scale += std::log(magnitude);
    }
  }

  if (!finite(lower.y) || !finite(lower.d))
    return false;
  if (lower.y == 0.0L) {
    // A rounded modal zero is valid. Recover its derivative from the
    // adjacent order instead of dividing by the zero value or rejecting
    // the regular solution and attempting phase propagation.
    const real coupling = std::sqrt((nu - static_cast<real>(target)) *
        (nu + static_cast<real>(target) + 1.0L)) /
        (target == 0 ? std::sqrt(2.0L) : 1.0L);
    const real scaled_derivative = -coupling * upper.y;
    if (!finite(scaled_derivative) || scaled_derivative == 0.0L)
      return false;
    out = EvalResult{};
    out.log_abs_dtheta = log_scale + std::log(std::fabs(scaled_derivative));
    out.dtheta_sign = scaled_derivative > 0.0L ? 1 : -1;
    out.dtheta = out.dtheta_sign * safe_exp(out.log_abs_dtheta);
    out.d2theta = -cot * out.dtheta;
    if (cot != 0.0L) {
      out.log_abs_d2theta = out.log_abs_dtheta + std::log(std::fabs(cot));
      out.d2theta_sign = -out.dtheta_sign * (cot > 0.0L ? 1 : -1);
    }
    return finite(out.dtheta) && finite(out.d2theta);
  }
  const real log_value = log_scale + std::log(std::fabs(lower.y));
  const real value = std::copysign(safe_exp(log_value), lower.y);
  // Recover d/dtheta from the exact raising identity after long recurrences.
  const real At = (nu - static_cast<real>(target)) *
                  (nu + static_cast<real>(target) + 1.0L);
  const real longitude_ratio = (target == 0) ? std::sqrt(2.0L) : 1.0L;
  const real derivative_ratio =
      static_cast<real>(target) * cot -
      std::sqrt(At) / longitude_ratio * (upper.y / lower.y);
  const real derivative = value * derivative_ratio;
  const real coeff =
      nu * (nu + 1.0L) - static_cast<real>(target) * target / (st * st);
  out.value = value;
  out.dtheta = derivative;
  out.d2theta = -cot * derivative - coeff * value;
  out.dtheta_over_value = derivative_ratio;
  out.d2theta_over_value = -cot * derivative_ratio - coeff;
  out.log_abs_value = log_value;
  out.sign = (lower.y > 0.0L) ? 1 : -1;
  out.status = status_ok;
  return finite(derivative) || derivative == 0.0L;
}

// Accurate log Gamma(x+a)-log Gamma(x+b) without subtracting two O(x log x)
// values.  This is needed by the near-order normalization on Windows, where
// the Intel compiler's long double has the same precision as double.
inline real shifted_log_gamma_ratio(real x, real a, real b) noexcept {
  if (x < 64.0L)
    return std::lgamma(x + a) - std::lgamma(x + b);
  auto stirling_tail = [](real y) {
    const real inverse = 1.0L / y;
    const real inverse2 = inverse * inverse;
    return inverse * (1.0L / 12.0L + inverse2 *
        (-1.0L / 360.0L + inverse2 * (1.0L / 1260.0L + inverse2 *
        (-1.0L / 1680.0L + inverse2 / 1188.0L))));
  };
  return (a - b) * std::log(x) +
         (x + a - 0.5L) * std::log1p(a / x) -
         (x + b - 0.5L) * std::log1p(b / x) - (a - b) +
         stirling_tail(x + a) - stirling_tail(x + b);
}

inline real stable_near_order_log_scale(real nu, int m) noexcept {
  if (m == 0)
    return legacy_log_scale(nu, m);
  const real gap = nu - static_cast<real>(m);
  const real a = 0.5L * (gap + 1.0L);
  const real mr = static_cast<real>(m);
  return 0.5L * gap * log2_value - 0.25L * std::log(pi) -
         0.5L * std::lgamma(gap + 1.0L) + 0.5L * log2_value +
         0.5L * (shifted_log_gamma_ratio(mr, a, 1.0L) +
                 shifted_log_gamma_ratio(mr, a + 0.5L, 1.0L));
}

struct NearOrderDDState {
  double r = 0.0;
  CRNear::DoubleDouble log_value;
};

struct NearOrderPointKey {
  std::uint64_t degree_bits = 0;
  std::uint64_t theta_bits = 0;
  int order = 0;
  bool operator==(const NearOrderPointKey &other) const noexcept {
    return degree_bits == other.degree_bits && theta_bits == other.theta_bits &&
           order == other.order;
  }
};

struct NearOrderPointKeyHash {
  std::size_t operator()(const NearOrderPointKey &key) const noexcept {
    std::size_t value = std::hash<std::uint64_t>{}(key.degree_bits);
    value ^= std::hash<std::uint64_t>{}(key.theta_bits) +
             0x9e3779b97f4a7c15ULL + (value << 6) + (value >> 2);
    value ^= std::hash<int>{}(key.order) + 0x9e3779b97f4a7c15ULL +
             (value << 6) + (value >> 2);
    return value;
  }
};

inline CRNear::DoubleDouble near_dd_add(CRNear::DoubleDouble a,
                                         const CRNear::DoubleDouble &b) {
  a.add(b);
  return a;
}

// Integrate q=d(log H)/dt after removing the exact sin(theta)^m carrier.
// This avoids cancellation of O(m*log(theta_start)) amplitudes and removes
// the amplification by m of the rounded logarithmic angle coordinates.
inline bool near_order_sdirk_step(double kappa,int m,double t,double h,
    const NearOrderDDState& input,NearOrderDDState& output) noexcept {
  constexpr double gamma=0.292893218813452475599155637895;
  const double degree=kappa-0.5,gap=degree-m;
  const double eigenvalue=gap*(degree+m+1.0);
  auto stage=[&](double base,double time,double& result) {
    const double tangent=std::tanh(time),sech2=1.0/(std::cosh(time)*std::cosh(time));
    const double diagonal=gamma*h,b=1.0-2.0*m*diagonal*tangent;
    const double shifted=base-diagonal*eigenvalue*sech2;
    const double discriminant=b*b+4.0*diagonal*shifted;
    if(!(discriminant>0.0))return false;
    result=2.0*shifted/(b+std::sqrt(discriminant));return std::isfinite(result);
  };
  double first=0,second=0;
  if(!stage(input.r,t+gamma*h,first))return false;
  const double slope=(first-input.r)/(gamma*h);
  if(!stage(input.r+h*(1.0-gamma)*slope,t+h,second))return false;
  output.r=second;
  output.log_value=near_dd_add(input.log_value,
      CRNear::DoubleDouble(h*((1.0-gamma)*first+gamma*second)));
  return std::isfinite(output.log_value.value());
}

inline bool near_order_fixed_mesh(double kappa, int m, double t0, double t1,
                                  int count, const NearOrderDDState &initial,
                                  NearOrderDDState &result) noexcept {
  result = initial;
  const double h = (t1 - t0) / static_cast<double>(count);
  for (int index = 0; index < count; ++index) {
    const double t = t0 + static_cast<double>(index) * h;
    NearOrderDDState next;
    if (!near_order_sdirk_step(kappa, m, t, h, result, next))
      return false;
    result = next;
  }
  return true;
}

inline bool evaluate_fractional_near_order_evanescent(
    real nu, int m, real theta, EvalResult &out) noexcept {
  const real gap = nu - static_cast<real>(m);
  const real kappa = nu + 0.5L;
  if (!(nu >= 512.0L && gap >= 0.0L && gap <= 64.0L &&
        m > 0 && static_cast<real>(m) / kappa > 0.99L &&
        theta >= 85.0L * pi / 180.0L && theta <= half_pi))
    return false;
  const real eta2 = static_cast<real>(m) * m - 0.25L;
  if (!(kappa * kappa - eta2 / sqr(std::sin(theta)) < 0.0L))
    return false;
  const auto local_map = crli::detail::canonical_map(
      static_cast<double>(static_cast<real>(m) / kappa),
      static_cast<double>(theta));
  if (static_cast<double>(kappa) * std::abs(local_map.action) < 8.0)
    return false;

  const double degree_key = static_cast<double>(nu);
  const double theta_key = static_cast<double>(theta);
  NearOrderPointKey key;
  std::memcpy(&key.degree_bits, &degree_key, sizeof(degree_key));
  std::memcpy(&key.theta_bits, &theta_key, sizeof(theta_key));
  key.order = m;
  thread_local std::unordered_map<NearOrderPointKey, EvalResult,
                                  NearOrderPointKeyHash> cache;
  const auto cached = cache.find(key);
  if (cached != cache.end()) {
    out = cached->second;
    return true;
  }

  const real start=choose_series_start(kappa,m,half_pi);
  if(!(start<theta))return false;
  const real z=sqr(std::sin(start/2));
  real term=1,sum=1,weighted=0;bool initialized=false;
  for(int j=1;j<=128;++j) {
    term*=((j-1.0L-gap)*(j+m+nu)/((j+m)*static_cast<real>(j)))*z;
    sum+=term;weighted+=j*term;
    if(std::fabs(term)*j<2.0e-19L*std::fabs(sum)){initialized=true;break;}
  }
  if(!initialized || !(sum>0))return false;
  const double t0=std::log(std::tan(0.5*static_cast<double>(start)));
  const double t1=std::log(std::tan(0.5*static_cast<double>(theta)));
  const NearOrderDDState initial{
      static_cast<double>(2*sqr(std::cos(start/2))*weighted/sum),
      CRNear::DoubleDouble(static_cast<double>(std::log(sum)))};
  int mesh = std::max(2048, static_cast<int>(std::ceil(
      static_cast<double>(kappa) / 200.0)));
  NearOrderDDState coarse;
  if (!near_order_fixed_mesh(static_cast<double>(kappa), m, t0, t1,
                             mesh, initial, coarse))
    return false;

  bool converged = false;
  double previous_r=0.0,extrapolated_r=0.0,log_error=0.0;
  crli::detail::DoubleDouble previous_log{0.0,0.0},extrapolated_log{0.0,0.0};
  for(int refinement=0;refinement<12;++refinement) {
    if (mesh > 25000000)
      return false;
    NearOrderDDState fine;
    if (!near_order_fixed_mesh(static_cast<double>(kappa), m, t0, t1,
                               2 * mesh, initial, fine))
      return false;
    extrapolated_r = fine.r + (fine.r - coarse.r) / 3.0;
    const crli::detail::DoubleDouble fine_log{fine.log_value.hi,fine.log_value.lo},
        coarse_log{coarse.log_value.hi,coarse.log_value.lo};
    extrapolated_log=fine_log+(fine_log-coarse_log)/crli::detail::DoubleDouble{3.0,0.0};
    const auto change=extrapolated_log-previous_log;
    log_error=std::abs(change.hi+change.lo);
    if(refinement>0 && log_error<=1.0e-9 &&
        std::abs(extrapolated_r-previous_r)<=2.0e-10*(1.0+std::abs(extrapolated_r))) {
      converged = true;
      break;
    }
    previous_r = extrapolated_r;
    previous_log = extrapolated_log;
    coarse = fine;
    mesh *= 2;
  }
  if (!converged)
    return false;

  const real st = std::sin(theta), ct = std::cos(theta);
  using namespace crli::detail;
  const DoubleDouble c{static_cast<double>(ct),0.0};
  const auto log_sine=dd_log_positive(DoubleDouble{1.0,0.0}-c*c)*0.5;
  const auto logarithm=log_sine*static_cast<double>(m)+extrapolated_log+
      DoubleDouble{static_cast<double>(stable_near_order_log_scale(nu,m)),0.0};
  out.log_abs_value=logarithm.hi;out.log_abs_value_lo=logarithm.lo;
  out.sign=1;out.value=safe_exp(out.log_abs_value+out.log_abs_value_lo);
  out.dtheta_over_value=(m*ct+static_cast<real>(extrapolated_r))/st;
  const real coefficient=gap*(nu+m+1.0L)+m-static_cast<real>(m)*m*sqr(ct/st);
  out.d2theta_over_value = -(ct / st) * out.dtheta_over_value - coefficient;
  out.dtheta = out.value * out.dtheta_over_value;
  out.d2theta = out.value * out.d2theta_over_value;
  out.residual=static_cast<real>(log_error);
  out.status = status_ok;
  const bool valid = finite(out.log_abs_value) &&
                     finite(out.dtheta_over_value) &&
                     finite(out.d2theta_over_value);
  if (valid) {
    if (cache.size() >= 512U)
      cache.clear();
    cache.emplace(key, out);
  }
  return valid;
}

// Use exact equatorial data across the ill-conditioned turning transition.
inline bool evaluate_coalescent_wkb(real nu,int m,real theta,EvalResult& out) {
  const real kappa=nu+0.5L,gap=nu-m;
  if(!(kappa>=1000000.0L && gap>=64.0L &&
       theta>=85.0L*pi/180.0L && theta<=half_pi && m/kappa>0.95L))return false;
  const real span=(half_pi-theta)*std::sqrt((kappa-m)*(kappa+m)+0.25L);
  const auto pt=wkb_density(kappa,m,theta,true),p0=wkb_density(kappa,m,half_pi,true);
  if(!pt.valid || !p0.valid ||
      std::max(pt.relative_defect,p0.relative_defect)*std::max(1.0L,span)>1.0e-8L)return false;
  // Integrate the physical density from the exact Gamma anchor. A phase
  // difference avoids cancellation of two large canonical actions near
  // the turning-point/equator coalescence.
  const auto density=[=](real x)->real {
    const auto p=wkb_density(kappa,m,x,true);
    return p.valid && p.relative_defect*std::max(1.0L,span)<=1.0e-8L
        ? p.p : std::numeric_limits<real>::quiet_NaN();
  };
  real phase=0;
  if(!gauss_phase_integral(density,theta,half_pi,phase,5.0e-10L))return false;
  phase=-phase-p0.p*6.12323399573676603587e-17L;
  real y0=0,d0=0;equator_legacy_values(nu,m,y0,d0);
  const real amplitude=std::sqrt(p0.p/pt.p),s=std::sin(theta),cot=std::cos(theta)/s;
  const real sn=std::sin(phase),cs=std::cos(phase);
  const real shape=y0*cs+d0/p0.p*sn;
  const real y=amplitude*shape;
  const real dy=amplitude*(pt.p*(-y0*sn+d0/p0.p*cs)-0.5L*pt.ptheta/pt.p*shape);
  out.value=y/std::sqrt(s);out.dtheta=(dy-0.5L*cot*y)/std::sqrt(s);
  const real coefficient=(nu-m)*(nu+m)+nu-static_cast<real>(m)*m*cot*cot;
  out.d2theta=-cot*out.dtheta-coefficient*out.value;
  out.status=status_ok;
  return finite(out.value)&&finite(out.dtheta)&&finite(out.d2theta);
}

inline bool evaluate_equator_ivp(real nu, int m, real theta,
                                 EvalResult &out) noexcept {
  if (m < 1 || !(theta > 0.75L * half_pi && theta <= half_pi))
    return false;
  const real gap = nu - static_cast<real>(m);
  if (!(gap >= 0.0L))
    return false;
  const real kappa = nu + 0.5L;
  const real eta2 = static_cast<real>(m) * m - 0.25L;
  const real local_cot = std::cos(theta)/std::sin(theta);
  const real qtheta = (kappa-m)*(kappa+m)+0.25L-eta2*local_cot*local_cot;
  real scaled_forbidden_action = 0.0L;
  if (qtheta < 0.0L) {
    const auto local_map = crli::detail::canonical_map(
        static_cast<double>(static_cast<real>(m) / kappa),
        static_cast<double>(theta));
    scaled_forbidden_action =
        kappa * static_cast<real>(std::abs(local_map.action));
    if (scaled_forbidden_action > 8.0L)
      return false;
  }
  const real equator_wavenumber = std::sqrt(
      std::max<real>(0.0L, nu * (nu + 1.0L) - static_cast<real>(m) * m));
  const real phase_span = (half_pi - theta) * equator_wavenumber;
  const bool original_near_order = gap < 64.0L;
  const bool near_turning_equator =
      qtheta > 0.0L && qtheta < 0.04L * kappa * kappa && phase_span < 64.0L;
  const bool accurate_high_equator = kappa >= 1000000.0L &&
      theta >= 89.9L * pi / 180.0L && qtheta > 0.0L &&
      qtheta < 0.04L * kappa * kappa && phase_span <= 512.0L;
  const bool resolved_high_near_order=kappa>=1000000.0L && gap<=4096.0L &&
      qtheta>0.0L && phase_span<=512.0L;
  if (!(original_near_order || near_turning_equator || accurate_high_equator || resolved_high_near_order))
    return false;

  real y0 = 0.0L, d0 = 0.0L;
  equator_legacy_values(nu, m, y0, d0);
  if (!finite(y0) || !finite(d0))
    return false;
  // Integrate in carrier-scaled distance from the representable center.
  // This avoids accumulating ulp-sized angular steps at very high degree.
  constexpr real center_delta = -6.1232339957367660359e-17L;
  const real center_coefficient = (nu-m)*(nu+m)+nu;
  const real integration_scale = gap <= 20000.0L ? std::sqrt(kappa) : kappa;
  const real old_y0 = y0;
  y0 += center_delta * d0;
  d0 -= center_delta * center_coefficient * old_y0;
  std::array<real, 2> state{y0, d0 / integration_scale};
  auto rhs = [=](real x, const std::array<real, 2> &y, std::array<real, 2> &f) {
    const real offset = center_delta + x / integration_scale;
    const real st = std::cos(offset), ct = -std::sin(offset);
    const real cot=ct/st;
    const real coeff = center_coefficient-static_cast<real>(m)*m*cot*cot;
    f[0] = y[1];
    f[1] = -(ct / st) / integration_scale * y[1] - coeff / (integration_scale * integration_scale) * y[0];
  };
  const real target = integration_scale * (theta - half_pi);
  if (!dopri5_integrate(0.0L, target, state, rhs, 2.0e-15L, 2.0e-18L, 0.02L)) {
    state = {y0, d0 / integration_scale};
    if (!(near_turning_equator || accurate_high_equator || resolved_high_near_order ||
          (qtheta < 0.0L && scaled_forbidden_action <= 8.0L)) ||
        !dopri5_integrate(0.0L, target, state, rhs, 2.0e-12L, 2.0e-15L, 0.02L))
      return false;
  }
  state[1] *= integration_scale;
  const real st = std::sin(theta), ct = std::cos(theta);
  const real coeff = nu * (nu + 1.0L) -
                     static_cast<real>(m) * static_cast<real>(m) / (st * st);
  out.value = state[0];
  out.dtheta = state[1];
  out.d2theta = -(ct / st) * state[1] - coeff * state[0];
  out.log_abs_value = (state[0] == 0.0L)
                          ? -std::numeric_limits<real>::infinity()
                          : std::log(std::fabs(state[0]));
  out.sign = (state[0] > 0.0L) - (state[0] < 0.0L);
  out.status = status_ok;
  return finite(out.value) && finite(out.dtheta) && finite(out.d2theta);
}

struct CoreState {
  real kappa = 0.0L;
  int m = 0;
  real theta_start = 0.0L;
  real theta_match = 0.0L;
  PolarData polar_start;
  PolarData polar_match;
  PhaseState phase_match;
  MatchData match;
  int status = status_ok;
};

inline CoreState prepare_core(real nu, int m, real target) {
  CoreState core;
  core.kappa = nu + 0.5L;
  core.m = m;
  if (!(nu >= static_cast<real>(m)) || m < 0 ||
      !(target > 0.0L && target < pi)) {
    core.status = status_domain;
    return core;
  }

  real start = choose_series_start(core.kappa, m, target);
  PolarData polar;
  bool got = false;
  for (int attempt = 0; attempt < 12; ++attempt) {
    polar = polar_series(core.kappa, m, start);
    if (polar.converged) {
      got = true;
      break;
    }
    start *= 0.5L;
    if (start < 1.0e-12L)
      break;
  }
  if (!got) {
    core.status = status_series;
    return core;
  }

  const real match_theta = choose_match_angle(core.kappa, m, start);
  const PolarData polar_start = polar;
  if (match_theta > start) {
    if (!propagate_logarithmic(core.kappa, m, start, match_theta, polar)) {
      core.status = status_integrator;
      return core;
    }
  }

  PhaseState phase = phase_anchor(core.kappa, m);
  if (!integrate_phase(core.kappa, m, half_pi, match_theta, phase)) {
    core.status = status_integrator;
    return core;
  }

  MatchData match = build_match(match_theta, polar, phase);
  if (match.status != status_ok) {
    core.status = match.status;
    return core;
  }

  core.theta_start = start;
  core.theta_match = match_theta;
  core.polar_start = polar_start;
  core.polar_match = polar;
  core.phase_match = phase;
  core.match = match;
  core.status = status_ok;
  return core;
}

inline bool phase_state_at(real kappa, int m, real theta, PhaseState &s) {
  s = phase_anchor(kappa, m);
  return integrate_phase(kappa, m, half_pi, theta, s);
}

struct ModeKey {
  std::uint64_t nu_bits = 0;
  int m = 0;
  bool operator==(const ModeKey &other) const noexcept {
    return nu_bits == other.nu_bits && m == other.m;
  }
};

struct ModeKeyHash {
  std::size_t operator()(const ModeKey &key) const noexcept {
    const std::size_t h1 = std::hash<std::uint64_t>{}(key.nu_bits);
    const std::size_t h2 = std::hash<int>{}(key.m);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

inline ModeKey make_mode_key(real nu, int m) noexcept {
  const double nud = static_cast<double>(nu);
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(nud), "unexpected double size");
  std::memcpy(&bits, &nud, sizeof(bits));
  return ModeKey{bits, m};
}

inline CoreState prepare_core_cached(real nu, int m) {

  thread_local std::unordered_map<ModeKey, CoreState, ModeKeyHash> cache;
  const ModeKey key = make_mode_key(nu, m);
  const auto found = cache.find(key);
  if (found != cache.end())
    return found->second;

  CoreState core = prepare_core(nu, m, half_pi);
  // Preparation depends only on (nu, m), including deterministic failures.
  // Preserve failed states so repeated calls take the existing rescue paths
  // without retrying the same unsuccessful integration. The first attempt and
  // the rescue algorithms are unchanged; the thread-local cache remains bounded.
  if (cache.size() >= 8192U)
    cache.clear();
  cache.emplace(key, core);
  return core;
}

inline EvalResult evaluate_legacy(real nu, int m, real theta) {
  EvalResult out;
  if (m < 0 || !(nu >= static_cast<real>(m)) ||
      !(theta >= 0.0L && theta < pi)) {
    out.status = status_domain;
    out.value = out.dtheta = out.d2theta =
        std::numeric_limits<real>::quiet_NaN();
    return out;
  }

  const real kappa = nu + 0.5L;
  const real scale_log = legacy_log_scale(nu, m);

  if (theta == 0.0L) {
    if (m == 0) {
      out.value = 1.0L;
      out.dtheta = 0.0L;
      out.d2theta = -0.5L * nu * (nu + 1.0L);
    } else if (m == 1) {
      const real amp = safe_exp(scale_log);
      out.value = 0.0L;
      out.dtheta = amp;
      out.d2theta = 0.0L;
    } else if (m == 2) {
      const real amp = safe_exp(scale_log);
      out.value = 0.0L;
      out.dtheta = 0.0L;
      out.d2theta = 2.0L * amp;
    } else {
      out.value = out.dtheta = out.d2theta = 0.0L;
    }
    out.sign = (out.value > 0.0L) - (out.value < 0.0L);
    out.log_abs_value = (out.value == 0.0L)
                            ? -std::numeric_limits<real>::infinity()
                            : std::log(std::fabs(out.value));
    return out;
  }

  const real mr = static_cast<real>(m);
  const real turning_sine =
      (m == 0)
          ? 0.0L
          : std::min<real>(
                1.0L, std::sqrt(std::max<real>(0.0L, mr * mr - 0.25L)) / kappa);
  const real turning_angle = std::asin(turning_sine);
  EvalResult equator;
  if (theta >= turning_angle && evaluate_equator_ivp(nu, m, theta, equator))
    return equator;
  if ((kappa <= 512.0L || theta <= turning_angle ||
       (kappa < 10000.0L && theta <= pi/6.0L)) &&
      evaluate_order_downward(nu, m, theta, out))
    return out;

  if (evaluate_polar_log(nu, m, theta, scale_log, out))
    return out;

  // Skip cancellation-dominated O(n) series when the gated LG path applies.
  if (kappa * theta > 24.0L && wkb_evaluate_legacy(nu, m, theta, out))
    return out;

  PolarData direct = polar_series(kappa, m, theta);
  if (direct.converged) {
    const real st_direct = std::sin(theta);
    const real ct_direct = std::cos(theta);
    const real total_log = direct.log_abs_L + scale_log;
    const real mag = safe_exp(total_log);
    out.value = static_cast<real>(direct.sign) * mag;
    out.dtheta = out.value * direct.r_t / st_direct;
    const real coeff =
        nu * (nu + 1.0L) - static_cast<real>(m) * m / (st_direct * st_direct);
    out.dtheta_over_value = direct.r_t / st_direct;
    out.d2theta_over_value = finite(direct.d2theta_over_value)
                                 ? direct.d2theta_over_value
                                 : -(ct_direct / st_direct) *
                                       out.dtheta_over_value - coeff;
    out.d2theta = out.d2theta_over_value * out.value;
    out.sign = (out.value > 0.0L) - (out.value < 0.0L);
    out.log_abs_value = (out.value == 0.0L)
                            ? -std::numeric_limits<real>::infinity()
                            : std::log(std::fabs(out.value));
    out.status = status_ok;
    return out;
  }

  // The zero-free polar branch is best propagated in logarithmic form.
  const real local_start = choose_series_start(kappa, m, half_pi);
  const real local_match = choose_match_angle(kappa, m, local_start);
  if (theta < local_match) {
    PolarData polar = polar_series(kappa, m, local_start);
    if (polar.converged) {
      bool ok = true;
      if (theta <= local_start)
        polar = polar_series(kappa, m, theta);
      else
        ok = propagate_logarithmic(kappa, m, local_start, theta, polar);
      if (ok && polar.converged) {
        const real st = std::sin(theta), ct = std::cos(theta);
        const real total_log = polar.log_abs_L + scale_log;
        out.value = static_cast<real>(polar.sign) * safe_exp(total_log);
        out.dtheta = out.value * polar.r_t / st;
        const real coeff =
            nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
        out.dtheta_over_value = polar.r_t / st;
        out.d2theta_over_value = finite(polar.d2theta_over_value)
                                     ? polar.d2theta_over_value
                                     : -(ct / st) * out.dtheta_over_value - coeff;
        out.d2theta = out.d2theta_over_value * out.value;
        out.log_abs_value = total_log;
        out.sign = polar.sign;
        out.status = status_ok;
        return out;
      }
    }
  }

  // Accept the table-free oscillatory backend only after its defect check.
  if (wkb_evaluate_legacy(nu, m, theta, out))
    return out;

  // Bounded degree recurrence avoids an unnecessary phase IVP at regular
  // points where the accepted WKB approximation is unavailable. The same
  // recurrence already served as the rescue after a failed phase IVP.
  if (nu < 10000.0L && evaluate_order_downward(nu, m, theta, out))
    return out;

  CoreState core = prepare_core_cached(nu, m);
  if (core.status != status_ok) {
    EvalResult equator_rescue;
    if (evaluate_equator_ivp(nu, m, theta, equator_rescue))
      return equator_rescue;
    // Last-resort evaluator when phase and WKB paths reject a regular point.
    EvalResult recurrence_rescue;
    if (evaluate_order_downward(nu, m, theta, recurrence_rescue))
      return recurrence_rescue;
    out.status = core.status;
    return out;
  }

  const real st = std::sin(theta);
  const real ct = std::cos(theta);
  const real cot = ct / st;

  if (theta < core.theta_match - 1.0e-18L) {
    PolarData polar = core.polar_start;
    if (theta <= core.theta_start) {
      polar = polar_series(kappa, m, theta);
      if (!polar.converged) {
        out.status = status_series;
        return out;
      }
    } else if (!propagate_logarithmic(kappa, m, core.theta_start, theta,
                                      polar)) {
      out.status = status_integrator;
      return out;
    }
    const real total_log = polar.log_abs_L + scale_log;
    const real mag = safe_exp(total_log);
    out.value = static_cast<real>(polar.sign) * mag;
    out.dtheta = out.value * polar.r_t / st;
    out.dtheta_over_value = polar.r_t / st;
    out.phase = 0.0L;
    out.phase_derivative = 0.0L;
    out.phase_kappa = 0.0L;
  } else {
    PhaseState phase;
    if (!phase_state_at(kappa, m, theta, phase)) {
      out.status = status_integrator;
      return out;
    }
    const real Phi = phase.alpha + core.match.delta;
    const real Phi_k = phase.alpha_kappa + core.match.delta_kappa;
    const real pphase = safe_exp(phase.q);
    const real total_log =
        core.match.log_C - 0.5L * phase.q - 0.5L * std::log(st) + scale_log;
    const real amp = safe_exp(total_log);
    const real sn = std::sin(Phi);
    const real cs = std::cos(Phi);
    out.value = amp * sn;
    out.dtheta = amp * (pphase * cs + (phase.a - 0.5L * cot) * sn);
    out.phase = Phi;
    out.phase_derivative = pphase;
    out.phase_kappa = Phi_k;
  }

  const real coeff = nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
  out.d2theta = -cot * out.dtheta - coeff * out.value;
  if (finite(out.dtheta_over_value))
    out.d2theta_over_value = -cot * out.dtheta_over_value - coeff;
  out.sign = (out.value > 0.0L) - (out.value < 0.0L);
  out.log_abs_value = (out.value == 0.0L)
                          ? -std::numeric_limits<real>::infinity()
                          : std::log(std::fabs(out.value));
  const real R = out.d2theta + cot * out.dtheta + coeff * out.value;
  const real denom = std::fabs(out.d2theta) + std::fabs(cot * out.dtheta) +
                     std::fabs(coeff * out.value) +
                     std::numeric_limits<real>::min();
  out.residual = std::fabs(R) / denom;
  out.status =
      (finite(out.value) || out.value == 0.0L) ? status_ok : status_range;
  return out;
}

inline bool evaluate_legacy_degree_consistent(real nu, int m, real theta,
                                              EvalResult &out) {
  if (!(nu >= static_cast<real>(m + 2)))
    return false;
  EvalResult center = evaluate_legacy(nu, m, theta);
  EvalResult adjacent = evaluate_legacy(nu - 1.0L, m, theta);
  EvalResult previous = evaluate_legacy(nu - 2.0L, m, theta);
  if (center.status != status_ok || adjacent.status != status_ok ||
      previous.status != status_ok || center.sign == 0 ||
      adjacent.sign == 0 || previous.sign == 0 ||
      !finite(center.log_abs_value) || !finite(adjacent.log_abs_value) ||
      !finite(previous.log_abs_value))
    return false;
  const real st = std::sin(theta), ct = std::cos(theta);
  const real coupling = std::sqrt(
      (nu - static_cast<real>(m)) * (nu + static_cast<real>(m)));
  const real adjacent_ratio = static_cast<real>(adjacent.sign * center.sign) *
      safe_exp(adjacent.log_abs_value - center.log_abs_value);
  const real adjacent_coupling = std::sqrt(
      (nu - 1.0L - static_cast<real>(m)) *
      (nu - 1.0L + static_cast<real>(m)));
  const real previous_ratio = static_cast<real>(previous.sign * adjacent.sign) *
      safe_exp((previous.log_abs_value - adjacent.log_abs_value)+
              (previous.log_abs_value_lo-adjacent.log_abs_value_lo));
  const real adjacent_derivative =
      ((nu - 1.0L) * ct - adjacent_coupling * previous_ratio) / st;
  out = center;
  out.dtheta_over_value = (nu * ct - coupling * adjacent_ratio) / st;
  out.d2theta_over_value =
      ((nu - 1.0L) * ct * out.dtheta_over_value - nu * st -
       coupling * adjacent_ratio * adjacent_derivative) / st;
  out.dtheta = out.value * out.dtheta_over_value;
  out.d2theta = out.value * out.d2theta_over_value;
  return finite(out.dtheta_over_value) && finite(out.d2theta_over_value);
}

struct DegreeSensitivityData {
  EvalResult center;
  real value_nu = 0.0L;
  real dtheta_nu = 0.0L;
  real log_amplitude = 0.0L;
  real scaled_value = 0.0L;
  real scaled_dtheta = 0.0L;
  real scaled_value_nu = 0.0L;
  real scaled_dtheta_nu = 0.0L;
  // A scaled record may be valid even when ordinary values are out of range.
  bool valid = false;
  bool materialized = false;
  int failure_stage = 0;
};

// Differentiate the stable recurrence to keep value and sensitivity on one
// branch.
inline DegreeSensitivityData
degree_sensitivity_order_downward(real nu, int target, real theta) {
  DegreeSensitivityData out;
  if (!(theta > 0.0L && theta <= half_pi) || target < 0)
    return out;
  const int top = static_cast<int>(std::floor(nu));
  if (top < target + 1)
    return out;

  struct Pair {
    real y = 0.0L;
    real y_nu = 0.0L;
  };
  auto accurate_scale = [&](int order) {
    real sum = 0.0L, compensation = 0.0L;
    const real base = nu - static_cast<real>(order);
    for (int k = 1; k <= order; ++k) {
      const real kr = static_cast<real>(k);
      const real term = 0.5L * (std::log(base + kr) + std::log(nu + kr)) -
                        std::log(2.0L * kr);
      const real adjusted = term - compensation;
      const real updated = sum + adjusted;
      compensation = (updated - sum) - adjusted;
      sum = updated;
    }
    return sum + ((order == 0) ? 0.0L : 0.5L * std::log(2.0L));
  };
  auto scale_nu = [&](int order) {
    return 0.5L * (digamma_positive(nu + static_cast<real>(order) + 1.0L) -
                   digamma_positive(nu - static_cast<real>(order) + 1.0L));
  };
  auto top_value = [&](int order, Pair &value, real &log_magnitude,
                       real &log_magnitude_nu) {
    const PolarData polar =
        polar_series(nu + 0.5L, order, theta, 5.0e-19L, 20000);
    if (!polar.converged)
      return false;
    log_magnitude = polar.log_abs_L + accurate_scale(order);
    log_magnitude_nu = polar.log_abs_L_kappa + scale_nu(order);
    value.y = static_cast<real>(polar.sign);
    value.y_nu = 0.0L;
    return finite(log_magnitude) && finite(log_magnitude_nu);
  };

  Pair upper, lower;
  real log_upper = 0.0L, log_lower = 0.0L;
  real log_upper_nu = 0.0L, log_lower_nu = 0.0L;
  if (!top_value(top, upper, log_upper, log_upper_nu) ||
      !top_value(top - 1, lower, log_lower, log_lower_nu))
    return out;
  real log_scale = std::max(log_upper, log_lower);
  real log_scale_nu = (log_upper >= log_lower) ? log_upper_nu : log_lower_nu;
  const real upper_factor = safe_exp(log_upper - log_scale);
  const real lower_factor = safe_exp(log_lower - log_scale);
  upper.y *= upper_factor;
  upper.y_nu = upper.y * (log_upper_nu - log_scale_nu);
  lower.y *= lower_factor;
  lower.y_nu = lower.y * (log_lower_nu - log_scale_nu);

  const real st = std::sin(theta), ct = std::cos(theta);
  const real cot = ct / st;
  for (int order = top - 2; order >= target; --order) {
    const real A0 = (nu - static_cast<real>(order)) *
                    (nu + static_cast<real>(order) + 1.0L);
    const real A1 = (nu - static_cast<real>(order) - 1.0L) *
                    (nu + static_cast<real>(order) + 2.0L);
    if (!(A0 > 0.0L && A1 > 0.0L))
      return out;
    const real derivative_A = 2.0L * nu + 1.0L;
    const real longitude = (order == 0) ? std::sqrt(2.0L) : 1.0L;
    const real B = 2.0L * static_cast<real>(order + 1) * cot / std::sqrt(A1);
    const real B_nu = -0.5L * B * derivative_A / A1;
    const real C = longitude * std::sqrt(A0 / A1);
    const real C_nu_over_C = 0.5L * derivative_A * (1.0L / A0 - 1.0L / A1);
    Pair next;
    next.y = (B * lower.y - upper.y) / C;
    next.y_nu = (B_nu * lower.y + B * lower.y_nu - upper.y_nu) / C -
                next.y * C_nu_over_C;
    upper = lower;
    lower = next;

    const real magnitude = std::max(std::fabs(lower.y), std::fabs(upper.y));
    if (!finite(magnitude) || magnitude == 0.0L)
      return out;
    if (magnitude > 1.0e100L || magnitude < 1.0e-100L) {
      upper.y /= magnitude;
      upper.y_nu /= magnitude;
      lower.y /= magnitude;
      lower.y_nu /= magnitude;
      log_scale += std::log(magnitude);
    }
  }

  if (!finite(lower.y) || !finite(lower.y_nu) || !finite(upper.y) ||
      !finite(upper.y_nu))
    return out;
  const real A = (nu - static_cast<real>(target)) *
                 (nu + static_cast<real>(target) + 1.0L);
  if (!(A > 0.0L))
    return out;
  const real root_A = std::sqrt(A);
  const real root_A_nu = (2.0L * nu + 1.0L) / (2.0L * root_A);
  const real longitude_ratio = (target == 0) ? std::sqrt(2.0L) : 1.0L;
  out.log_amplitude = log_scale;
  out.scaled_value = lower.y;
  out.scaled_value_nu = lower.y_nu + lower.y * log_scale_nu;
  out.scaled_dtheta = static_cast<real>(target) * cot * lower.y -
                     root_A * upper.y / longitude_ratio;
  out.scaled_dtheta_nu = static_cast<real>(target) * cot * out.scaled_value_nu -
      (root_A_nu * upper.y + root_A * (upper.y_nu + upper.y * log_scale_nu)) /
          longitude_ratio;
  const auto materialize = [log_scale](real x) {
    return x == 0.0L ? 0.0L : std::copysign(safe_exp(log_scale + std::log(std::fabs(x))), x);
  };
  const real value = materialize(out.scaled_value);
  const real derivative = materialize(out.scaled_dtheta);
  const real value_nu = materialize(out.scaled_value_nu);
  const real derivative_nu = materialize(out.scaled_dtheta_nu);

  const real coefficient =
      nu * (nu + 1.0L) - static_cast<real>(target) * target / (st * st);
  out.center.value = value;
  out.center.dtheta = derivative;
  out.center.d2theta = -cot * derivative - coefficient * value;
  out.center.status = status_ok;
  out.value_nu = value_nu;
  out.dtheta_nu = derivative_nu;
  out.valid = finite(out.log_amplitude) && finite(out.scaled_value) &&
              finite(out.scaled_dtheta) && finite(out.scaled_value_nu) &&
              finite(out.scaled_dtheta_nu);
  out.materialized = finite(value) && finite(derivative) && finite(value_nu) &&
                     finite(derivative_nu) &&
                     (out.scaled_value == 0.0L || value != 0.0L) &&
                     (out.scaled_dtheta == 0.0L || derivative != 0.0L) &&
                     (out.scaled_value_nu == 0.0L || value_nu != 0.0L) &&
                     (out.scaled_dtheta_nu == 0.0L || derivative_nu != 0.0L);
  return out;
}

// Seed in the forbidden region, then continue the same normalized branch
// and its exact degree sensitivity through the turning point.
inline DegreeSensitivityData degree_sensitivity_turning_ivp(real nu, int m,
                                                            real theta) {
  DegreeSensitivityData out;
  if (m < 1 || !(nu >= static_cast<real>(m)) ||
      !(theta > 0.0L && theta <= half_pi))
    return out;

  const real kappa = nu + 0.5L;
  const real eta2 = static_cast<real>(m) * m - 0.25L;
  if (!(eta2 > 0.0L) || !(eta2 < kappa * kappa))
    return out;
  const real turning = std::asin(std::sqrt(eta2) / kappa);

  // Stay several Airy widths from turning without underflowing the seed.
  const real margin = std::max<real>(0.02L, 12.0L / kappa);
  const real match = std::max<real>(1.0e-8L, turning - margin);
  if (theta <= match)
    return degree_sensitivity_order_downward(nu, m, theta);

  DegreeSensitivityData seed = degree_sensitivity_order_downward(nu, m, match);
  if (!seed.valid || !seed.materialized)
    return out;

  std::array<real, 4> state{seed.center.value, seed.center.dtheta,
                            seed.value_nu, seed.dtheta_nu};
  real magnitude = std::max({std::fabs(state[0]), std::fabs(state[1]),
                             std::fabs(state[2]), std::fabs(state[3])});
  if (!(magnitude > 0.0L) || !finite(magnitude))
    return out;
  for (real &component : state)
    component /= magnitude;
  real log_amplitude = std::log(magnitude);

  auto rhs = [=](real x, const std::array<real, 4> &y, std::array<real, 4> &f) {
    const real st = std::sin(x);
    const real cot = std::cos(x) / st;
    const real coefficient =
        nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
    f[0] = y[1];
    f[1] = -cot * y[1] - coefficient * y[0];
    f[2] = y[3];
    f[3] = -cot * y[3] - coefficient * y[2] - (2.0L * nu + 1.0L) * y[0];
  };

  real current = match;
  while (current < theta) {
    const real next = std::min(theta, current + 0.002L);
    if (!dopri5_integrate(current, next, state, rhs, 8.0e-16L, 2.0e-18L,
                          0.0002L)) {
      out.failure_stage = 4;
      return out;
    }
    magnitude = std::max({std::fabs(state[0]), std::fabs(state[1]),
                          std::fabs(state[2]), std::fabs(state[3])});
    if (!(magnitude > 0.0L) || !finite(magnitude)) {
      out.failure_stage = 4;
      return out;
    }
    for (real &component : state)
      component /= magnitude;
    log_amplitude += std::log(magnitude);
    current = next;
  }

  out.log_amplitude = log_amplitude;
  out.scaled_value = state[0];
  out.scaled_dtheta = state[1];
  out.scaled_value_nu = state[2];
  out.scaled_dtheta_nu = state[3];
  out.valid = finite(log_amplitude) && finite(state[0]) && finite(state[1]) &&
              finite(state[2]) && finite(state[3]);
  if (!out.valid)
    return out;

  const real amplitude = safe_exp(log_amplitude);
  if (!(amplitude > 0.0L) || !finite(amplitude)) {
    out.center.status = status_range;
    return out;
  }
  out.center.value = amplitude * state[0];
  out.center.dtheta = amplitude * state[1];
  out.value_nu = amplitude * state[2];
  out.dtheta_nu = amplitude * state[3];
  const real st = std::sin(theta);
  const real cot = std::cos(theta) / st;
  const real coefficient =
      nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
  out.center.d2theta =
      -cot * out.center.dtheta - coefficient * out.center.value;
  out.center.status = status_ok;
  out.materialized = finite(out.center.value) && finite(out.center.dtheta) &&
                     finite(out.value_nu) && finite(out.dtheta_nu);
  return out;
}

inline real legacy_log_scale_nu(real nu, int m) {
  return 0.5L * (digamma_positive(nu + static_cast<real>(m) + 1.0L) -
                 digamma_positive(nu - static_cast<real>(m) + 1.0L));
}

inline DegreeSensitivityData degree_sensitivity_equator(real nu, int m,
                                                         real theta) {
  DegreeSensitivityData out;
  if (m < 1 || !(theta > 0.75L * half_pi && theta <= half_pi))
    return out;
  // The legacy equatorial normalization is only the selected value branch
  // once at least one full order separates nu and m.  For 0 <= nu-m < 1 the
  // regular polar branch must be continued instead.
  if (!(nu - static_cast<real>(m) >= 1.0L))
    return out;

  real value = 0.0L, derivative = 0.0L;
  real value_nu = 0.0L, derivative_nu = 0.0L;
  if (!equator_legacy_sensitivities(nu, m, value, derivative, value_nu,
                                    derivative_nu))
    return out;
  std::array<real, 4> state{value, derivative, value_nu, derivative_nu};

  auto rhs = [=](real x, const std::array<real, 4> &y, std::array<real, 4> &f) {
    const real st = std::sin(x);
    const real cot = std::cos(x) / st;
    const real coefficient =
        nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
    f[0] = y[1];
    f[1] = -cot * y[1] - coefficient * y[0];
    f[2] = y[3];
    f[3] = -cot * y[3] - coefficient * y[2] - (2.0L * nu + 1.0L) * y[0];
  };
  if (theta < half_pi && !dopri5_integrate(half_pi, theta, state, rhs, 2.0e-15L,
                                           2.0e-18L, 0.0002L))
    return out;

  out.center.value = state[0];
  out.center.dtheta = state[1];
  out.value_nu = state[2];
  out.dtheta_nu = state[3];
  const real st = std::sin(theta);
  const real cot = std::cos(theta) / st;
  const real coefficient =
      nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
  out.center.d2theta = -cot * state[1] - coefficient * state[0];
  out.center.status = status_ok;
  out.log_amplitude = 0.0L;
  out.scaled_value = state[0];
  out.scaled_dtheta = state[1];
  out.scaled_value_nu = state[2];
  out.scaled_dtheta_nu = state[3];
  out.valid = finite(state[0]) && finite(state[1]) && finite(state[2]) &&
              finite(state[3]);
  out.materialized = out.valid;
  return out;
}

// Integrate exact analytic degree sensitivity; no finite differences in nu.
inline DegreeSensitivityData degree_sensitivity_ode(real nu, int m,
                                                    real theta) {
  DegreeSensitivityData out;
  if (m < 0 || !(nu >= static_cast<real>(m)) || !(theta > 0.0L && theta < pi))
    return out;

  auto order_rescue = [&](int failure_stage) {
    DegreeSensitivityData rescued =
        degree_sensitivity_order_downward(nu, m, theta);
    if (rescued.valid)
      return rescued;
    DegreeSensitivityData failed;
    failed.failure_stage = failure_stage;
    return failed;
  };

  const real kappa = nu + 0.5L;
  if (kappa > 512.0L && theta <= 0.02L) {
    out = degree_sensitivity_order_downward(nu, m, theta);
    if (out.valid)
      return out;
  }
  if ((kappa <= 512.0L ||
       (m >= 2 && static_cast<real>(m) >= 0.75L * kappa * std::sin(theta)))) {
    out = degree_sensitivity_turning_ivp(nu, m, theta);
    if (out.valid)
      return out;
    out = degree_sensitivity_order_downward(nu, m, theta);
    if (out.valid)
      return out;
  }

  if (m >= 512) {
    const real eta2 = static_cast<real>(m) * m - 0.25L;
    if (eta2 > 0.0L && eta2 < kappa * kappa) {
      const real turning = std::asin(std::sqrt(eta2) / kappa);
      if (theta <= turning + 0.12L) {
        out = degree_sensitivity_order_downward(nu, m, theta);
        if (out.valid)
          return out;
      }
    }
  }

  out = degree_sensitivity_equator(nu, m, theta);
  if (out.valid)
    return out;

  real start = choose_series_start(kappa, m, theta);
  PolarData polar;
  bool got = false;
  for (int attempt = 0; attempt < 12; ++attempt) {
    polar = polar_series(kappa, m, start);
    if (polar.converged) {
      got = true;
      break;
    }
    start *= 0.5L;
    if (start < 1.0e-12L)
      break;
  }
  if (!got)
    return order_rescue(1);

  real requested_match = choose_match_angle(kappa, m, start);
  if (m >= 1) {
    const real turning_ratio =
        std::sqrt(static_cast<real>(m) * m - 0.25L) / kappa;
    if (turning_ratio < 1.0L) {
      const real turning = std::asin(turning_ratio);
      requested_match = std::max(start, turning - 0.12L);
    }
  }
  const real match = std::min(theta, requested_match);
  if (match > start && !propagate_logarithmic(kappa, m, start, match, polar)) {
    return order_rescue(2);
  }

  const real scale_log = legacy_log_scale(nu, m);
  const real scale_log_nu = legacy_log_scale_nu(nu, m);
  real log_amplitude = polar.log_abs_L + scale_log;
  const real log_amplitude_nu = polar.log_abs_L_kappa + scale_log_nu;
  const real sm = std::sin(match);
  if (!(sm > 0.0L))
    return order_rescue(3);

  const real sign = static_cast<real>(polar.sign);
  const real logarithmic_theta = polar.r_t / sm;
  std::array<real, 4> state{
      sign, sign * logarithmic_theta, sign * log_amplitude_nu,
      sign * (polar.r_t_kappa / sm + logarithmic_theta * log_amplitude_nu)};

  if (theta > match) {
    auto rhs = [=](real x, const std::array<real, 4> &y,
                   std::array<real, 4> &f) {
      const real st = std::sin(x);
      const real cot = std::cos(x) / st;
      const real coefficient =
          nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
      f[0] = y[1];
      f[1] = -cot * y[1] - coefficient * y[0];
      f[2] = y[3];
      f[3] = -cot * y[3] - coefficient * y[2] - (2.0L * nu + 1.0L) * y[0];
    };
    real current = match;
    while (current < theta) {
      const real next = std::min(theta, current + 0.002L);
      if (!dopri5_integrate(current, next, state, rhs, 2.0e-15L, 2.0e-18L,
                            0.0005L)) {
        return order_rescue(4);
      }
      const real magnitude =
          std::max({std::fabs(state[0]), std::fabs(state[1]),
                    std::fabs(state[2]), std::fabs(state[3])});
      if (!(magnitude > 0.0L) || !finite(magnitude)) {
        return order_rescue(4);
      }
      for (real &component : state)
        component /= magnitude;
      log_amplitude += std::log(magnitude);
      current = next;
    }
  }

  out.log_amplitude = log_amplitude;
  out.scaled_value = state[0];
  out.scaled_dtheta = state[1];
  out.scaled_value_nu = state[2];
  out.scaled_dtheta_nu = state[3];
  out.valid = finite(log_amplitude) && finite(state[0]) && finite(state[1]) &&
              finite(state[2]) && finite(state[3]);
  if (!out.valid)
    return order_rescue(5);

  const real amplitude = safe_exp(log_amplitude);
  if (!(amplitude > 0.0L) || !finite(amplitude)) {
    out.center.status = status_range;
    return out;
  }
  out.center.value = amplitude * state[0];
  out.center.dtheta = amplitude * state[1];
  out.value_nu = amplitude * state[2];
  out.dtheta_nu = amplitude * state[3];
  const real st = std::sin(theta);
  const real cot = std::cos(theta) / st;
  const real coefficient =
      nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
  out.center.d2theta =
      -cot * out.center.dtheta - coefficient * out.center.value;
  out.center.status = status_ok;
  out.materialized = finite(out.center.value) && finite(out.center.dtheta) &&
                     finite(out.value_nu) && finite(out.dtheta_nu);
  return out.materialized ? out : order_rescue(6);
}

struct BoundaryPhase {
  real Phi = 0.0L;
  real Phi_kappa = 0.0L;
  real log_C = 0.0L;
  real phase_q = 0.0L;
  real phase_q_kappa = 0.0L;
  real phase_a_kappa = 0.0L;
  real boundary_residual = 0.0L;
  int status = status_ok;
};

inline BoundaryPhase boundary_phase(real nu, int m, real theta0) {
  BoundaryPhase out;
  CoreState core = prepare_core(nu, m, theta0);
  if (core.status != status_ok) {
    out.status = core.status;
    return out;
  }
  PhaseState phase;
  if (!phase_state_at(nu + 0.5L, m, theta0, phase)) {
    out.status = status_integrator;
    return out;
  }
  out.Phi = phase.alpha + core.match.delta;
  out.Phi_kappa = phase.alpha_kappa + core.match.delta_kappa;
  out.log_C = core.match.log_C;
  out.phase_q = phase.q;
  out.phase_q_kappa = phase.q_kappa;
  out.phase_a_kappa = phase.a_kappa;
  out.boundary_residual = std::fabs(std::sin(out.Phi));
  if (!(out.Phi_kappa > 0.0L) || !finite(out.Phi_kappa))
    out.status = status_phase;
  return out;
}

struct BoundaryDegreeData {
  EvalResult center;
  real value_nu = 0.0L;
  real dtheta_nu = 0.0L;
  bool valid = false;
};

inline BoundaryDegreeData endpoint_degree_data(real nu, int m, real theta) {
  BoundaryDegreeData out;
  if (m < 0 || nu < static_cast<real>(m) ||
      !(theta > 0.0L && theta <= half_pi))
    return out;
  const real z = sqr(std::sin(0.5L * theta));
  const real parameter = sqr((nu + 0.5L) * std::sin(0.5L * theta)) /
                         (1.0L + static_cast<real>(m));
  if (parameter > 4.0L)
    return out;

  real term = 1.0L, term_n = 0.0L;
  real sum = 1.0L, sum_n = 0.0L, sum_z = 0.0L, sum_nz = 0.0L;
  real c_sum = 0.0L, c_n = 0.0L, c_z = 0.0L, c_nz = 0.0L;
  const auto add = [](real value, real &total, real &compensation) {
    const real adjusted = value - compensation;
    const real updated = total + adjusted;
    compensation = (updated - total) - adjusted;
    total = updated;
  };
  bool converged = false;
  for (int k = 1; k <= 300; ++k) {
    const real kr = static_cast<real>(k);
    const real denominator = kr * (static_cast<real>(m) + kr);
    const real ratio =
        ((kr - 1.0L - nu) * (kr + nu) / denominator) * z;
    const real ratio_n = -(2.0L * nu + 1.0L) / denominator * z;
    const real previous = term;
    term = previous * ratio;
    term_n = term_n * ratio + previous * ratio_n;
    add(term, sum, c_sum);
    add(term_n, sum_n, c_n);
    add(kr * term / z, sum_z, c_z);
    add(kr * term_n / z, sum_nz, c_nz);
    const real tolerance = 16.0L * std::numeric_limits<real>::epsilon();
    if (std::fabs(term) <= tolerance * std::max<real>(1.0L, std::fabs(sum)) &&
        std::fabs(term_n) <=
            tolerance * std::max<real>(1.0L, std::fabs(sum_n))) {
      converged = true;
      break;
    }
  }
  if (!converged || !finite(sum) || !finite(sum_n) || !finite(sum_z) ||
      !finite(sum_nz))
    return out;

  const real mr = static_cast<real>(m);
  const real longitude_factor = m == 0 ? 1.0L : 2.0L;
  const real log_prefactor =
      0.5L * (std::log(longitude_factor) + std::lgamma(nu + mr + 1.0L) -
              std::lgamma(nu - mr + 1.0L)) +
      mr * std::log(std::tan(0.5L * theta)) - std::lgamma(mr + 1.0L);
  const real prefactor = safe_exp(log_prefactor);
  const real prefactor_n = 0.5L *
      (digamma_positive(nu + mr + 1.0L) -
       digamma_positive(nu - mr + 1.0L));
  const real z_theta = 0.5L * std::sin(theta);
  const real polar_rate = mr / std::sin(theta);
  const real angular_sum = polar_rate * sum + z_theta * sum_z;
  const real angular_sum_n = polar_rate * sum_n + z_theta * sum_nz;

  out.center.value = prefactor * sum;
  out.center.dtheta = prefactor * angular_sum;
  out.value_nu = prefactor * (prefactor_n * sum + sum_n);
  out.dtheta_nu =
      prefactor * (prefactor_n * angular_sum + angular_sum_n);
  const real st = std::sin(theta), ct = std::cos(theta);
  const real coefficient = nu * (nu + 1.0L) - mr * mr / (st * st);
  out.center.d2theta =
      -ct / st * out.center.dtheta - coefficient * out.center.value;
  out.center.status = status_ok;
  out.valid = finite(out.center.value) && finite(out.center.dtheta) &&
              finite(out.value_nu) && finite(out.dtheta_nu);
  return out;
}

inline BoundaryDegreeData boundary_degree_data(real nu, int m, real theta0) {
  const BoundaryDegreeData endpoint = endpoint_degree_data(nu, m, theta0);
  if (endpoint.valid)
    return endpoint;
  BoundaryDegreeData out;
  const DegreeSensitivityData sensitivity =
      degree_sensitivity_ode(nu, m, theta0);
  if (!sensitivity.valid || !sensitivity.materialized)
    return out;
  out.center = sensitivity.center;
  out.value_nu = sensitivity.value_nu;
  out.dtheta_nu = sensitivity.dtheta_nu;
  out.valid = true;
  return out;
}

// Independent cap integral with Lagrange initialization and segment rescaling.
inline real cap_norm_turning_ivp_impl(int m, real nu, real theta0) {
  if (m < 1 || !(nu >= static_cast<real>(m)) ||
      !(theta0 > 0.0L && theta0 <= half_pi))
    return std::numeric_limits<real>::quiet_NaN();
  const real kappa = nu + 0.5L;
  const real eta2 = static_cast<real>(m) * m - 0.25L;
  if (!(eta2 > 0.0L) || !(eta2 < kappa * kappa))
    return std::numeric_limits<real>::quiet_NaN();
  const real turning = std::asin(std::sqrt(eta2) / kappa);
  const real match =
      std::max<real>(1.0e-8L, turning - std::max<real>(0.02L, 12.0L / kappa));
  if (!(theta0 > match))
    return std::numeric_limits<real>::quiet_NaN();

  const DegreeSensitivityData seed =
      degree_sensitivity_order_downward(nu, m, match);
  if (!seed.valid || !seed.materialized)
    return std::numeric_limits<real>::quiet_NaN();
  std::array<real, 5> state{seed.center.value, seed.center.dtheta,
                            seed.value_nu, seed.dtheta_nu, 0.0L};
  real amplitude = std::max({std::fabs(state[0]), std::fabs(state[1]),
                             std::fabs(state[2]), std::fabs(state[3])});
  if (!(amplitude > 0.0L) || !finite(amplitude))
    return std::numeric_limits<real>::quiet_NaN();
  for (int i = 0; i < 4; ++i)
    state[i] /= amplitude;
  real norm = std::sin(match) * (state[2] * state[1] - state[0] * state[3]) /
              (2.0L * nu + 1.0L) * amplitude * amplitude;

  auto rhs = [=](real x, const std::array<real, 5> &y, std::array<real, 5> &f) {
    const real st = std::sin(x);
    const real cot = std::cos(x) / st;
    const real coefficient =
        nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
    f[0] = y[1];
    f[1] = -cot * y[1] - coefficient * y[0];
    f[2] = y[3];
    f[3] = -cot * y[3] - coefficient * y[2] - (2.0L * nu + 1.0L) * y[0];
    f[4] = y[0] * y[0] * st;
  };

  real current = match;
  while (current < theta0) {
    const real next = std::min(theta0, current + 0.001L);
    state[4] = 0.0L;
    if (!dopri5_integrate(current, next, state, rhs, 5.0e-16L, 1.0e-20L,
                          0.0001L))
      return std::numeric_limits<real>::quiet_NaN();
    norm += amplitude * amplitude * state[4];
    const real local_scale =
        std::max({std::fabs(state[0]), std::fabs(state[1]), std::fabs(state[2]),
                  std::fabs(state[3])});
    if (!(local_scale > 0.0L) || !finite(local_scale))
      return std::numeric_limits<real>::quiet_NaN();
    for (int i = 0; i < 4; ++i)
      state[i] /= local_scale;
    amplitude *= local_scale;
    if (!(amplitude > 0.0L) || !finite(amplitude))
      return std::numeric_limits<real>::quiet_NaN();
    current = next;
  }
  return (norm > 0.0L && finite(norm)) ? norm
                                       : std::numeric_limits<real>::quiet_NaN();
}

inline real cap_norm_order_downward_log(int m, real nu, real theta0) {
  const DegreeSensitivityData data =
      degree_sensitivity_order_downward(nu, m, theta0);
  if (!data.valid || !(data.scaled_dtheta != 0.0L))
    return std::numeric_limits<real>::quiet_NaN();
  const real wronskian = data.scaled_dtheta * data.scaled_value_nu -
                         data.scaled_value * data.scaled_dtheta_nu;
  if (!(wronskian > 0.0L) || !finite(wronskian))
    return std::numeric_limits<real>::quiet_NaN();
  const real log_norm = std::log(std::sin(theta0)) -
                        std::log(2.0L * nu + 1.0L) + 2.0L * data.log_amplitude +
                        std::log(wronskian);
  return safe_exp(log_norm);
}

inline real cap_norm_polar(int m, real nu, real theta0) {
  const real z = sqr(std::sin(theta0 / 2.0L));
  if (z == 0.0L)
    return 0.0L;
  std::array<real, 64> h{}, b{};
  h[0] = b[0] = 1.0L;
  int last = 0;
  for (int k = 1; k < 64; ++k) {
    const real kr = k;
    h[k] = h[k-1] * ((m + kr - 1.0L - nu) * (m + kr + nu)) *
           z / (kr * (m + kr));
    b[k] = b[k-1] * (kr - 1.0L - m) * z / kr;
    last = k;
    if (k >= 8 && std::fabs(h[k]) + std::fabs(b[k]) < 1.0e-20L)
      break;
    if (k == 63)
      return std::numeric_limits<real>::quiet_NaN();
  }
  real sum = 0.0L, correction = 0.0L;
  for (int i = 0; i <= last; ++i)
    for (int j = 0; j <= last; ++j)
      for (int k = 0; k <= last; ++k) {
        const real term = h[i] * h[j] * b[k] / (m + i + j + k + 1.0L);
        const real adjusted = term - correction;
        const real updated = sum + adjusted;
        correction = (updated - sum) - adjusted;
        sum = updated;
      }
  if (!(sum > 0.0L) || !finite(sum))
    return std::numeric_limits<real>::quiet_NaN();
  const auto amplitude = pole_log_scale(nu, m);
  return safe_exp(2.0L * (amplitude.hi + amplitude.lo) +
      m * (2.0L * log2_value + std::log(z)) +
      log2_value + std::log(z) + std::log(sum));
}

inline real cap_norm_legacy_impl(int m, real nu, real theta0) {
  if ((nu + 0.5L) * theta0 <= 1.5L)
    return cap_norm_polar(m, nu, theta0);
  BoundaryPhase b = boundary_phase(nu, m, theta0);
  if (b.status == status_ok && std::fabs(std::sin(b.Phi)) <= 2.0e-7L) {
    const real kappa = nu + 0.5L;
    // Full Lagrange identity for U=C*exp(-q/2)*sin(Phi).  The
    // amplitude sensitivity cancels; the q/a sensitivity terms vanish
    // only at an exact Dirichlet zero and must be retained near a root.
    const real sine = std::sin(b.Phi), cosine = std::cos(b.Phi);
    const real norm_factor = b.Phi_kappa -
        b.phase_q_kappa * sine * cosine +
        0.5L * b.phase_a_kappa * safe_exp(-b.phase_q) * sine * sine;
    if (norm_factor > 0.0L && finite(norm_factor)) {
      const real log_norm_amp =
          2.0L * b.log_C + std::log(norm_factor) - std::log(2.0L * kappa);
      const real log_norm_legacy = log_norm_amp + 2.0L * legacy_log_scale(nu, m);
      return safe_exp(log_norm_legacy);
    }
  }

  const BoundaryDegreeData data = boundary_degree_data(nu, m, theta0);
  if (!data.valid && m >= 512) {
    const real unified = cap_norm_turning_ivp_impl(m, nu, theta0);
    if (finite(unified) && unified > 0.0L)
      return unified;
  }
  if (!data.valid) {
    const real logarithmic = cap_norm_order_downward_log(m, nu, theta0);
    return (logarithmic > 0.0L && finite(logarithmic))
               ? logarithmic
               : std::numeric_limits<real>::quiet_NaN();
  }
  const real norm = std::sin(theta0) *
                    (data.center.dtheta * data.value_nu -
                     data.center.value * data.dtheta_nu) /
                    (2.0L * nu + 1.0L);
  if (norm > 0.0L && finite(norm))
    return norm;
  const real logarithmic = cap_norm_order_downward_log(m, nu, theta0);
  return (logarithmic > 0.0L && finite(logarithmic))
             ? logarithmic
             : std::numeric_limits<real>::quiet_NaN();
}

inline real boundary_degree_derivative_legacy_impl(int m, real nu,
                                                   real theta0) {
  const BoundaryDegreeData data = boundary_degree_data(nu, m, theta0);
  return data.valid ? data.value_nu : std::numeric_limits<real>::quiet_NaN();
}

inline bool evaluate_table_legacy(real nu, int m, real theta,
                                  EvalResult &out) {
#ifndef RSCHA_CR_TABLE_BUILD
  double p_tilde = 0.0, p_theta = 0.0, p_theta2 = 0.0;
  if (!CRTable::legendre(static_cast<double>(m), static_cast<double>(nu),
                         static_cast<double>(theta), p_tilde, p_theta,
                         p_theta2))
    return false;
  const real st = std::sin(theta);
  const real ct = std::cos(theta);
  const real longitude_factor = m == 0 ? 1.0L : 2.0L;
  const real factor =
      std::sqrt(longitude_factor / ((nu + 0.5L) * st));
  out.value = factor * static_cast<real>(p_tilde);
  out.dtheta = factor *
               (static_cast<real>(p_theta) - 0.5L * ct / st * p_tilde);
  const real coefficient =
      nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
  out.d2theta = -ct / st * out.dtheta - coefficient * out.value;
  out.residual = 0.0L;
  out.status = status_ok;
  return finite(out.value) && finite(out.dtheta) && finite(out.d2theta);
#else
  (void)nu;
  (void)m;
  (void)theta;
  (void)out;
  return false;
#endif
}

inline bool evaluate_global_airy_legacy(real nu, int m, real theta,
                                        EvalResult &out,
                                        bool value_only = false) {
#ifndef RSCHA_CR_TABLE_BUILD
  const double degree = static_cast<double>(nu);
  const double angle = static_cast<double>(theta);
  const double carrier = degree + 0.5;
  if (!(degree >= 1000.0 && degree <= 100000000.0 && m >= 0 && m <= degree &&
        angle > 0.0 && angle <= 0.5 * static_cast<double>(pi)))
    return false;
  struct ExtendedLogValue {
    int sign = 0;
    real log_abs = 0.0L;
    real log_abs_lo = 0.0L;
  };
  auto evanescent_p = [&](int order, ExtendedLogValue &p) {
    const real carrier_extended = nu + 0.5L;
    const real alpha_extended = static_cast<real>(order) / carrier_extended;
    const double alpha = static_cast<double>(alpha_extended);
    if (!((alpha >= 0.1 || order >= 32) && alpha <= 0.9999 &&
          angle < std::asin(alpha)))
      return false;
    const auto map = crli::detail::canonical_map(alpha, angle);
    if (map.region != crli::Region::evanescent ||
        !(map.zeta_over_q > 0.0))
      return false;
    const long double theta_ld = theta;
    const long double sine = std::sin(theta_ld);
    const long double cosine = std::cos(theta_ld);
    const long double order_extended = static_cast<long double>(order);
    using crli::detail::DoubleDouble;
    const DoubleDouble carrier_dd{static_cast<double>(carrier_extended),0.0};
    const DoubleDouble order_dd{static_cast<double>(order_extended),0.0};
    const DoubleDouble sine_dd{static_cast<double>(sine),0.0};
    const DoubleDouble cosine_dd{static_cast<double>(cosine),0.0};
    const DoubleDouble one_dd{1.0,0.0};
    const DoubleDouble k_dd=(sine_dd*sine_dd*(carrier_dd-order_dd)*
        (carrier_dd+order_dd))/(order_dd*order_dd*cosine_dd*cosine_dd);
    const DoubleDouble one_minus_t_dd=k_dd/
        (one_dd+crli::detail::dd_sqrt_positive(one_dd-k_dd));
    const DoubleDouble atanh_t_dd=(
        crli::detail::dd_log_positive(DoubleDouble{2.0,0.0}-one_minus_t_dd)-
        crli::detail::dd_log_positive(one_minus_t_dd))*0.5;
    const DoubleDouble at_denominator_dd=
        (carrier_dd-order_dd)+order_dd*one_minus_t_dd;
    const DoubleDouble at_numerator_dd=
        (carrier_dd+order_dd)-order_dd*one_minus_t_dd;
    const DoubleDouble atanh_at_dd=(
        crli::detail::dd_log_positive(at_numerator_dd)-
        crli::detail::dd_log_positive(at_denominator_dd))*0.5;
    const DoubleDouble action_dd=order_dd*atanh_t_dd-
                                 carrier_dd*atanh_at_dd;
    const long double action_extended=action_dd.hi;
    const long double q_extended =
        (carrier_extended * sine - order_extended) *
        (carrier_extended * sine + order_extended) /
        (carrier_extended * carrier_extended * sine * sine);
    const long double zeta_extended = -std::pow(
        1.5L * action_extended / carrier_extended, 2.0L / 3.0L);
    const long double zeta_over_q_extended = zeta_extended / q_extended;
    if (!(action_extended >= 0.0L && zeta_over_q_extended > 0.0L) ||
        !finite(action_extended) || !finite(zeta_over_q_extended))
      return false;
    const double action = static_cast<double>(action_extended);
    const double argument = std::pow(1.5 * action, 2.0 / 3.0);
    const auto airy = crli::detail::airy(argument, action);
    double near_correction = 0.0;
    if (CRTable::legendre_near_evanescent(
            1.0 / carrier, degree - order, std::log1p(action),
            near_correction)) {
      const double prefactor =
          std::sqrt(2.0) * std::pow(carrier, 1.0 / 6.0) *
          std::pow(static_cast<double>(zeta_over_q_extended), 0.25);
      const double scaled = prefactor * airy.ai_scaled;
      p.sign = (scaled > 0.0) - (scaled < 0.0);
      const real local = std::log(std::fabs(static_cast<real>(scaled))) +
                         static_cast<real>(near_correction);
      const DoubleDouble log_value=
          DoubleDouble{static_cast<double>(local),0.0}-action_dd;
      p.log_abs = log_value.hi;
      p.log_abs_lo = log_value.lo;
      return p.sign != 0 && finite(p.log_abs);
    }
    const double b0 = crli::detail::airy_transport_b0_fast(alpha, angle, map);
    const double relative_a1 = carrier < 4096.0
        ? crli::detail::airy_transport_a1_fast(alpha, angle) /
              (carrier * carrier)
        : 0.0;
    const double scale = 1.0 + (std::abs(relative_a1) <= 0.01
                              ? relative_a1 : 0.0);
    // B0 remains material when mu/carrier is small, even at large
    // carrier. Carrier alone is not a valid asymptotic truncation test.
    const double derivative_scale=std::pow(carrier,-4.0/3.0)*b0;
    const double corrected =
        scale * airy.ai_scaled - derivative_scale * airy.aip_scaled;
    const double prefactor =
        std::sqrt(2.0) * std::pow(carrier, 1.0 / 6.0) *
        std::pow(static_cast<double>(zeta_over_q_extended), 0.25);
    const double scaled = prefactor * corrected;
    p.sign = (scaled > 0.0) - (scaled < 0.0);
    const real log_scaled = std::log(std::fabs(static_cast<real>(scaled)));
    const DoubleDouble log_value=
        DoubleDouble{static_cast<double>(log_scaled),0.0}-action_dd;
    p.log_abs = log_value.hi;
    p.log_abs_lo = log_value.lo;
    return p.sign != 0 && finite(p.log_abs);
  };
  if (angle < std::asin(static_cast<double>(m) / carrier)) {
    ExtendedLogValue primary;
    if (!evanescent_p(m, primary))
      return false;
    const real st = std::sin(theta);
    const real ct = std::cos(theta);
    const real longitude = m == 0 ? 1.0L : 2.0L;
    const real log_factor = 0.5L *
        (std::log(longitude) - std::log(nu + 0.5L) - std::log(st));
    const crli::detail::DoubleDouble normalized_log=
        crli::detail::DoubleDouble{static_cast<double>(primary.log_abs),
                                   static_cast<double>(primary.log_abs_lo)}+
        crli::detail::DoubleDouble{static_cast<double>(log_factor),0.0};
    out.log_abs_value = normalized_log.hi;
    out.log_abs_value_lo = normalized_log.lo;
    out.sign = primary.sign;
    out.value = static_cast<real>(out.sign) * safe_exp(out.log_abs_value);
    if (value_only) {
      out.status = status_ok;
      return finite(out.log_abs_value);
    }
    if (m >= 10000 && nu >= 90000.0L &&
        theta < 0.5L * std::asin(static_cast<real>(m)/(nu+0.5L))) {
      // In a deeply recessive region, subtracting adjacent-degree logs
      // amplifies their rounding by nu/m. Differentiate the smooth log
      // amplitude instead, retaining both log parts in each difference.
      const real h=theta*0.001L;
      crli::detail::DoubleDouble differences[3];
      for(int j=1;j<=3;++j) {
        EvalResult left,right;
        if(!evaluate_global_airy_legacy(nu,m,theta-j*h,left,true) ||
           !evaluate_global_airy_legacy(nu,m,theta+j*h,right,true) ||
           left.sign!=out.sign || right.sign!=out.sign) return false;
        differences[j-1]=crli::detail::DoubleDouble{static_cast<double>(right.log_abs_value),static_cast<double>(right.log_abs_value_lo)}-
            crli::detail::DoubleDouble{static_cast<double>(left.log_abs_value),static_cast<double>(left.log_abs_value_lo)};
      }
      const auto numerator=differences[0]*45.0-differences[1]*9.0+differences[2];
      out.dtheta_over_value=(numerator.hi+numerator.lo)/(60.0L*h);
      const real coefficient=nu*(nu+1.0L)-static_cast<real>(m)*m/(st*st);
      out.d2theta_over_value=-(ct/st)*out.dtheta_over_value-coefficient;
      out.dtheta=out.value*out.dtheta_over_value;
      out.d2theta=out.value*out.d2theta_over_value;
      out.status=status_ok;
      return finite(out.dtheta_over_value)&&finite(out.d2theta_over_value);
    }
    if (nu >= static_cast<real>(m + 2)) {
      EvalResult adjacent, previous;
      if (!evaluate_global_airy_legacy(nu - 1.0L, m, theta,
                                       adjacent, true) ||
          !evaluate_global_airy_legacy(nu - 2.0L, m, theta,
                                       previous, true))
        return false;
      const real coupling = std::sqrt(
          (nu - static_cast<real>(m)) * (nu + static_cast<real>(m)));
      const real adjacent_ratio = static_cast<real>(adjacent.sign * out.sign) *
          safe_exp([&]() { const auto delta = crli::detail::DoubleDouble{static_cast<double>(adjacent.log_abs_value), static_cast<double>(adjacent.log_abs_value_lo)} - crli::detail::DoubleDouble{static_cast<double>(out.log_abs_value), static_cast<double>(out.log_abs_value_lo)}; return delta.hi + delta.lo; }());
      const real adjacent_coupling = std::sqrt(
          (nu - 1.0L - static_cast<real>(m)) *
          (nu - 1.0L + static_cast<real>(m)));
      const real previous_ratio = static_cast<real>(previous.sign * adjacent.sign) *
          safe_exp([&]() { const auto delta = crli::detail::DoubleDouble{static_cast<double>(previous.log_abs_value), static_cast<double>(previous.log_abs_value_lo)} - crli::detail::DoubleDouble{static_cast<double>(adjacent.log_abs_value), static_cast<double>(adjacent.log_abs_value_lo)}; return delta.hi + delta.lo; }());
      const real adjacent_derivative =
          ((nu - 1.0L) * ct - adjacent_coupling * previous_ratio) / st;
      out.dtheta_over_value = (nu * ct - coupling * adjacent_ratio) / st;
      out.d2theta_over_value =
          ((nu - 1.0L) * ct * out.dtheta_over_value - nu * st -
           coupling * adjacent_ratio * adjacent_derivative) / st;
      out.dtheta = out.dtheta_over_value * out.value;
      out.d2theta = out.d2theta_over_value * out.value;
      out.residual = 0.0L;
      out.status = status_ok;
      return finite(out.dtheta_over_value) && finite(out.d2theta_over_value);
    }
    ExtendedLogValue adjacent;
    if (!evanescent_p(m + 1, adjacent))
      return false;
    const real adjacent_ratio =
        static_cast<real>(adjacent.sign * primary.sign) *
        std::exp(static_cast<real>(adjacent.log_abs - primary.log_abs) +
                 0.5L * (std::log(2.0L) - std::log(longitude)));
    const real coupling = std::sqrt(
        (nu - static_cast<real>(m)) * (nu + static_cast<real>(m) + 1.0L) *
        longitude / 2.0L);
    out.dtheta_over_value = static_cast<real>(m) * ct / st -
                            coupling * adjacent_ratio;
    const real coefficient =
        nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
    out.d2theta_over_value =
        -ct / st * out.dtheta_over_value - coefficient;
    out.dtheta = out.dtheta_over_value * out.value;
    out.d2theta = out.d2theta_over_value * out.value;
    out.residual = 0.0L;
    out.status = status_ok;
    return finite(out.log_abs_value) && finite(out.dtheta_over_value) &&
           finite(out.d2theta_over_value);
  }
  auto tabulated_p = [&](double local_degree, int order, double &p) {
    const double local_carrier = local_degree + 0.5;
    const double alpha = order / local_carrier;
    if (!((alpha >= 0.1 || order >= 32) && alpha <= 0.9999))
      return false;
    if (angle < std::asin(alpha)) {
      if (local_degree != degree)
        return false;
      ExtendedLogValue evanescent;
      if (!evanescent_p(order, evanescent) ||
          evanescent.log_abs < std::log(std::numeric_limits<double>::min()))
        return false;
      p = evanescent.sign * std::exp(evanescent.log_abs);
      return std::isfinite(p);
    }
    const auto map = crli::detail::canonical_map(alpha, angle);
    const double action = local_carrier * map.action;
    if (!(action >= 0.0))
      return false;
    double real_correction = 0.0, imag_correction = 0.0;
    const double table_inverse =
        std::max(1.0 / local_carrier, 1.0 / 100000000.5);
    double scaled_zeta = 0.0;
    const double correction_alpha_max = local_carrier < 9999.0 ? 0.9 : 0.95;
    const double correction_action_max = local_carrier < 9999.0 ? 100.0 : 400.0;
    if (action <= 7.5) {
      if (!(alpha >= 0.01 && alpha <= correction_alpha_max))
        return false;
      scaled_zeta = std::pow(1.5 * action, 2.0 / 3.0);
      const bool table_ok = local_carrier > 1000000.5
          ? CRTable::legendre_airy_high(
                table_inverse, alpha, scaled_zeta, real_correction,
                imag_correction)
          : CRTable::legendre_airy(
                table_inverse, alpha, scaled_zeta, real_correction,
                imag_correction);
      if (!table_ok)
        return false;
    } else {
      if (!(alpha >= 0.1 && alpha <= correction_alpha_max) ||
          action > correction_action_max)
        return false;
      bool table_ok = local_carrier > 1000000.5
          ? CRTable::legendre_outer_high(
                table_inverse, alpha, std::log(action), real_correction,
                imag_correction)
          : CRTable::legendre_outer(
                table_inverse, alpha, std::log(action), real_correction,
                imag_correction);
      if (!table_ok && local_carrier >= 10000.5 &&
          local_carrier <= 1000000.5)
        table_ok = CRTable::legendre_outer_patch(
            table_inverse, alpha, std::log(action), real_correction,
            imag_correction);
      if (!table_ok)
        return false;
    }
    const auto base =
        crli::detail::airy_basis_value_fast(local_degree, order, angle);
    if (!base.ok)
      return false;
    const double real =
        1.0 + real_correction / (local_carrier * local_carrier);
    const double imag =
        imag_correction / (local_carrier * local_carrier);
    p = real * base.a - imag * base.b;
    return std::isfinite(p);
  };
  double p_tilde = 0.0;
  if (!tabulated_p(degree, m, p_tilde))
    return false;
  const real st = std::sin(theta);
  const real ct = std::cos(theta);
  const real longitude_factor = m == 0 ? 1.0L : 2.0L;
  const real factor =
      std::sqrt(longitude_factor / ((nu + 0.5L) * st));
  out.value = factor * static_cast<real>(p_tilde);
  if (value_only) {
    out.log_abs_value = out.value == 0.0L
        ? -std::numeric_limits<real>::infinity()
        : std::log(std::fabs(out.value));
    out.sign = (out.value > 0.0L) - (out.value < 0.0L);
    out.status = status_ok;
    return finite(out.log_abs_value);
  }
  if (carrier > 1000000.5 && nu >= static_cast<real>(m + 1)) {
    double adjacent_p = 0.0;
    if (tabulated_p(degree, m + 1, adjacent_p)) {
      const real adjacent_value =
          std::sqrt(2.0L / ((nu + 0.5L) * st)) *
          static_cast<real>(adjacent_p);
      const real coupling = std::sqrt(
          (nu - static_cast<real>(m)) *
          (nu + static_cast<real>(m) + 1.0L) * longitude_factor / 2.0L);
      out.dtheta = static_cast<real>(m) * ct / st * out.value -
                   coupling * adjacent_value;
      const real coefficient =
          nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
      out.d2theta = -ct / st * out.dtheta - coefficient * out.value;
      out.dtheta_over_value = out.value == 0.0L
          ? std::numeric_limits<real>::quiet_NaN()
          : out.dtheta / out.value;
      out.d2theta_over_value = out.value == 0.0L
          ? std::numeric_limits<real>::quiet_NaN()
          : out.d2theta / out.value;
      out.log_abs_value = out.value == 0.0L
          ? -std::numeric_limits<real>::infinity()
          : std::log(std::fabs(out.value));
      out.sign = (out.value > 0.0L) - (out.value < 0.0L);
      out.residual = 0.0L;
      out.status = status_ok;
      return finite(out.log_abs_value) && finite(out.dtheta) &&
             finite(out.d2theta);
    }
  }
  bool differentiated_degree_identity = false;
  if (nu == static_cast<real>(m)) {
    out.dtheta = static_cast<real>(m) * ct / st * out.value;
  } else {
    EvalResult adjacent_result;
    if (!evaluate_global_airy_legacy(nu - 1.0L, m, theta,
                                     adjacent_result, true))
      return false;
    const real adjacent_value = adjacent_result.value;
    const real coupling = std::sqrt(
        (nu - static_cast<real>(m)) * (nu + static_cast<real>(m)));
    out.dtheta =
        (nu * ct * out.value - coupling * adjacent_value) / st;
    real adjacent_derivative = 0.0L;
    if (nu - 1.0L == static_cast<real>(m)) {
      adjacent_derivative =
          static_cast<real>(m) * ct / st * adjacent_value;
    } else if (nu - 1.0L > static_cast<real>(m)) {
      EvalResult previous_result;
      if (!evaluate_global_airy_legacy(nu - 2.0L, m, theta,
                                       previous_result, true))
        return false;
      const real previous_value = previous_result.value;
      const real adjacent_coupling = std::sqrt(
          (nu - 1.0L - static_cast<real>(m)) *
          (nu - 1.0L + static_cast<real>(m)));
      adjacent_derivative =
          ((nu - 1.0L) * ct * adjacent_value -
           adjacent_coupling * previous_value) / st;
    } else {
      return false;
    }
    out.d2theta = ((nu - 1.0L) * ct * out.dtheta -
                   nu * st * out.value -
                   coupling * adjacent_derivative) / st;
    differentiated_degree_identity = true;
  }
  if (!differentiated_degree_identity) {
    const real coefficient =
        nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
    out.d2theta = -ct / st * out.dtheta - coefficient * out.value;
  }
  out.residual = 0.0L;
  out.status = status_ok;
  return finite(out.value) && finite(out.dtheta) && finite(out.d2theta);
#else
  (void)nu;
  (void)m;
  (void)theta;
  (void)out;
  (void)value_only;
  return false;
#endif
}

inline bool evaluate_pole_residual_fast(real nu, int m, real theta,
                                        EvalResult &out) {
  const double degree = static_cast<double>(nu);
  const double angle = static_cast<double>(theta);
  const double carrier = degree + 0.5;
  if (!(degree >= 512.0 && degree <= 100000000.0 && m >= 0 &&
        angle > 0.0 && angle <= 0.5 * static_cast<double>(pi)))
    return false;
  const double alpha = m / carrier;
  if (!(angle > std::asin(alpha) && carrier * std::sin(angle) >= 8.0))
    return false;
  const double parameter = (static_cast<double>(m) * m + 0.25) /
                           (2.0 * carrier * std::sin(angle));
  if (parameter > (m <= 33 ? 8.0 : 2.6))
    return false;
  double p_theta = 0.0, truncation_error = 0.0;
  const auto pair = crli::detail::pole_residual_pair(
      degree, static_cast<double>(m), angle, &p_theta, nullptr,
      &truncation_error);
  if (!pair.ok || pair.p_log_scale != 0.0 || truncation_error > 5.0e-11)
    return false;
  const real st = std::sin(theta), ct = std::cos(theta);
  const real longitude = m == 0 ? 1.0L : 2.0L;
  const real factor = std::sqrt(longitude / ((nu + 0.5L) * st));
  out.value = factor * static_cast<real>(pair.p);
  out.dtheta = factor *
      (static_cast<real>(p_theta) - 0.5L * ct / st * pair.p);
  const real coefficient = nu * (nu + 1.0L) -
                           static_cast<real>(m) * m / (st * st);
  out.d2theta = -ct / st * out.dtheta - coefficient * out.value;
  out.residual = static_cast<real>(truncation_error);
  out.status = status_ok;
  return finite(out.value) && finite(out.dtheta) && finite(out.d2theta);
}

inline bool evaluate_endpoint_p_fast(real nu, int m, real theta,
                                     EvalResult &out) {
  const double degree = static_cast<double>(nu);
  const double angle = static_cast<double>(theta);
  const double carrier = degree + 0.5;
  if (!(degree >= 512.0 && degree <= 100000000.0 && m >= 0 &&
        degree >= m && angle > 0.0 && angle <= 0.5 * static_cast<double>(pi)))
    return false;
  const double endpoint_parameter =
      std::pow(carrier * std::sin(0.5 * angle), 2) / (1.0 + m);
  const double endpoint_limit = m / carrier > 0.9 ? 12.0 : 8.0;
  const bool transformed_near_order =
      degree - m <= 64.0 && angle < std::asin(m / carrier);
  if (endpoint_parameter > endpoint_limit && !transformed_near_order)
    return false;
  const auto primary = crli::detail::endpoint_normalized_p(degree, m, angle);
  if (!primary.ok)
    return false;
  const real st = std::sin(theta), ct = std::cos(theta);
  const real longitude = m == 0 ? 1.0L : 2.0L;
  const real log_factor = 0.5L *
      (std::log(longitude) - std::log(nu + 0.5L) - std::log(st));
  out.log_abs_value = static_cast<real>(primary.log_abs) + log_factor;
  out.sign = primary.sign;
  out.value = static_cast<real>(out.sign) * safe_exp(out.log_abs_value);
  if (nu == static_cast<real>(m)) {
    out.dtheta_over_value = static_cast<real>(m) * ct / st;
  } else {
    const auto adjacent =
        crli::detail::endpoint_normalized_p(degree, m + 1, angle);
    if (!adjacent.ok)
      return false;
    const real ratio = static_cast<real>(adjacent.sign * primary.sign) *
        std::exp(static_cast<real>(adjacent.log_abs - primary.log_abs) +
                 0.5L * (std::log(2.0L) - std::log(longitude)));
    const real coupling = std::sqrt(
        (nu - static_cast<real>(m)) * (nu + static_cast<real>(m) + 1.0L) *
        longitude / 2.0L);
    out.dtheta_over_value = static_cast<real>(m) * ct / st - coupling * ratio;
  }
  const real coefficient = nu * (nu + 1.0L) -
                           static_cast<real>(m) * m / (st * st);
  out.d2theta_over_value = -ct / st * out.dtheta_over_value - coefficient;
  out.dtheta = out.dtheta_over_value * out.value;
  out.d2theta = out.d2theta_over_value * out.value;
  out.residual = 0.0L;
  out.status = status_ok;
  return finite(out.log_abs_value) && finite(out.dtheta_over_value) &&
         finite(out.d2theta_over_value);
}

inline bool evaluate_near_table_legacy(real nu, int m, real theta,
                                       EvalResult &out) {
  const double degree = static_cast<double>(nu);
  const double angle = static_cast<double>(theta);
  const double carrier = degree + 0.5;
  if (!(degree >= 1000.0 && m >= 0 && m + 1 <= degree))
    return false;
  if (degree < 10000.0 && CRTable::legendre_near_enabled()) {
    const double alpha = m / carrier;
    if (degree - m <= 1024.0 && angle < std::asin(alpha)) {
      out = evaluate_legacy(nu, m, theta);
      return out.status == status_ok;
    }
    return false;
  }
  EvalResult raw_primary;
  bool raw_primary_ready = false;
  auto recurrence_quality = [&](int order) {
    EvalResult minus, center, plus;
    if (!evaluate_global_airy_legacy(nu - 1.0L, order, theta, minus) ||
        !evaluate_global_airy_legacy(nu, order, theta, center) ||
        !evaluate_global_airy_legacy(nu + 1.0L, order, theta, plus))
      return 2;
    if (order == m) {
      raw_primary = center;
      raw_primary_ready = true;
    }
    const real a = std::sqrt((nu - order + 1.0L) *
                             (nu + order + 1.0L));
    const real b = (2.0L * nu + 1.0L) * std::cos(theta);
    const real c = std::sqrt((nu + order) * (nu - order));
    const real shift = std::max({std::log(a) + plus.log_abs_value,
                                 std::log(b) + center.log_abs_value,
                                 std::log(c) + minus.log_abs_value});
    const real ta = plus.sign * std::exp(std::log(a) + plus.log_abs_value - shift);
    const real tb = -center.sign * std::exp(std::log(b) + center.log_abs_value - shift);
    const real tc = minus.sign * std::exp(std::log(c) + minus.log_abs_value - shift);
    const real scale = std::fabs(ta) + std::fabs(tb) + std::fabs(tc);
    return scale > 0.0L && std::fabs(ta + tb + tc) / scale > 5.0e-11L
               ? 2
               : 0;
  };
  auto quality = [&](int order) {
    const double alpha = order / carrier;
    if (!(angle < std::asin(alpha)))
      return 0;
    const auto map = crli::detail::canonical_map(alpha, angle);
    // Sampled interpolation discrepancy; source and derivative errors are separate.
    double correction = 0.0, sampled_interpolation_estimate = 0.0;
    if (map.region != crli::Region::evanescent)
      return 0;
    if (order / carrier > 0.2 &&
        (degree - order <= 100000.0 ||
         (degree <= 2000000.0 && degree - order <= 300000.0 &&
          angle <= 0.03 * static_cast<double>(pi) / 180.0)))
      return 2;
    const bool table_hit = CRTable::legendre_near_evanescent_raw(
            1.0 / carrier, degree - order,
            std::log1p(carrier * std::abs(map.action)), correction,
            sampled_interpolation_estimate);
    if (!table_hit) {
      if (!(CRTable::legendre_near_enabled() && degree <= 100000000.0 &&
            degree - order <= 450000.0 && order / carrier > 0.2))
        return 0;
      return recurrence_quality(order);
    }
    if (!(std::isfinite(sampled_interpolation_estimate) && sampled_interpolation_estimate <= 4.0e-11))
      return 2;
    if (degree - order <= 450000.0 && order / carrier > 0.2 &&
        recurrence_quality(order) == 2)
      return 2;
    return 1;
  };
  const int primary = quality(m), adjacent = quality(m + 1);
  if (primary == 0 && adjacent == 0 && raw_primary_ready) {
    out = raw_primary;
    return true;
  }
  if (primary == 1 && adjacent == 1 && raw_primary_ready) {
    out = raw_primary;
    return true;
  }
  if (primary == 1 && adjacent == 1)
    return evaluate_global_airy_legacy(nu, m, theta, out);
  if (primary == 2 || adjacent == 2) {
    const auto p0 = CRNear::recessive_log(degree, m, angle);
    const auto p1 = CRNear::recessive_log(degree, m + 1, angle);
    if (!p0.ok || !p1.ok)
      return false;
    const real st = std::sin(theta), ct = std::cos(theta);
    const real longitude = m == 0 ? 1.0L : 2.0L;
    const real log_factor = 0.5L *
        (std::log(longitude) - std::log(nu + 0.5L) - std::log(st));
    const auto normalized_log=crli::detail::DoubleDouble{
        static_cast<double>(p0.log_abs),static_cast<double>(p0.log_abs_lo)}+
        crli::detail::DoubleDouble{static_cast<double>(log_factor),0.0};
    out.log_abs_value=normalized_log.hi;out.log_abs_value_lo=normalized_log.lo;
    out.sign = p0.sign;
    out.value = static_cast<real>(out.sign) * safe_exp(out.log_abs_value);
    const real ratio = static_cast<real>(p1.sign * p0.sign) *
        std::exp((p1.log_abs-p0.log_abs)+(p1.log_abs_lo-p0.log_abs_lo)+
                 0.5L*(std::log(2.0L)-std::log(longitude)));
    const real coupling = std::sqrt(
        (nu - static_cast<real>(m)) * (nu + static_cast<real>(m) + 1.0L) *
        longitude / 2.0L);
    out.dtheta_over_value = static_cast<real>(m) * ct / st - coupling * ratio;
    const real coefficient = nu * (nu + 1.0L) -
                             static_cast<real>(m) * m / (st * st);
    out.d2theta_over_value = -ct / st * out.dtheta_over_value - coefficient;
    out.dtheta = out.dtheta_over_value * out.value;
    out.d2theta = out.d2theta_over_value * out.value;
    out.residual = 0.0L;
    out.status = status_ok;
    return finite(out.log_abs_value) && finite(out.dtheta_over_value) &&
           finite(out.d2theta_over_value);
  }
  return false;
}

// Laplace moment expansion for the recessive, almost-diagonal solution.
// The expansion is used only well inside its small-moment domain. All three
// derivatives are obtained from the same analytic terms; the ODE is not used.
inline bool evaluate_near_order_moments(real nu,int m,real theta,EvalResult& out) {
  const real g=nu-static_cast<real>(m);
  if(nu<90000.0L || m<1000 || g<0.0L || g>64.0L ||
      !(theta>0.0L && theta<half_pi))return false;
  const real st=std::sin(theta),ct=std::cos(theta),tan=st/ct,cot=ct/st;
  if((g+2)*(g+2)*tan*tan/(m+1.0L)>0.05L)return false;
  real term=1,sum=1,first=0,second=0,previous=1;int quiet=0;
  for(int j=1;j<=256;++j) {
    term*=-(g-2*j+2)*(g-2*j+1)*tan*tan/(4*j*(m+static_cast<real>(j)));
    if(std::fabs(term)>previous)return false;
    previous=std::fabs(term);
    const real d=2*j/(st*ct);
    sum+=term;first+=term*d;
    second+=term*(d*d-2*j*(ct*ct-st*st)/(st*st*ct*ct));
    if(std::fabs(term)*j*j<2.0e-19L*std::fabs(sum))++quiet;else quiet=0;
    if(term==0 || quiet>=3)break;
    if(j==256)return false;
  }
  if(!(sum>0))return false;
  using namespace crli::detail;
  const DoubleDouble c{static_cast<double>(ct),0.0};
  const DoubleDouble log_sine=theta>0.5L
      ? dd_log_positive(DoubleDouble{1.0,0.0}-c*c)*0.5
      : dd_log_positive(DoubleDouble{static_cast<double>(st),0.0});
  const DoubleDouble logarithm=log_sine*static_cast<double>(m)+
      DoubleDouble{static_cast<double>(stable_near_order_log_scale(nu,m)+
          g*std::log(ct)+std::log(sum)),0.0};
  const real b=m*cot-g*tan,h1=first/sum;
  out.log_abs_value=logarithm.hi;out.log_abs_value_lo=logarithm.lo;out.sign=1;
  out.dtheta_over_value=b+h1;
  out.d2theta_over_value=b*b-m/(st*st)-g/(ct*ct)+2*b*h1+second/sum;
  out.value=safe_exp(out.log_abs_value+out.log_abs_value_lo);
  out.dtheta=out.value*out.dtheta_over_value;
  out.d2theta=out.value*out.d2theta_over_value;
  out.status=status_ok;return true;
}


// The low-degree fractional anchors have a nonoscillatory Laplace integral.
// Fixed Gaussian rules avoid the endpoint-series failure between the pole
// moment expansion and the equatorial initial-value chart. Both derivatives
// are differentiated under the integral, independently of the Legendre ODE.
struct NearAnchorGaussRule {
  int count=0;
  std::array<double,64> node{},weight{};
};
inline NearAnchorGaussRule near_anchor_gauss_rule(int count) {
  NearAnchorGaussRule rule;rule.count=count;
  for(int i=0;i<(count+1)/2;++i) {
    double x=std::cos(static_cast<double>(pi)*(i+0.75)/(count+0.5));
    auto polynomial=[&](double t,double& derivative) {
      double previous=1,current=t;
      for(int j=2;j<=count;++j) {
        const double next=((2*j-1)*t*current-(j-1)*previous)/j;
        previous=current;current=next;
      }
      derivative=count*(t*current-previous)/(t*t-1);return current;
    };
    for(int iteration=0;iteration<30;++iteration) {
      double derivative=0;const double value=polynomial(x,derivative);
      const double change=value/derivative;
      x-=change;if(std::abs(change)<2.0e-16)break;
    }
    double derivative=0;polynomial(x,derivative);
    const double weight=2/((1-x*x)*derivative*derivative);
    rule.node[i]=-x;rule.node[count-1-i]=x;
    rule.weight[i]=rule.weight[count-1-i]=weight;
  }
  return rule;
}
inline bool evaluate_fractional_anchor_quadrature(
    real nu,int m,real theta,EvalResult& out) {
  const real g=nu-static_cast<real>(m);
  if(nu<90000 || m<1000 || g<0 || g>2 || !(theta>0&&theta<half_pi))
    return false;
  const double s=std::sin(static_cast<double>(theta)),c=std::cos(static_cast<double>(theta));
  const double lambda=m*c*c/(s*s);
  if(!(lambda>=1&&lambda<=400))return false;
  static const std::array<NearAnchorGaussRule,3> rules{
      near_anchor_gauss_rule(16),near_anchor_gauss_rule(32),near_anchor_gauss_rule(64)};
  const double root=std::sqrt(static_cast<double>(m)),gap=static_cast<double>(g);
  const double normalization=2*std::exp(static_cast<double>(
      shifted_log_gamma_ratio(static_cast<real>(m),1,0.5L)))/
      std::sqrt(static_cast<double>(pi)*m);
  const double edges[]{0,1,2,4,8,12};
  const double cot=c/s,kappa=static_cast<double>(nu)+0.5;
  std::array<double,3> previous{},current{};double discrepancy=0;bool converged=false;
  for(const auto& rule:rules) {
    std::array<CRNear::DoubleDouble,3> sums;
    for(int interval=0;interval<5;++interval) {
      const double half=(edges[interval+1]-edges[interval])/2,
                   center=(edges[interval+1]+edges[interval])/2;
      for(int j=0;j<rule.count;++j) {
        const double u=center+half*rule.node[j];
        const std::complex<double> z(1,s*u/(c*root)),dz(-s/c,u/root);
        const auto power=std::pow(z,gap),lower=power/z;
        const double weight=normalization*half*rule.weight[j]*
            std::exp((m-0.5)*std::log1p(-u*u/m));
        const double h0=power.real(),h1=gap*(lower*dz).real(),
                     h2=gap*(gap-1)*(lower/z*dz*dz).real()-gap*h0;
        sums[0].add(CRNear::DoubleDouble(weight*h0));
        sums[1].add(CRNear::DoubleDouble(weight*h1));
        sums[2].add(CRNear::DoubleDouble(weight*h2));
      }
    }
    const double h0=sums[0].value(),h1=sums[1].value(),h2=sums[2].value();
    current={h0,m*cot*h0+h1,
        (static_cast<double>(m)*m*cot*cot-m/(s*s))*h0+2*m*cot*h1+h2};
    if(!(current[0]>0)||!std::isfinite(current[2]))return false;
    const double scale=std::abs(current[0])+std::abs(current[1])/kappa+
                       std::abs(current[2])/(kappa*kappa);
    if(rule.count>16) {
      discrepancy=std::max({std::abs(current[0]-previous[0]),
          std::abs(current[1]-previous[1])/kappa,
          std::abs(current[2]-previous[2])/(kappa*kappa)})/scale;
      if(discrepancy<=1.0e-12){converged=true;break;}
    }
    previous=current;
  }
  if(!converged)return false;
  using namespace crli::detail;
  const DoubleDouble cosine{c,0};
  const auto logarithm=dd_log_positive(DoubleDouble{1,0}-cosine*cosine)*
      (0.5*static_cast<double>(m))+
      DoubleDouble{static_cast<double>(stable_near_order_log_scale(nu,m))+
          gap*std::log(c)+std::log(current[0]),0};
  out.log_abs_value=logarithm.hi;out.log_abs_value_lo=logarithm.lo;out.sign=1;
  out.dtheta_over_value=current[1]/current[0];
  out.d2theta_over_value=current[2]/current[0];
  out.value=safe_exp(out.log_abs_value+out.log_abs_value_lo);
  out.dtheta=out.value*out.dtheta_over_value;out.d2theta=out.value*out.d2theta_over_value;
  out.residual=discrepancy;out.status=status_ok;return true;
}

inline bool evaluate_near_order_recurrence(real nu, int m, real theta,
                                           EvalResult &out) {
  const real gap=nu-static_cast<real>(m);
  const real carrier=nu+0.5L;
  if(!(nu>=512.0L&&gap>=0.0L&&gap<=8192.0L&&m/carrier>0.9L
       &&theta>0.0L&&theta<=half_pi))
    return false;
  const int steps=static_cast<int>(std::floor(gap));
  const real lower=nu-static_cast<real>(steps);
  const real st=std::sin(theta),ct=std::cos(theta),cot=ct/st;
  auto anchor=[&](real degree,real& value,real& derivative,
                  real& second,crli::detail::DoubleDouble& log_value)->bool{
    EvalResult moment_anchor;
    if(evaluate_near_order_moments(degree,m,theta,moment_anchor) ||
       evaluate_fractional_anchor_quadrature(degree,m,theta,moment_anchor)) {
      log_value={static_cast<double>(moment_anchor.log_abs_value),static_cast<double>(moment_anchor.log_abs_value_lo)};
      value=1;derivative=moment_anchor.dtheta_over_value;
      second=moment_anchor.d2theta_over_value;return true;
    }
    const auto endpoint=crli::detail::endpoint_normalized_p_with_derivative(
        static_cast<double>(degree),static_cast<double>(m),
        static_cast<double>(theta));
    const real longitude=m==0?1.0L:2.0L;
    if(!endpoint.ok||endpoint.value.sign==0){
      auto direct=evaluate_legacy(degree,m,theta);
      finalize_metadata(direct);
      if(direct.status!=status_ok||direct.sign==0||
         !finite(direct.log_abs_value)||
         !finite(direct.dtheta_over_value)||
         !finite(direct.d2theta_over_value))return false;
      log_value={static_cast<double>(direct.log_abs_value),static_cast<double>(direct.log_abs_value_lo)};
      value=static_cast<real>(direct.sign);
      derivative=value*direct.dtheta_over_value;
      second=value*direct.d2theta_over_value;
      return true;
    }
    log_value=crli::detail::DoubleDouble{endpoint.value.log_abs,0.0}+
        crli::detail::DoubleDouble{static_cast<double>(0.5L*(std::log(longitude)-std::log(degree+0.5L)-std::log(st))),0.0};
    const real normalized_ratio=endpoint.derivative.sign==0?0.0L:
        static_cast<real>(endpoint.derivative.sign*endpoint.value.sign)
        *std::exp(static_cast<real>(endpoint.derivative.log_abs
                                    -endpoint.value.log_abs));
    value=static_cast<real>(endpoint.value.sign);
    derivative=value*(normalized_ratio-0.5L*cot);
    const real coefficient=degree*(degree+1.0L)
        -static_cast<real>(m)*m/(st*st);
    second=-cot*derivative-coefficient*value;
    return finite(log_value.hi)&&finite(derivative)&&finite(second);
  };
  real previous=0.0L,previous_d=0.0L,previous_d2=0.0L;
  real value=0.0L,derivative=0.0L,second=0.0L;
  crli::detail::DoubleDouble log0{0.0,0.0},log1{0.0,0.0};
  if(!anchor(lower,previous,previous_d,previous_d2,log0))return false;
  if(steps==0){
    value=previous;derivative=previous_d;second=previous_d2;log1=log0;
  }else{
    if(!anchor(lower+1.0L,value,derivative,second,log1))return false;
    const auto common=log0.hi>log1.hi?log0:log1;
    const auto delta0=log0-common,delta1=log1-common;
    const real scale0=std::exp(delta0.hi+delta0.lo),scale1=std::exp(delta1.hi+delta1.lo);
    previous*=scale0;previous_d*=scale0;previous_d2*=scale0;
    value*=scale1;derivative*=scale1;second*=scale1;
    auto recurrence_log=common;
    for(int index=1;index<steps;++index){
      const real degree=lower+static_cast<real>(index);
      const real denominator=(degree-m+1.0L)*(degree+m+1.0L);
      const real a=(2.0L*degree+1.0L)/std::sqrt(denominator);
      const real b=std::sqrt((degree+m)*(degree-m)/denominator);
      const real next=a*ct*value-b*previous;
      const real next_d=a*(-st*value+ct*derivative)-b*previous_d;
      const real next_d2=a*(-ct*value-2.0L*st*derivative+ct*second)
          -b*previous_d2;
      previous=value;previous_d=derivative;previous_d2=second;
      value=next;derivative=next_d;second=next_d2;
      const real magnitude=std::max({std::fabs(previous),std::fabs(value),
          std::fabs(previous_d),std::fabs(derivative),
          std::fabs(previous_d2),std::fabs(second)});
      if(!(magnitude>0.0L)||!finite(magnitude))return false;
      if(magnitude>1.0e100L||magnitude<1.0e-100L){
        previous/=magnitude;previous_d/=magnitude;previous_d2/=magnitude;
        value/=magnitude;derivative/=magnitude;second/=magnitude;
        recurrence_log=recurrence_log+crli::detail::DoubleDouble{static_cast<double>(std::log(magnitude)),0.0};
      }
    }
    log1=recurrence_log;
  }
  if(value==0.0L)return false;
  const auto retained_log=log1+crli::detail::DoubleDouble{static_cast<double>(std::log(std::fabs(value))),0.0};
  out.log_abs_value=retained_log.hi;out.log_abs_value_lo=retained_log.lo;
  out.sign=value>0.0L?1:-1;
  out.dtheta_over_value=derivative/value;
  out.d2theta_over_value=second/value;
  out.value=static_cast<real>(out.sign)*safe_exp(out.log_abs_value+out.log_abs_value_lo);
  out.dtheta=out.dtheta_over_value*out.value;
  out.d2theta=out.d2theta_over_value*out.value;
  out.residual=0.0L;out.status=status_ok;
  return finite(out.log_abs_value)&&finite(out.dtheta_over_value)
      &&finite(out.d2theta_over_value);
}

inline bool evaluate_canonical_legacy(real nu, int m, real theta,
                                      EvalResult &out) {
  if (!(theta > 0.0L && theta < pi) || m < 0 || nu < static_cast<real>(m) ||
      nu > 1.0e8L)
    return false;
  const crli::Result primary = crli::evaluate(
      static_cast<double>(nu), static_cast<double>(m),
      static_cast<double>(theta), true);
  if (primary.status != crli::Status::ok)
    return false;
  const real st = std::sin(theta);
  const real ct = std::cos(theta);
  const real longitude_factor = m == 0 ? 1.0L : 2.0L;
  const real log_factor =
      0.5L * (std::log(longitude_factor) - std::log(nu + 0.5L) - std::log(st));
  out.log_abs_value = static_cast<real>(primary.p_tilde.log_abs) + log_factor;
  out.sign = primary.p_tilde.sign;
  out.value = static_cast<real>(out.sign) * safe_exp(out.log_abs_value);

  if (primary.p_theta_available) {
    const real normal_ratio = primary.p_theta.sign == 0 ? 0.0L :
        static_cast<real>(primary.p_theta.sign * out.sign) * safe_exp(
            static_cast<real>(primary.p_theta.log_abs-primary.p_tilde.log_abs));
    out.dtheta_over_value=normal_ratio-0.5L*ct/st;
    const real coefficient=nu*(nu+1.0L)-static_cast<real>(m)*m/(st*st);
    out.d2theta_over_value=-(ct/st)*out.dtheta_over_value-coefficient;
    out.dtheta=out.value*out.dtheta_over_value;
    out.d2theta=out.value*out.d2theta_over_value;
    out.status=status_ok;
    return finite(out.log_abs_value)&&finite(out.dtheta_over_value)&&finite(out.d2theta_over_value);
  }
  if (nu == static_cast<real>(m)) {
    out.dtheta_over_value = static_cast<real>(m) * ct / st;
  } else {
    if (nu < static_cast<real>(m + 1))
      return false;
    const crli::Result adjacent = crli::evaluate(
        static_cast<double>(nu - 1.0L), static_cast<double>(m),
        static_cast<double>(theta));
    if (adjacent.status != crli::Status::ok)
      return false;
    const real adjacent_log =
        static_cast<real>(adjacent.p_tilde.log_abs) +
        0.5L * (std::log(longitude_factor) - std::log(nu - 0.5L) -
                std::log(st));
    const real adjacent_ratio =
        static_cast<real>(adjacent.p_tilde.sign * out.sign) *
        safe_exp(adjacent_log - out.log_abs_value);
    const real coupling = std::sqrt(
        (nu - static_cast<real>(m)) * (nu + static_cast<real>(m)));
    out.dtheta_over_value =
        (nu * ct - coupling * adjacent_ratio) / st;
    real adjacent_log_derivative = 0.0L;
    if (nu - 1.0L == static_cast<real>(m)) {
      adjacent_log_derivative = static_cast<real>(m) * ct / st;
    } else {
      const crli::Result previous = crli::evaluate(
          static_cast<double>(nu - 2.0L), static_cast<double>(m),
          static_cast<double>(theta));
      if (previous.status != crli::Status::ok)
        return false;
      const real previous_log =
          static_cast<real>(previous.p_tilde.log_abs) +
          0.5L * (std::log(longitude_factor) - std::log(nu - 1.5L) -
                  std::log(st));
      const real previous_ratio =
          static_cast<real>(previous.p_tilde.sign * adjacent.p_tilde.sign) *
          safe_exp(previous_log - adjacent_log);
      const real adjacent_coupling = std::sqrt(
          (nu - 1.0L - static_cast<real>(m)) *
          (nu - 1.0L + static_cast<real>(m)));
      adjacent_log_derivative =
          ((nu - 1.0L) * ct - adjacent_coupling * previous_ratio) / st;
    }
    out.d2theta_over_value =
        ((nu - 1.0L) * ct * out.dtheta_over_value - nu * st -
         coupling * adjacent_ratio * adjacent_log_derivative) / st;
  }
  if (nu == static_cast<real>(m)) {
    const real coefficient =
        nu * (nu + 1.0L) - static_cast<real>(m) * m / (st * st);
    out.d2theta_over_value = -ct / st * out.dtheta_over_value - coefficient;
  }
  out.dtheta = out.dtheta_over_value * out.value;
  out.d2theta = out.d2theta_over_value * out.value;
  out.residual = 0.0L;
  out.status = status_ok;
  return finite(out.log_abs_value) && finite(out.dtheta) && finite(out.d2theta);
}

}

namespace detail {

// At extremely small positive angles the leading polar term is sufficient:
// all omitted terms are O((nu*theta)^2 + m*theta^2), below binary64 accuracy.
// Construct derivatives in signed-log form; their ratios may overflow even
// when the actual second derivative is finite (in particular at m=2).
inline bool evaluate_tiny_pole_schmidt(real nu,int m,real theta,EvalResult& out) {
  if (!(theta>0.0L && theta<1.0e-150L) || nu>1.0e8L) return false;
  using crli::detail::DoubleDouble;
  using crli::detail::dd_log_positive;
  const DoubleDouble degree{static_cast<double>(nu),0.0};
  const DoubleDouble order{static_cast<double>(m),0.0};
  const DoubleDouble one{1.0,0.0};
  const auto logtheta=dd_log_positive(DoubleDouble{static_cast<double>(theta),0.0});
  const auto logtwo=dd_log_positive(DoubleDouble{2.0,0.0});
  const auto logamp=(pole_log_gamma(degree+order+one)-
      pole_log_gamma(degree-order+one))*0.5-
      pole_log_gamma(order+one)-logtwo*static_cast<double>(m)+
      (m==0 ? DoubleDouble{0.0,0.0} : logtwo*0.5);
  auto component=[](DoubleDouble log,int sign,real& value,real& hi,real& lo,int& output_sign) {
    hi=log.hi;lo=log.lo;output_sign=sign;
    value=sign*safe_exp(static_cast<real>(log.hi)+static_cast<real>(log.lo));
  };
  const auto logarithm=logamp+logtheta*static_cast<double>(m);
  component(logarithm,1,out.value,out.log_abs_value,out.log_abs_value_lo,out.sign);
  const real lambda=nu*(nu+1.0L);
  if(m==0) {
    if(lambda>0.0L) {
      const auto logfactor=dd_log_positive(DoubleDouble{static_cast<double>(lambda/2.0L),0.0});
      component(logfactor+logtheta,-1,out.dtheta,out.log_abs_dtheta,out.log_abs_dtheta_lo,out.dtheta_sign);
      component(logfactor,-1,out.d2theta,out.log_abs_d2theta,out.log_abs_d2theta_lo,out.d2theta_sign);
    }
  } else {
    component(logamp+dd_log_positive(order)+logtheta*static_cast<double>(m-1),1,
        out.dtheta,out.log_abs_dtheta,out.log_abs_dtheta_lo,out.dtheta_sign);
    if(m==1) {
      const auto logfactor=dd_log_positive(DoubleDouble{static_cast<double>(0.75L*lambda-0.5L),0.0});
      component(logamp+logfactor+logtheta,-1,out.d2theta,out.log_abs_d2theta,out.log_abs_d2theta_lo,out.d2theta_sign);
    } else {
      component(logamp+dd_log_positive(order)+dd_log_positive(order-one)+
          logtheta*static_cast<double>(m-2),1,out.d2theta,out.log_abs_d2theta,out.log_abs_d2theta_lo,out.d2theta_sign);
    }
  }
  out.status=status_ok;out.residual=0.0L;
  out.evaluation_method="evaluate_tiny_pole_schmidt";
  return true;
}

inline EvalResult evaluate_legacy_schmidt_raw(double m, double n,
                                              double theta_radians) {
  int mi = 0;
  if (!detail::parse_order(m, mi) ||
      !detail::valid_mode_domain(mi, n, theta_radians)) {
    detail::EvalResult invalid;
    invalid.status = detail::status_domain;
    return invalid;
  }
  detail::EvalResult tiny_pole;
  if (detail::evaluate_tiny_pole_schmidt(n,mi,theta_radians,tiny_pole)) {
    detail::finalize_metadata(tiny_pole);return tiny_pole;
  }
  detail::EvalResult near_order;
  if ((detail::evaluate_integer_near_order(static_cast<detail::real>(n), mi,
                                          static_cast<detail::real>(theta_radians),
                                          near_order) && (near_order.evaluation_method = "evaluate_integer_near_order", true))) {
    detail::finalize_metadata(near_order);
    return near_order;
  }
  detail::EvalResult convergent_pole;
  if(detail::evaluate_convergent_pole_schmidt(n,mi,theta_radians,convergent_pole)) {
    convergent_pole.evaluation_method="evaluate_convergent_pole_schmidt";
    detail::finalize_metadata(convergent_pole);return convergent_pole;
  }
  detail::EvalResult near_pole;
  if ((static_cast<detail::real>(n) + 0.5L) *
          static_cast<detail::real>(theta_radians) <= 2.0L &&
      (detail::evaluate_polar_log(
          static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians),
          detail::legacy_log_scale(static_cast<detail::real>(n), mi),
          near_pole) && (near_pole.evaluation_method = "evaluate_polar_log", true))) {
    detail::finalize_metadata(near_pole);
    return near_pole;
  }
  detail::EvalResult moments;
  if(detail::evaluate_near_order_moments(n,mi,theta_radians,moments)) {
    moments.evaluation_method="evaluate_near_order_moments";
    detail::finalize_metadata(moments);return moments;
  }
  detail::EvalResult accurate_equator;
  if(detail::evaluate_coalescent_wkb(static_cast<detail::real>(n),mi,
      static_cast<detail::real>(theta_radians),accurate_equator)) {
    accurate_equator.evaluation_method="evaluate_coalescent_wkb";
    detail::finalize_metadata(accurate_equator);
    return accurate_equator;
  }
  if (n >= 1000000.0 && theta_radians >= 85.0 * static_cast<double>(detail::pi) / 180.0 &&
      (detail::evaluate_equator_ivp(static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians), accurate_equator) && (accurate_equator.evaluation_method = "evaluate_equator_ivp", true))) {
    detail::finalize_metadata(accurate_equator);
    return accurate_equator;
  }
  detail::EvalResult tabulated;
  {
    const detail::real local_carrier = static_cast<detail::real>(n) + 0.5L;
    if ((local_carrier < 90000.0L &&
         theta_radians > 85.0 * static_cast<double>(detail::pi) / 180.0) &&
        (detail::evaluate_legacy_degree_consistent(
            static_cast<detail::real>(n), mi,
            static_cast<detail::real>(theta_radians), tabulated) && (tabulated.evaluation_method = "evaluate_legacy_degree_consistent", true))) {
      detail::finalize_metadata(tabulated);
      return tabulated;
    }
  }
  detail::EvalResult fractional_near_order;
  const detail::real local_degree = static_cast<detail::real>(n);
  const detail::real local_angle = static_cast<detail::real>(theta_radians);
  const detail::real local_gap = local_degree - static_cast<detail::real>(mi);
  const detail::real local_carrier = local_degree + 0.5L;
  const bool fractional_near_order_domain =
      local_degree >= 512.0L && local_gap >= 0.0L && local_gap <= 64.0L &&
      std::fabs(local_degree - std::round(local_degree)) > 1.0e-12L &&
      static_cast<detail::real>(mi) / local_carrier > 0.99L &&
      local_angle >= 85.0L * detail::pi / 180.0L;
  // A short degree recurrence is accurate and bounded in cost when both
  // fractional anchors satisfy the analytic moment expansion's own guard.
  const detail::real anchor_fraction=local_gap-std::floor(local_gap);
  const bool accurate_moment_anchors=fractional_near_order_domain &&
      local_degree-std::floor(local_gap)>=90000.0L && mi>=1000 &&
      local_angle<detail::half_pi &&
      (anchor_fraction+3.0L)*(anchor_fraction+3.0L)*
          detail::sqr(std::tan(local_angle))/(mi+1.0L)<=0.05L;
  // When the moment guard fails close to the equator, the removed sine
  // carrier has a bounded exponent. Try the endpoint-anchored short
  // recurrence there; an unavailable anchor still falls through to SDIRK.
  const bool bounded_endpoint_anchors=fractional_near_order_domain &&
      static_cast<detail::real>(mi)*detail::sqr(std::cos(local_angle))<=400.0L;
  const bool primary_near_order_domain =
      fractional_near_order_domain &&
      (local_degree <= 100000.0L || accurate_moment_anchors ||
       bounded_endpoint_anchors);
  if (primary_near_order_domain &&
      (detail::evaluate_near_order_recurrence(
          local_degree, mi, local_angle, tabulated) && (tabulated.evaluation_method = "evaluate_near_order_recurrence", true))) {
    detail::finalize_metadata(tabulated);
    return tabulated;
  }
  if (fractional_near_order_domain &&
      (detail::evaluate_fractional_near_order_evanescent(
          local_degree, mi, local_angle, fractional_near_order) && (fractional_near_order.evaluation_method = "evaluate_fractional_near_order_evanescent", true))) {
    detail::finalize_metadata(fractional_near_order);
    return fractional_near_order;
  }
  if (!primary_near_order_domain && fractional_near_order_domain &&
      (detail::evaluate_equator_ivp(
          local_degree, mi, local_angle, fractional_near_order) && (fractional_near_order.evaluation_method = "evaluate_equator_ivp", true))) {
    detail::finalize_metadata(fractional_near_order);
    return fractional_near_order;
  }
  if (!primary_near_order_domain &&
      (detail::evaluate_near_order_recurrence(
          static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians), tabulated) && (tabulated.evaluation_method = "evaluate_near_order_recurrence", true))) {
    detail::finalize_metadata(tabulated);
    return tabulated;
  }
  if (primary_near_order_domain &&
      (detail::evaluate_equator_ivp(
          local_degree, mi, local_angle, fractional_near_order) && (fractional_near_order.evaluation_method = "evaluate_equator_ivp", true))) {
    detail::finalize_metadata(fractional_near_order);
    return fractional_near_order;
  }
  {
    const detail::real local_carrier = static_cast<detail::real>(n) + 0.5L;
    const detail::real local_alpha = static_cast<detail::real>(mi) /
                                     local_carrier;
    const double half_degree = 0.5 * static_cast<double>(detail::pi) / 180.0;
    const double one_degree = static_cast<double>(detail::pi) / 180.0;
    if (local_carrier >= 25000.0L && local_carrier < 70000.0L &&
        local_alpha >= 0.2L && theta_radians < std::asin(local_alpha) &&
        theta_radians >= half_degree && theta_radians <= one_degree) {
      const detail::real degree = static_cast<detail::real>(n);
      const detail::real angle = static_cast<detail::real>(theta_radians);
      std::array<detail::EvalResult, 7> values;
      bool coherent = true;
      for (int offset = -3; offset <= 3; ++offset)
        coherent = coherent && (detail::evaluate_canonical_legacy(
            degree + offset, mi, angle, values[offset + 3]) && (values[offset + 3].evaluation_method = "evaluate_canonical_legacy", true));
      for (int index = 1; coherent && index <= 5; ++index) {
        const detail::real center_degree = degree + index - 3;
        const detail::real a = std::sqrt(
            (center_degree - mi + 1.0L) *
            (center_degree + mi + 1.0L));
        const detail::real b =
            (2.0L * center_degree + 1.0L) * std::cos(angle);
        const detail::real c = std::sqrt(
            (center_degree + mi) * (center_degree - mi));
        const detail::real shift = std::max({
            std::log(a) + values[index + 1].log_abs_value,
            std::log(b) + values[index].log_abs_value,
            std::log(c) + values[index - 1].log_abs_value});
        const detail::real ta = values[index + 1].sign *
            std::exp(std::log(a) + values[index + 1].log_abs_value - shift);
        const detail::real tb = -values[index].sign *
            std::exp(std::log(b) + values[index].log_abs_value - shift);
        const detail::real tc = values[index - 1].sign *
            std::exp(std::log(c) + values[index - 1].log_abs_value - shift);
        const detail::real scale =
            std::fabs(ta) + std::fabs(tb) + std::fabs(tc);
        coherent = scale > 0.0L &&
                   std::fabs(ta + tb + tc) / scale <= 5.0e-11L;
      }
      if (coherent) {
          tabulated = values[3];
          detail::finalize_metadata(tabulated);
          return tabulated;
      }
      if ((detail::evaluate_legacy_degree_consistent(
              degree, mi, angle, tabulated) && (tabulated.evaluation_method = "evaluate_legacy_degree_consistent", true))) {
        detail::finalize_metadata(tabulated);
        return tabulated;
      }
      if ((detail::evaluate_near_table_legacy(
              degree, mi, angle, tabulated) && (tabulated.evaluation_method = "evaluate_near_table_legacy", true))) {
        detail::finalize_metadata(tabulated);
        return tabulated;
      }
    }
  }
  {
    const detail::real local_carrier = static_cast<detail::real>(n) + 0.5L;
    const detail::real local_alpha = static_cast<detail::real>(mi) /
                                     local_carrier;
    if (local_carrier >= 10000.0L && local_alpha >= 0.2L &&
        theta_radians < std::asin(local_alpha) &&
        theta_radians <= 0.5 * static_cast<double>(detail::pi) / 180.0 &&
        (detail::evaluate_near_table_legacy(
            static_cast<detail::real>(n), mi,
            static_cast<detail::real>(theta_radians), tabulated) && (tabulated.evaluation_method = "evaluate_near_table_legacy", true))) {
      detail::finalize_metadata(tabulated);
      return tabulated;
    }
  }
  {
    const detail::real local_carrier = static_cast<detail::real>(n) + 0.5L;
    const detail::real local_alpha = static_cast<detail::real>(mi) /
                                     local_carrier;
    if (local_carrier < 10000.0L) {
      tabulated = detail::evaluate_legacy(
          static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians));
      tabulated.evaluation_method = "evaluate_legacy";
      if (tabulated.status == detail::status_ok) {
        detail::finalize_metadata(tabulated);
        return tabulated;
      }
    }
    if (local_carrier < 10000.0L &&
        (detail::evaluate_endpoint_p_fast(
            static_cast<detail::real>(n), mi,
            static_cast<detail::real>(theta_radians), tabulated) && (tabulated.evaluation_method = "evaluate_endpoint_p_fast", true))) {
      detail::finalize_metadata(tabulated);
      return tabulated;
    }
    const bool low_degree = local_carrier < 10000.0L;
    const bool low_degree_evanescent = local_carrier < 25000.0L &&
        local_alpha >= 0.2L && theta_radians < std::asin(local_alpha) &&
        theta_radians > 0.5 * static_cast<double>(detail::pi) / 180.0;
    const bool resolved_evanescent = local_carrier <= 100001.0L &&
        local_alpha >= 0.1L && theta_radians < std::asin(local_alpha) &&
        theta_radians > static_cast<double>(detail::pi) / 180.0;
    const detail::real local_turning = local_alpha < 1.0L
        ? std::asin(local_alpha) : detail::half_pi;
    const bool resolved_turning_band = local_carrier >= 10000.0L &&
        local_carrier <= 100001.0L && local_alpha >= 0.01L &&
        theta_radians < local_turning &&
        theta_radians >= (local_alpha < 0.05L ? 0.5L : 0.95L) * local_turning;
    if (resolved_turning_band) {
      tabulated = detail::evaluate_legacy(
          static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians));
      tabulated.evaluation_method = "evaluate_legacy";
      if (tabulated.status == detail::status_ok) {
        detail::finalize_metadata(tabulated);
        return tabulated;
      }
    }
    if ((low_degree || low_degree_evanescent || resolved_evanescent) &&
        (detail::evaluate_legacy_degree_consistent(
            static_cast<detail::real>(n), mi,
            static_cast<detail::real>(theta_radians), tabulated) && (tabulated.evaluation_method = "evaluate_legacy_degree_consistent", true))) {
      detail::finalize_metadata(tabulated);
      return tabulated;
    }
  }
  if ((detail::evaluate_table_legacy(static_cast<detail::real>(n), mi,
                                    static_cast<detail::real>(theta_radians),
                                    tabulated) && (tabulated.evaluation_method = "evaluate_table_legacy", true))) {
    detail::finalize_metadata(tabulated);
    return tabulated;
  }
  detail::EvalResult canonical;
  if ((detail::evaluate_near_table_legacy(
          static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians), canonical) && (canonical.evaluation_method = "evaluate_near_table_legacy", true))) {
    detail::finalize_metadata(canonical);
    return canonical;
  }
  {
    const detail::real local_carrier = static_cast<detail::real>(n) + 0.5L;
    const detail::real local_alpha = static_cast<detail::real>(mi) / local_carrier;
    // The canonical chart is independently normalized and agrees with the
    // direct legacy solution in this moderate-carrier Airy core.  Leaving the
    // interval (100001,1e6) to the outer global-Airy fallback can lose the
    // connection constant by O(1) near the turning point.
    const bool established_large_carrier = local_carrier >= 1000000.0L;
    const bool missing_transition_band =
        local_carrier > 100001.0L && local_carrier < 1000000.0L;
    if (established_large_carrier || missing_transition_band) {
      const auto map = crli::detail::canonical_map(
          static_cast<double>(local_alpha), theta_radians);
      const detail::real scaled_action =
          local_carrier * std::fabs(static_cast<detail::real>(map.action));
      const detail::real action_limit =
          established_large_carrier ? 400.0L : 8.0L;
      if (scaled_action <= action_limit &&
          (detail::evaluate_canonical_legacy(
              static_cast<detail::real>(n), mi,
              static_cast<detail::real>(theta_radians), canonical) && (canonical.evaluation_method = "evaluate_canonical_legacy", true))) {
        const detail::real st = std::sin(
            static_cast<detail::real>(theta_radians));
        const detail::real ct = std::cos(
            static_cast<detail::real>(theta_radians));
        const detail::real coefficient = static_cast<detail::real>(n) *
            (static_cast<detail::real>(n) + 1.0L) -
            static_cast<detail::real>(mi) * mi / (st * st);
        canonical.d2theta_over_value =
            -(ct / st) * canonical.dtheta_over_value - coefficient;
        canonical.d2theta = canonical.value * canonical.d2theta_over_value;
        detail::finalize_metadata(canonical);
        return canonical;
      }
    }
  }
  {
    const detail::real local_carrier = static_cast<detail::real>(n) + 0.5L;
    const detail::real local_alpha = static_cast<detail::real>(mi) / local_carrier;
    if (local_carrier >= 1000000.0L &&
        theta_radians < std::asin(local_alpha) &&
        (detail::evaluate_global_airy_legacy(
            static_cast<detail::real>(n), mi,
            static_cast<detail::real>(theta_radians), canonical) && (canonical.evaluation_method = "evaluate_global_airy_legacy", true))) {
      detail::finalize_metadata(canonical);
      return canonical;
    }
  }
  {
    const detail::real local_carrier = static_cast<detail::real>(n) + 0.5L;
    const detail::real local_alpha = static_cast<detail::real>(mi) /
                                     local_carrier;
    if (local_carrier >= 1000.0L && local_carrier < 10000.0L &&
        theta_radians > std::asin(local_alpha)) {
      const auto map = crli::detail::canonical_map(
          static_cast<double>(local_alpha), theta_radians);
      detail::EvalResult transition;
      if (local_carrier * static_cast<detail::real>(map.action) < 1.0L &&
          (detail::evaluate_canonical_legacy(
              static_cast<detail::real>(n), mi,
              static_cast<detail::real>(theta_radians), transition) && (transition.evaluation_method = "evaluate_canonical_legacy", true))) {
        detail::finalize_metadata(transition);
        return transition;
      }
    }
  }
  if ((detail::evaluate_endpoint_p_fast(
          static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians), canonical) && (canonical.evaluation_method = "evaluate_endpoint_p_fast", true))) {
    detail::finalize_metadata(canonical);
    return canonical;
  }
  if ((detail::evaluate_pole_residual_fast(
          static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians), canonical) && (canonical.evaluation_method = "evaluate_pole_residual_fast", true))) {
    detail::finalize_metadata(canonical);
    return canonical;
  }
  if (mi == 0 && n >= 10000.0 &&
      theta_radians <= 0.1 * static_cast<double>(detail::pi) / 180.0 &&
      (detail::evaluate_canonical_legacy(
          static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians), canonical) && (canonical.evaluation_method = "evaluate_canonical_legacy", true))) {
    detail::finalize_metadata(canonical);
    return canonical;
  }
  const detail::real carrier = static_cast<detail::real>(n) + 0.5L;
  const detail::real alpha = static_cast<detail::real>(mi) / carrier;
  if (carrier >= 1000.0L && carrier < 10000.0L && alpha >= 0.01L) {
    const auto map = crli::detail::canonical_map(
        static_cast<double>(alpha), theta_radians);
    if (carrier * std::fabs(static_cast<detail::real>(map.action)) <= 8.0L &&
        (detail::evaluate_legacy_degree_consistent(
            static_cast<detail::real>(n), mi,
            static_cast<detail::real>(theta_radians), canonical) && (canonical.evaluation_method = "evaluate_legacy_degree_consistent", true))) {
      detail::finalize_metadata(canonical);
      return canonical;
    }
  }
  if (carrier >= 10000.0L && carrier <= 100001.0L && alpha > 0.0L) {
    const auto map = crli::detail::canonical_map(
        static_cast<double>(alpha), theta_radians);
    if (carrier * std::fabs(static_cast<detail::real>(map.action)) <= 17.0L) {
      detail::EvalResult transition = detail::evaluate_legacy(
          static_cast<detail::real>(n), mi,
          static_cast<detail::real>(theta_radians));
      transition.evaluation_method = "evaluate_legacy";
      if (transition.status == detail::status_ok) {
        detail::finalize_metadata(transition);
        return transition;
      }
    }
  }
  bool prefer_wkb = false;
  detail::real dispatch_action = 0.0L;
  if (theta_radians > std::asin(alpha)) {
    const auto map = crli::detail::canonical_map(
        static_cast<double>(alpha), theta_radians);
    dispatch_action = carrier * static_cast<detail::real>(map.action);
    prefer_wkb = dispatch_action >= 40.0L;
  }
  if (prefer_wkb) {
    const bool precise_equator = carrier >= 1000000.0L &&
        theta_radians >= 88.0 * static_cast<double>(detail::pi) / 180.0;
    if ((!precise_equator && carrier >= 1000000.0L && dispatch_action > 400.0L &&
         (detail::evaluate_canonical_legacy(
             static_cast<detail::real>(n), mi,
             static_cast<detail::real>(theta_radians), canonical) && (canonical.evaluation_method = "evaluate_canonical_legacy", true))) ||
        (detail::wkb_evaluate_fast(static_cast<detail::real>(n), mi,
                                  static_cast<detail::real>(theta_radians),
                                  canonical) && (canonical.evaluation_method = "wkb_evaluate_fast", true)) ||
        (detail::wkb_evaluate_legacy(static_cast<detail::real>(n), mi,
                                    static_cast<detail::real>(theta_radians),
                                    canonical) && (canonical.evaluation_method = "wkb_evaluate_legacy", true))) {
      detail::finalize_metadata(canonical);
      return canonical;
    }
  }
  if ((detail::evaluate_global_airy_legacy(static_cast<detail::real>(n), mi,
                                          static_cast<detail::real>(theta_radians),
                                          canonical) && (canonical.evaluation_method = "evaluate_global_airy_legacy", true))) {
    detail::finalize_metadata(canonical);
    return canonical;
  }
  if (!prefer_wkb &&
      ((detail::wkb_evaluate_fast(static_cast<detail::real>(n), mi,
                                 static_cast<detail::real>(theta_radians),
                                 canonical) && (canonical.evaluation_method = "wkb_evaluate_fast", true)) ||
       (detail::wkb_evaluate_legacy(static_cast<detail::real>(n), mi,
                                   static_cast<detail::real>(theta_radians),
                                   canonical) && (canonical.evaluation_method = "wkb_evaluate_legacy", true)))) {
    detail::finalize_metadata(canonical);
    return canonical;
  }
  if ((detail::evaluate_canonical_legacy(static_cast<detail::real>(n), mi,
                                        static_cast<detail::real>(theta_radians),
                                        canonical) && (canonical.evaluation_method = "evaluate_canonical_legacy", true))) {
    detail::finalize_metadata(canonical);
    return canonical;
  }
  detail::EvalResult result =
      detail::evaluate_legacy(static_cast<detail::real>(n), mi,
                              static_cast<detail::real>(theta_radians));
  result.evaluation_method = "evaluate_legacy";
  detail::finalize_metadata(result);
  return result;
}

struct RootBoundary {
  int sign = 0;
  real log_abs = 0.0L;
  bool valid = false;
};

struct RootPhaseModel {
  std::array<real, 3> degree{};
  std::array<real, 3> correction{};
  int points = 0;

  RootBoundary evaluate(int order, double nu_value,
                        double theta_value) noexcept {
    const real nu = static_cast<real>(nu_value);
    const real angle = static_cast<real>(theta_value);
    const real kappa = nu + 0.5L;
    if (!(kappa >= 500000.0L && kappa < 1000000.0L &&
          angle <= 6.0L * pi / 180.0L))
      return {};
    const real eta2 = static_cast<real>(order) * order - 0.25L;
    const real sine = std::sin(angle);
    const real qtheta = kappa * kappa - eta2 / (sine * sine);
    if (!(qtheta > 0.01L * kappa * kappa))
      return {};
    const WkbDensity endpoint = wkb_density(kappa, order, angle);
    const WkbDensity equator = wkb_density(kappa, order, half_pi);
    if (!endpoint.valid || !equator.valid)
      return {};
    const auto integrand = [=](real x) {
      const WkbDensity density = wkb_density(kappa, order, x);
      if (!density.valid)
        return std::numeric_limits<real>::quiet_NaN();
      const real st = std::sin(x);
      return density.p -
             std::sqrt(kappa * kappa - eta2 / (st * st));
    };
    const real radius = 0.6L * pi / angle;
    const bool reusable = points > 0 &&
                          std::fabs(nu - degree[0]) <= radius;
    real phase_correction = 0.0L;
    if (reusable && points == 3) {
      for (int i = 0; i < 3; ++i) {
        real basis = 1.0L;
        for (int j = 0; j < 3; ++j)
          if (i != j)
            basis *= (nu - degree[j]) / (degree[i] - degree[j]);
        phase_correction += correction[i] * basis;
      }
    } else {
      if (!gauss_phase_integral(integrand, angle, half_pi,
                                phase_correction, 1.0e-10L))
        return {};
      if (reusable) {
        degree[points] = nu;
        correction[points] = phase_correction;
        ++points;
      } else {
        degree[0] = nu;
        correction[0] = phase_correction;
        points = 1;
      }
    }
    const real phase = reduced_leading_phase(kappa, order, angle,
                                             phase_correction);
    real value0 = 0.0L, derivative0 = 0.0L;
    equator_legacy_values(nu, order, value0, derivative0);
    const real rootp0 = std::sqrt(equator.p);
    const real c1 = value0 * rootp0;
    const real c2 = derivative0 / rootp0;
    const real shape = c1 * std::cos(phase) + c2 * std::sin(phase);
    if (!finite(shape) || shape == 0.0L)
      return {};
    const real log_abs = std::log(std::fabs(shape)) -
                         0.5L * std::log(endpoint.p) -
                         0.5L * std::log(sine);
    return {shape > 0.0L ? 1 : -1, log_abs, finite(log_abs)};
  }
};

inline RootBoundary root_boundary(int m, double n, double theta,
                                  RootPhaseModel &model) {
  const RootBoundary accelerated = model.evaluate(m, n, theta);
  if (accelerated.valid)
    return accelerated;
  const EvalResult value = evaluate_legacy_schmidt_raw(m, n, theta);
  return {value.sign, value.log_abs_value,
          value.status == status_ok &&
              (value.sign == 0 || finite(value.log_abs_value))};
}

inline real root_action(int m, real kappa, real theta) {
  if (m == 0)
    return kappa * theta;
  const real beta = std::sqrt(static_cast<real>(m) * m - 0.25L);
  const real alpha = beta / kappa;
  if (!(alpha < std::sin(theta)))
    return 0.0L;
  const auto map = crli::detail::canonical_map(
      static_cast<double>(alpha), static_cast<double>(theta));
  return kappa * static_cast<real>(map.action);
}

inline real root_predictor(int m, int root, real theta) {
  const real target = pi * (static_cast<real>(root) - 0.25L);
  if (m == 0)
    return target / theta - 0.5L;
  const real beta = std::sqrt(static_cast<real>(m) * m - 0.25L);
  real lo = beta / std::sin(theta) * (1.0L + 1.0e-15L);
  real hi = lo + target / theta + pi / theta;
  while (hi <= 100000000.5L && root_action(m, hi, theta) < target)
    hi = lo + 2.0L * (hi - lo);
  if (hi > 100000000.5L)
    return std::numeric_limits<real>::quiet_NaN();
  for (int i = 0; i < 80; ++i) {
    const real mid = 0.5L * (lo + hi);
    if (root_action(m, mid, theta) < target)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5L * (lo + hi) - 0.5L;
}

inline double root_secant(int m, double theta, real predictor) {
  const double gap = static_cast<double>(pi) / theta;
  double lo = std::max(static_cast<double>(m),
                       static_cast<double>(predictor) - 0.2 * gap);
  double hi = static_cast<double>(predictor) + 0.2 * gap;
  RootPhaseModel model;
  RootBoundary flo = root_boundary(m, lo, theta, model);
  RootBoundary fhi = root_boundary(m, hi, theta, model);
  for (int expand = 0; expand < 4 &&
       (!flo.valid || !fhi.valid || flo.sign == fhi.sign); ++expand) {
    lo = std::max(static_cast<double>(m), lo - 0.1 * gap);
    hi += 0.1 * gap;
    flo = root_boundary(m, lo, theta, model);
    fhi = root_boundary(m, hi, theta, model);
  }
  if (!flo.valid || !fhi.valid || flo.sign == fhi.sign || hi > 1.0e8)
    return std::numeric_limits<double>::quiet_NaN();
  int previous_side = 0;
  for (int iteration = 0; iteration < 48; ++iteration) {
    if (std::nextafter(lo, hi) == hi)
      break;
    const real scale = std::max(flo.log_abs, fhi.log_abs);
    const real a = flo.sign * std::exp(flo.log_abs - scale);
    const real b = fhi.sign * std::exp(fhi.log_abs - scale);
    double candidate = static_cast<double>(
        (static_cast<real>(lo) * b - static_cast<real>(hi) * a) / (b - a));
    if (!std::isfinite(candidate) || candidate <= lo || candidate >= hi ||
        iteration % 10 == 9)
      candidate = lo + 0.5 * (hi - lo);
    const RootBoundary value = root_boundary(m, candidate, theta, model);
    if (value.valid && value.sign == 0)
      return candidate;
    if (!value.valid) {
      const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
          std::max({1.0, std::fabs(lo), std::fabs(hi)});
      if (hi - lo <= tolerance)
        return candidate;
      return std::numeric_limits<double>::quiet_NaN();
    }
    const int side = value.sign == flo.sign ? -1 : 1;
    if (side < 0) {
      lo = candidate;
      flo = value;
      if (previous_side == side)
        fhi.log_abs -= log2_value;
    } else {
      hi = candidate;
      fhi = value;
      if (previous_side == side)
        flo.log_abs -= log2_value;
    }
    previous_side = side;
  }
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::fabs(lo), std::fabs(hi)});
  return hi - lo <= tolerance ? lo + 0.5 * (hi - lo)
                             : std::numeric_limits<double>::quiet_NaN();
}

struct RootDegreeBoundary {
  real value = 0.0L;
  real derivative = 0.0L;
  int sign = 0;
  bool valid = false;
};

inline RootDegreeBoundary root_degree_boundary(int m, real n, real theta) {
  const BoundaryDegreeData data = boundary_degree_data(n, m, theta);
  if (!data.valid)
    return {};
  const real value = data.center.value;
  return {value, data.value_nu,
          value > 0.0L ? 1 : value < 0.0L ? -1 : 0,
          finite(value) && finite(data.value_nu)};
}

inline double root_newton(int m, real theta, real predictor) {
  const real gap = pi / theta;
  real lo = std::max(static_cast<real>(m), predictor - 0.48L * gap);
  real hi = predictor + 0.48L * gap;
  RootDegreeBoundary flo = root_degree_boundary(m, lo, theta);
  RootDegreeBoundary fhi = root_degree_boundary(m, hi, theta);
  for (int expand = 0; expand < 8 &&
       (!flo.valid || !fhi.valid || flo.sign == fhi.sign); ++expand) {
    lo = std::max(static_cast<real>(m), lo - 0.25L * gap);
    hi += 0.25L * gap;
    flo = root_degree_boundary(m, lo, theta);
    fhi = root_degree_boundary(m, hi, theta);
  }
  if (!flo.valid || !fhi.valid || flo.sign == fhi.sign || hi > 1.0e8L)
    return std::numeric_limits<double>::quiet_NaN();
  real x = RSCHACompat::clamp(predictor, lo, hi);
  real best = x;
  real best_step = std::numeric_limits<real>::infinity();
  for (int iteration = 0; iteration < 12; ++iteration) {
    RootDegreeBoundary value = root_degree_boundary(m, x, theta);
    if (!value.valid) {
      x = 0.5L * (lo + hi);
      value = root_degree_boundary(m, x, theta);
    }
    if (!value.valid)
      return std::numeric_limits<double>::quiet_NaN();
    const real step = value.derivative != 0.0L
                          ? value.value / value.derivative
                          : std::numeric_limits<real>::infinity();
    if (std::fabs(step) < best_step) {
      best = x;
      best_step = std::fabs(step);
    }
    if (best_step <= 16.0L * std::numeric_limits<real>::epsilon() *
                         std::max(1.0L, std::fabs(best)))
      break;
    if (value.sign == flo.sign) {
      lo = x;
      flo = value;
    } else {
      hi = x;
    }
    real candidate = value.derivative != 0.0L ? x - step
                                               : 0.5L * (lo + hi);
    if (!finite(candidate) || candidate <= lo || candidate >= hi)
      candidate = 0.5L * (lo + hi);
    x = candidate;
  }
  return best_step <= 1.0e-10L * std::max(1.0L, std::fabs(best))
      ? static_cast<double>(best)
      : std::numeric_limits<double>::quiet_NaN();
}

}

inline detail::EvalResult EvaluateLegacySchmidt(double m, double n,
                                                double theta_radians) {
  int axis_order = 0;
  if (theta_radians == 0.0 && detail::parse_order(m, axis_order) &&
      std::isfinite(n) && n >= axis_order && n <= 1.0e8) {
    detail::EvalResult axis;
    const detail::real degree = n;
    if (axis_order == 0) {
      axis.value = 1.0L;
      axis.d2theta = -0.5L * degree * (degree + 1.0L);
    } else if (axis_order == 1) {
      axis.dtheta = std::sqrt(0.5L * degree * (degree + 1.0L));
    } else if (axis_order == 2) {
      axis.d2theta = std::sqrt((degree - 1.0L) * degree *
          (degree + 1.0L) * (degree + 2.0L) / 8.0L);
    }
    axis.evaluation_method = "analytic_axis";
    detail::finalize_metadata(axis);
    return axis;
  }
  if (theta_radians > 0.0 && theta_radians * theta_radians == 0.0 &&
      detail::parse_order(m, axis_order) &&
      detail::valid_mode_domain(axis_order, n, theta_radians)) {
    detail::EvalResult tiny;
    const detail::real lambda = static_cast<detail::real>(n) * (n + 1.0L);
    const auto amplitude = detail::pole_log_scale(n, axis_order);
    const auto log_angle = crli::detail::dd_log_positive(
        crli::detail::DoubleDouble{theta_radians,0.0});
    tiny.sign = 1;
    const auto value_log = amplitude + log_angle*static_cast<double>(axis_order);
    tiny.log_abs_value = value_log.hi;
    tiny.log_abs_value_lo = value_log.lo;
    tiny.value = detail::safe_exp(tiny.log_abs_value + tiny.log_abs_value_lo);
    const detail::real first = axis_order == 0 ? -0.5L * lambda : axis_order;
    const detail::real second = axis_order == 0 ? -0.5L * lambda :
        axis_order == 1 ? (2.0L - 3.0L * lambda) / 4.0L :
        static_cast<detail::real>(axis_order) * (axis_order - 1.0L);
    if (first != 0.0L) {
      tiny.dtheta_sign = first > 0.0L ? 1 : -1;
      const auto first_log = amplitude + log_angle *
          static_cast<double>(axis_order == 0 ? 1 : axis_order - 1) +
          crli::detail::DoubleDouble{static_cast<double>(std::log(std::fabs(first))),0.0};
      tiny.log_abs_dtheta = first_log.hi;
      tiny.log_abs_dtheta_lo = first_log.lo;
      tiny.dtheta = tiny.dtheta_sign * detail::safe_exp(first_log.hi + first_log.lo);
    }
    if (second != 0.0L) {
      tiny.d2theta_sign = second > 0.0L ? 1 : -1;
      const auto second_log = amplitude + log_angle *
          static_cast<double>(axis_order < 2 ? axis_order : axis_order - 2) +
          crli::detail::DoubleDouble{static_cast<double>(std::log(std::fabs(second))),0.0};
      tiny.log_abs_d2theta = second_log.hi;
      tiny.log_abs_d2theta_lo = second_log.lo;
      tiny.d2theta = tiny.d2theta_sign * detail::safe_exp(second_log.hi + second_log.lo);
    }
    tiny.dtheta_over_value = axis_order == 0 ? first * theta_radians : axis_order / theta_radians;
    tiny.d2theta_over_value = axis_order < 2 ? second :
        (static_cast<detail::real>(axis_order) * (axis_order - 1.0L) / theta_radians) / theta_radians;
    tiny.evaluation_method = "analytic_small_axis";
    detail::finalize_metadata(tiny);
    return tiny;
  }
  detail::EvalResult result =
      detail::evaluate_legacy_schmidt_raw(m, n, theta_radians);
  if (result.status != detail::status_ok)
    return result;
  // Raw polar paths may retain a representable derivative only in log form.
  if (result.dtheta == 0.0L && result.dtheta_sign != 0 &&
      detail::finite(result.log_abs_dtheta)) {
    result.dtheta = result.dtheta_sign * detail::safe_exp(
        result.log_abs_dtheta + result.log_abs_dtheta_lo);
    if (result.dtheta != 0.0L)
      result.log_abs_dtheta_lo = 0.0L;
  }
  const detail::real theta = static_cast<detail::real>(theta_radians);
  const detail::real degree = static_cast<detail::real>(n);
  if ((degree + 0.5L) * std::min(theta, detail::pi - theta) <= 2.0L) {
    if (result.d2theta == 0.0L && result.d2theta_sign != 0 &&
        detail::finite(result.log_abs_d2theta))
      result.d2theta = result.d2theta_sign * detail::safe_exp(
          result.log_abs_d2theta + result.log_abs_d2theta_lo);
    if (!detail::finite(result.dtheta) || !detail::finite(result.d2theta))
      result.status = detail::status_range;
    detail::finalize_metadata(result);
    return result;
  }
  const detail::real order = static_cast<detail::real>(std::llround(m));
  const detail::real st = std::sin(theta);
  const detail::real ct = std::cos(theta);
  const detail::real coefficient =
      degree * (degree + 1.0L) - order * order / (st * st);
  if (result.value == 0.0L && result.sign != 0 &&
      detail::finite(result.log_abs_value) &&
      detail::finite(result.dtheta_over_value)) {
    // The ordinary value may underflow while its derivatives remain in
    // range.  Form the ODE ratio before independently materializing d2.
    const detail::real ratio = -(ct / st) * result.dtheta_over_value - coefficient;
    if (!detail::finite(ratio)) {
      result.status = detail::status_range;
      return result;
    }
    result.d2theta_over_value = ratio;
    result.d2theta_sign = result.sign * ((ratio > 0.0L) - (ratio < 0.0L));
    result.log_abs_d2theta = ratio == 0.0L
        ? -std::numeric_limits<detail::real>::infinity()
        : result.log_abs_value + std::log(std::fabs(ratio));
    result.log_abs_d2theta_lo = result.log_abs_value_lo;
    result.d2theta = result.d2theta_sign * detail::safe_exp(
        result.log_abs_d2theta + result.log_abs_d2theta_lo);
    if (result.d2theta != 0.0L)
      result.log_abs_d2theta_lo = 0.0L;
  } else {
    result.d2theta = -(ct / st) * result.dtheta - coefficient * result.value;
    result.d2theta_over_value = result.value != 0.0L
        ? result.d2theta / result.value
        : std::numeric_limits<detail::real>::quiet_NaN();
  }
  if (!detail::finite(result.d2theta)) {
    result.status = detail::status_range;
    return result;
  }
  detail::finalize_metadata(result);
  return result;
}

inline bool PrepareLegacySchmidtMode(double m, double n) {
  int mi = 0;
  return detail::parse_order(m, mi) && std::isfinite(n) &&
         n >= static_cast<double>(mi) && n <= 1.0e8;
}

inline double FindDirichletDegree(int m, int root, double theta0_radians) {
  if (m < 0 || m > 100000000 || root < 1 ||
      !(theta0_radians > 0.0 &&
        theta0_radians <= static_cast<double>(detail::half_pi)))
    return std::numeric_limits<double>::quiet_NaN();
  const detail::real theta = static_cast<detail::real>(theta0_radians);
  const detail::real predictor = detail::root_predictor(m, root, theta);
  if (!detail::finite(predictor) || predictor < static_cast<detail::real>(m) ||
      predictor > 1.0e8L)
    return std::numeric_limits<double>::quiet_NaN();
  if (m <= 1 && root >= 100 && theta <= 2.0L * detail::pi / 180.0L)
    return detail::root_newton(m, theta, predictor);
  return detail::root_secant(m, theta0_radians, predictor);
}

inline double CapNormLegacySchmidtTurningIVP(double m, double n,
                                             double theta0_radians) {
  int mi = 0;
  if (!detail::parse_order(m, mi) ||
      !detail::valid_mode_domain(mi, n, theta0_radians))
    return std::numeric_limits<double>::quiet_NaN();
  return detail::to_double(detail::cap_norm_turning_ivp_impl(
      mi, static_cast<detail::real>(n),
      static_cast<detail::real>(theta0_radians)));
}

inline double CapNormLegacySchmidt(double m, double n, double theta0_radians) {
  int mi = 0;
  if (!detail::parse_order(m, mi) ||
      !detail::valid_mode_domain(mi, n, theta0_radians))
    return std::numeric_limits<double>::quiet_NaN();
  return detail::to_double(
      detail::cap_norm_legacy_impl(mi, static_cast<detail::real>(n),
                                   static_cast<detail::real>(theta0_radians)));
}

inline double BoundaryDegreeDerivativeLegacySchmidt(double m, double n,
                                                    double theta0_radians) {
  int mi = 0;
  if (!detail::parse_order(m, mi) ||
      !detail::valid_mode_domain(mi, n, theta0_radians))
    return std::numeric_limits<double>::quiet_NaN();
#ifndef RSCHA_CR_TABLE_BUILD
  double tabulated = 0.0;
  if (CRTable::legendre_degree(mi, n, theta0_radians, tabulated))
    return tabulated;
#endif
  return detail::to_double(detail::boundary_degree_derivative_legacy_impl(
      mi, static_cast<detail::real>(n),
      static_cast<detail::real>(theta0_radians)));
}

}


#pragma float_control(pop)
