#pragma once

#include <algorithm>
#include <cmath>

#pragma float_control(precise, on, push)
namespace CRNear {

struct LogScalar {
  int sign = 0;
  long double log_abs = 0.0L;
  bool ok = false;
  long double log_abs_lo = 0.0L;
};

struct DoubleDouble {
  double hi = 0.0, lo = 0.0;
  explicit DoubleDouble(double value = 0.0) : hi(value) {}
  double value() const { return hi + lo; }
  void add(const DoubleDouble &other) {
    const double sum = hi + other.hi;
    const double error = std::abs(hi) >= std::abs(other.hi)
        ? (hi - sum) + other.hi : (other.hi - sum) + hi;
    const double tail = lo + other.lo + error;
    hi = sum + tail;
    lo = tail - (hi - sum);
  }
  void multiply(double factor) {
    const double product = hi * factor;
    const double error = std::fma(hi, factor, -product) + lo * factor;
    hi = product + error;
    lo = error - (hi - product);
  }
};

inline LogScalar endpoint(double degree, double order, double offset,
                          double theta) {
  const long double degree_ld = degree;
  const long double order_ld = order;
  const long double offset_ld = offset;
  const long double sine = std::sin(static_cast<long double>(theta));
  const long double z = std::pow(std::sin(0.5L * theta), 2);
  DoubleDouble term(1.0), sum(1.0);
  bool converged = false;
  for (int k = 0; k < 20000; ++k) {
    const long double index = static_cast<long double>(k);
    const double ratio = static_cast<double>(
        ((index - offset_ld) * (index + order_ld + degree_ld + 1.0L) /
         ((index + order_ld + 1.0L) * (index + 1.0L))) * z);
    term.multiply(ratio);
    sum.add(term);
    if (std::abs(term.value()) < 1.0e-28 * std::abs(sum.value())) {
      converged = true;
      break;
    }
  }
  const long double series = sum.value();
  if (!converged || series == 0.0 || !std::isfinite(series))
    return {};
  const long double a = degree_ld + order_ld + 1.0L;
  const long double b = order_ld + 1.0L;
  const long double inverse_a = 1.0L / a, inverse_b = 1.0L / b;
  const long double correction_a = inverse_a / 12.0L -
      inverse_a * inverse_a * inverse_a / 360.0L;
  const long double correction_b = inverse_b / 12.0L -
      inverse_b * inverse_b * inverse_b / 360.0L;
  using crli::detail::dd_log_positive;
  const auto ratio_log = dd_log_positive(
      crli::detail::DoubleDouble{2.0,0.0}+
      crli::detail::DoubleDouble{static_cast<double>(offset_ld-1.0L),0.0}/
      crli::detail::DoubleDouble{static_cast<double>(b),0.0});
  const auto balanced_log = dd_log_positive(
      crli::detail::DoubleDouble{static_cast<double>(sine),0.0})+
      ratio_log-dd_log_positive(crli::detail::DoubleDouble{2.0,0.0});
  DoubleDouble result_sum;
  auto accumulate = [&](long double value, long double factor = 1.0L) {
    DoubleDouble part(static_cast<double>(value));
    part.multiply(static_cast<double>(factor));
    result_sum.add(part);
  };
  accumulate(std::log(b), 0.5L * offset_ld - 0.25L);
  DoubleDouble balanced_part;
  balanced_part.hi=balanced_log.hi;balanced_part.lo=balanced_log.lo;
  balanced_part.multiply(static_cast<double>(order_ld));result_sum.add(balanced_part);
  DoubleDouble ratio_part;
  ratio_part.hi=ratio_log.hi;ratio_part.lo=ratio_log.lo;
  ratio_part.multiply(static_cast<double>(0.5L*offset_ld+0.25L));result_sum.add(ratio_part);
  accumulate(-0.5L * offset_ld + 0.5L -
             0.25L * std::log(2.0L * static_cast<long double>(crli::pi)));
  accumulate(correction_a, 0.5L);
  accumulate(correction_b, -1.0L);
  accumulate(std::log(degree_ld + 0.5L) -
             std::lgamma(offset_ld + 1.0L) + std::log(sine), 0.5L);
  accumulate(std::log(std::abs(series)));
  const long double result = result_sum.value();
  return {series < 0.0 ? -1 : 1, result_sum.hi, std::isfinite(result), result_sum.lo};
}

inline LogScalar recessive_log(double degree, int order, double theta) {
  const double gap = degree - order;
  if (!(degree >= order && gap <= 450000.0))
    return {};
  if (gap <= 4.0)
    return endpoint(degree, order, gap, theta);
  const int steps = static_cast<int>(std::floor(gap));
  const double fractional = gap - steps;
  const double lower = degree - steps;
  const auto p0 = endpoint(lower, order, fractional, theta);
  const auto p1 = endpoint(lower + 1.0, order, fractional + 1.0, theta);
  if (!p0.ok || !p1.ok || p0.sign == 0 || p1.sign == 0)
    return {};
  const long double sine = std::sin(static_cast<long double>(theta));
  using crli::detail::dd_log_positive;
  using DD=crli::detail::DoubleDouble;
  const auto log_sine=dd_log_positive(DD{static_cast<double>(sine),0.0});
  const auto log0=DD{static_cast<double>(p0.log_abs),static_cast<double>(p0.log_abs_lo)}-
      (dd_log_positive(DD{lower+0.5,0.0})+log_sine)*0.5;
  const auto log1=DD{static_cast<double>(p1.log_abs),static_cast<double>(p1.log_abs_lo)}-
      (dd_log_positive(DD{lower+1.5,0.0})+log_sine)*0.5;
  auto common=log0.hi>log1.hi?log0:log1;
  const auto delta0=log0-common,delta1=log1-common;
  DoubleDouble previous(p0.sign*std::exp(delta0.hi+delta0.lo));
  DoubleDouble value(p1.sign*std::exp(delta1.hi+delta1.lo));
  for (int index = 1; index < steps; ++index) {
    const long double local_degree = lower + index;
    const long double minus = fractional + index;
    const long double plus = local_degree + order;
    const long double denominator = (minus + 1.0L) * (plus + 1.0L);
    DoubleDouble next = value;
    next.multiply(static_cast<double>(
        (2.0L * local_degree + 1.0L) * std::cos(theta) /
        std::sqrt(denominator)));
    DoubleDouble backward = previous;
    backward.multiply(static_cast<double>(-
        std::sqrt(plus * minus / denominator)));
    next.add(backward);
    previous = value;
    value = next;
    const double magnitude = std::max(std::abs(previous.value()),
                                      std::abs(value.value()));
    if (!(magnitude > 0.0) || !std::isfinite(magnitude))
      return {};
    if (magnitude > 1.0e100L || magnitude < 1.0e-100L) {
      previous.multiply(1.0 / magnitude);
      value.multiply(1.0 / magnitude);
      common = common+dd_log_positive(DD{magnitude,0.0});
    }
  }
  if (!(value.value() > 0.0))
    return {};
  const auto result=common+dd_log_positive(DD{value.hi,value.lo})+
      (dd_log_positive(DD{degree+0.5,0.0})+log_sine)*0.5;
  return {1,result.hi,std::isfinite(result.hi),result.lo};
}

}
#pragma float_control(pop)
