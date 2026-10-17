#pragma once

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

// Real Schmidt semi-normalized harmonics, without the Condon--Shortley
// phase: P_nm(cos(theta)) * (g_nm cos(m phi) + h_nm sin(m phi)).
// Passive convention: f_cap(u) = f_geo(Rz(alpha) Ry(beta) Rz(gamma) u).
namespace SHRotation {

inline std::size_t triangular_size(int degree) {
  if (degree < 0 || degree > 10000)
    throw std::invalid_argument("SH rotation degree must be in [0,10000].");
  return (static_cast<std::size_t>(degree) + 1) * (degree + 2) / 2;
}

inline std::size_t index(int n, int m) {
  return static_cast<std::size_t>(n) * (n + 1) / 2 + m;
}

// Build d^(j) by coupling d^(j-1/2) with a spin-1/2 rotation.
// This recurrence uses bounded Clebsch--Gordan factors; no factorial sum,
// grid synthesis, quadrature, or coupling between different degrees.
// Only two adjacent matrices are retained. Total work is O(N^3).
inline void rotate_euler(int degree, double alpha_deg, double beta_deg,
                         double gamma_deg, std::vector<double>& g,
                         std::vector<double>& h, bool inverse = false) {
  const auto count = triangular_size(degree);
  if (g.size() != count || h.size() != count)
    throw std::invalid_argument("SH coefficients must be a complete triangular table.");
  for (const double angle : {alpha_deg, beta_deg, gamma_deg})
    if (!std::isfinite(angle))
      throw std::invalid_argument("SH rotation angles must be finite.");
  for (std::size_t i = 0; i < count; ++i)
    if (!std::isfinite(g[i]) || !std::isfinite(h[i]))
      throw std::invalid_argument("SH coefficients must be finite.");

  if (inverse) {
    const double old_alpha = alpha_deg;
    alpha_deg = -gamma_deg;
    beta_deg = -beta_deg;
    gamma_deg = -old_alpha;
  }
  const long double rad = std::acos(-1.0L) / 180.0L;
  const long double alpha = std::remainder(alpha_deg, 360.0) * rad;
  const long double beta = std::remainder(beta_deg, 720.0) * rad;
  const long double gamma = std::remainder(gamma_deg, 360.0) * rad;
  const long double c = std::cos(beta / 2);
  const long double s = std::sin(beta / 2);
  const long double invsqrt2 = 1 / std::sqrt(2.0L);
  std::vector<long double> previous(1, 1);

  for (int twice_j = 1; twice_j <= 2 * degree; ++twice_j) {
    const int width = twice_j + 1;
    const int old_width = twice_j;
    std::vector<long double> current(static_cast<std::size_t>(width) * width);
    std::vector<long double> minus(width), plus(width);
    for (int i = 0; i < width; ++i) {
      minus[i] = std::sqrt(static_cast<long double>(twice_j - i) / twice_j);
      plus[i] = std::sqrt(static_cast<long double>(i) / twice_j);
    }
    for (int row = 0; row < width; ++row) {
      for (int col = 0; col < width; ++col) {
        long double value = 0;
        if (row < old_width && col < old_width)
          value += c * minus[row] * minus[col] * previous[row * old_width + col];
        if (row < old_width && col > 0)
          value -= s * minus[row] * plus[col] * previous[row * old_width + col - 1];
        if (row > 0 && col < old_width)
          value += s * plus[row] * minus[col] * previous[(row - 1) * old_width + col];
        if (row > 0 && col > 0)
          value += c * plus[row] * plus[col] * previous[(row - 1) * old_width + col - 1];
        current[row * width + col] = value;
      }
    }
    previous.swap(current);
    if (twice_j % 2 != 0)
      continue;

    const int n = twice_j / 2;
    using Complex = std::complex<long double>;
    std::vector<Complex> input(width), output(width);
    input[n] = g[index(n, 0)];
    for (int m = 1; m <= n; ++m) {
      const long double sign = m % 2 ? -1 : 1;
      const auto ix = index(n, m);
      input[n - m] = Complex(g[ix], h[ix]) * invsqrt2;
      input[n + m] = sign * Complex(g[ix], -h[ix]) * invsqrt2;
    }
    for (int m = -n; m <= n; ++m)
      input[n + m] *= std::polar(1.0L, m * alpha);
    for (int row = 0; row < width; ++row) {
      Complex sum = 0;
      for (int col = 0; col < width; ++col)
        sum += previous[row * width + col] * input[col];
      output[row] = sum * std::polar(1.0L, (row - n) * gamma);
    }
    g[index(n, 0)] = static_cast<double>(output[n].real());
    h[index(n, 0)] = 0;
    for (int m = 1; m <= n; ++m) {
      const long double sign = m % 2 ? -1 : 1;
      g[index(n, m)] = static_cast<double>(
          (output[n - m].real() + sign * output[n + m].real()) * invsqrt2);
      h[index(n, m)] = static_cast<double>(
          (output[n - m].imag() - sign * output[n + m].imag()) * invsqrt2);
    }
  }
  h[0] = 0;
}

inline void geographic_to_cap(double longitude, double latitude, int degree,
                               std::vector<double>& g, std::vector<double>& h,
                               bool inverse = false, double meridian = 0) {
  if (!std::isfinite(longitude) || !std::isfinite(latitude) ||
      longitude < -180 || longitude > 180 || latitude < -90 || latitude > 90)
    throw std::invalid_argument("Invalid spherical cap centre.");
  rotate_euler(degree, longitude, 90 - latitude, meridian, g, h, inverse);
}
} // namespace SHRotation
