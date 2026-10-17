#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include "RSCHACompatibility.h"
#include <fstream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace CRTable {

enum class Kind : std::uint32_t { legendre = 1, mehler = 2, constant = 3 };

struct Cell {
  std::uint32_t branch = 0;
  std::array<std::uint16_t, 3> degree{};
  std::array<double, 3> lo{};
  std::array<double, 3> hi{};
  std::vector<double> coefficients;
};

class Atlas {
public:
  bool ready() const noexcept { return ready_; }

  bool evaluate_into(std::uint32_t branch, const std::array<double, 3> &point,
                     double *value, double *dlast = nullptr,
                     double *d2last = nullptr,
                     int derivative_axis = -1) const noexcept {
    if (!ready_ || !value)
      return false;
    const Cell *cell = locate(branch, point);
    if (!cell)
      return false;

    std::array<double, 3> x{};
    std::array<double, 3> derivative_scale{1.0, 1.0, 1.0};
    for (std::uint32_t axis = 0; axis < dimensions_; ++axis) {
      const double width = cell->hi[axis] - cell->lo[axis];
      x[axis] =
          (2.0 * point[axis] - cell->lo[axis] - cell->hi[axis]) / width;
      derivative_scale[axis] = 2.0 / width;
    }
    if (derivative_axis < 0)
      derivative_axis = static_cast<int>(dimensions_) - 1;
    if (derivative_axis < 0 ||
        derivative_axis >= static_cast<int>(dimensions_))
      return false;

    const std::size_t nx = static_cast<std::size_t>(cell->degree[0]) + 1;
    const std::size_t ny = dimensions_ > 1
                               ? static_cast<std::size_t>(cell->degree[1]) + 1
                               : 1;
    const std::size_t nz = dimensions_ > 2
                               ? static_cast<std::size_t>(cell->degree[2]) + 1
                               : 1;
    const std::size_t stride = nx * ny * nz;
    if (stride > workspace_capacity)
      return false;

    for (std::uint32_t field = 0; field < fields_; ++field) {
      const double *source = cell->coefficients.data() +
                             static_cast<std::size_t>(field) * stride;
      double v = 0.0, d1 = 0.0, d2 = 0.0;
      if (derivative_axis == 2) {
        std::array<double, workspace_capacity> yz;
        std::array<double, workspace_capacity> z;
        for (std::size_t k = 0; k < nz; ++k)
          for (std::size_t j = 0; j < ny; ++j)
            yz[j + ny * k] =
                clenshaw(source + nx * (j + ny * k), nx, 1, x[0]);
        for (std::size_t k = 0; k < nz; ++k)
          z[k] = clenshaw(yz.data() + ny * k, ny, 1, x[1]);
        clenshaw_derivatives(z.data(), nz, 1, x[2], v, d1, d2);
      } else if (derivative_axis == 1) {
        std::array<double, workspace_capacity> yz;
        std::array<double, workspace_capacity> y;
        for (std::size_t k = 0; k < nz; ++k)
          for (std::size_t j = 0; j < ny; ++j)
            yz[j + ny * k] =
                clenshaw(source + nx * (j + ny * k), nx, 1, x[0]);
        for (std::size_t j = 0; j < ny; ++j)
          y[j] = clenshaw(yz.data() + j, nz, ny, x[2]);
        clenshaw_derivatives(y.data(), ny, 1, x[1], v, d1, d2);
      } else {
        std::array<double, workspace_capacity> xz;
        std::array<double, workspace_capacity> xv;
        for (std::size_t k = 0; k < nz; ++k)
          for (std::size_t i = 0; i < nx; ++i)
            xz[i + nx * k] =
                clenshaw(source + i + nx * ny * k, ny, nx, x[1]);
        for (std::size_t i = 0; i < nx; ++i)
          xv[i] = clenshaw(xz.data() + i, nz, nx, x[2]);
        clenshaw_derivatives(xv.data(), nx, 1, x[0], v, d1, d2);
      }
      const double scale = derivative_scale[derivative_axis];
      value[field] = v;
      if (dlast)
        dlast[field] = d1 * scale;
      if (d2last)
        d2last[field] = d2 * scale * scale;
    }
    return true;
  }

  bool evaluate(std::uint32_t branch, const std::array<double, 3> &point,
                std::vector<double> &value, std::vector<double> *dlast = nullptr,
                std::vector<double> *d2last = nullptr,
                int derivative_axis = -1) const noexcept {
    value.assign(fields_, 0.0);
    if (dlast)
      dlast->assign(fields_, 0.0);
    if (d2last)
      d2last->assign(fields_, 0.0);
    return evaluate_into(branch, point, value.data(),
                         dlast ? dlast->data() : nullptr,
                         d2last ? d2last->data() : nullptr, derivative_axis);
  }

  bool load(const std::string &file, Kind expected,
            const std::string &section = {}) noexcept {
    try {
      std::ifstream in(file, std::ios::binary);
      if (!in)
        return false;
      char magic[8]{};
      in.read(magic, sizeof(magic));
      const std::string signature(magic, magic + 7);
      if (signature == "CRPACK1") {
        std::uint32_t version = 0, count = 0;
        in.read(reinterpret_cast<char *>(&version), sizeof(version));
        in.read(reinterpret_cast<char *>(&count), sizeof(count));
        if (!in || version != 1 || count == 0 || count > 64 || section.empty())
          return false;
        std::uint64_t selected_offset = 0, selected_size = 0,
                      selected_checksum = 0;
        for (std::uint32_t index = 0; index < count; ++index) {
          char name[16]{};
          std::uint64_t offset = 0, size = 0, digest = 0;
          in.read(name, sizeof(name));
          in.read(reinterpret_cast<char *>(&offset), sizeof(offset));
          in.read(reinterpret_cast<char *>(&size), sizeof(size));
          in.read(reinterpret_cast<char *>(&digest), sizeof(digest));
          if (!in)
            return false;
          if (section ==
              std::string(name, std::find(name, name + sizeof(name), '\0'))) {
            selected_offset = offset;
            selected_size = size;
            selected_checksum = digest;
          }
        }
        if (selected_size == 0 || selected_size > 2000000000ULL)
          return false;
        in.seekg(0, std::ios::end);
        const std::uint64_t file_size =
            static_cast<std::uint64_t>(in.tellg());
        if (selected_offset > file_size ||
            selected_size > file_size - selected_offset)
          return false;
        in.seekg(static_cast<std::streamoff>(selected_offset));
        std::string payload(static_cast<std::size_t>(selected_size), '\0');
        in.read((payload.empty()?nullptr:&payload[0]), static_cast<std::streamsize>(selected_size));
        if (!in)
          return false;
        std::uint64_t digest = 1469598103934665603ULL;
        for (const unsigned char byte : payload) {
          digest ^= byte;
          digest *= 1099511628211ULL;
        }
        if (digest != selected_checksum)
          return false;
        std::istringstream table(payload, std::ios::binary | std::ios::in);
        return load_table(table, expected);
      }
      if (!section.empty())
        return false;
      in.clear();
      in.seekg(0);
      return load_table(in, expected);
    } catch (...) {
      return false;
    }
  }

private:
  static constexpr std::size_t workspace_capacity = 8192;

  static double clenshaw(const double *coefficient, std::size_t count,
                         std::size_t step, double x) noexcept {
    double b1 = 0.0, b2 = 0.0;
    for (std::size_t k = count; k-- > 1;) {
      const double b0 =
          std::fma(2.0 * x, b1, coefficient[k * step] - b2);
      b2 = b1;
      b1 = b0;
    }
    return std::fma(x, b1, coefficient[0] - b2);
  }

  static void clenshaw_derivatives(const double *coefficient,
                                   std::size_t count, std::size_t step,
                                   double x, double &value, double &d1,
                                   double &d2) noexcept {
    double b1 = 0.0, b2 = 0.0;
    double db1 = 0.0, db2 = 0.0;
    double d2b1 = 0.0, d2b2 = 0.0;
    for (std::size_t k = count; k-- > 1;) {
      const double b0 =
          std::fma(2.0 * x, b1, coefficient[k * step] - b2);
      const double db0 = 2.0 * b1 + 2.0 * x * db1 - db2;
      const double d2b0 = 4.0 * db1 + 2.0 * x * d2b1 - d2b2;
      b2 = b1;
      b1 = b0;
      db2 = db1;
      db1 = db0;
      d2b2 = d2b1;
      d2b1 = d2b0;
    }
    value = std::fma(x, b1, coefficient[0] - b2);
    d1 = b1 + x * db1 - db2;
    d2 = 2.0 * db1 + x * d2b1 - d2b2;
  }

  bool contains(const Cell &cell,
                const std::array<double, 3> &point) const noexcept {
    for (std::uint32_t axis = 0; axis < dimensions_; ++axis) {
      const double tol = 32.0 * std::numeric_limits<double>::epsilon() *
                         (1.0 + std::fabs(point[axis]));
      if (point[axis] < cell.lo[axis] - tol ||
          point[axis] > cell.hi[axis] + tol)
        return false;
    }
    return true;
  }

  const Cell *locate(std::uint32_t branch,
                     const std::array<double, 3> &point) const noexcept {
    struct Cache {
      const Atlas *atlas = nullptr;
      const Cell *cell = nullptr;
      std::uint32_t branch = 0;
    };
    static thread_local Cache cache;
    if (cache.atlas == this && cache.branch == branch && cache.cell &&
        contains(*cache.cell, point))
      return cache.cell;
    const auto range = ranges_.find(branch);
    if (range == ranges_.end())
      return nullptr;
    for (std::uint32_t axis = 0; axis < dimensions_; ++axis)
      if (point[axis] < range->second.lo[axis] ||
          point[axis] > range->second.hi[axis])
        return nullptr;
    if (dimensions_ == 3) {
      const auto begin = cells_.begin() + range->second.first;
      const auto end = cells_.begin() + range->second.second;
      const double tolerance = 32.0 * std::numeric_limits<double>::epsilon() *
                               (1.0 + std::fabs(point[0]));
      const auto xbegin = std::lower_bound(
          begin, end, point[0] - tolerance,
          [](const Cell &cell, double value) { return cell.hi[0] < value; });
      if (xbegin != end && point[0] >= xbegin->lo[0] - tolerance) {
        auto xend = xbegin;
        while (xend != end && xend->lo[0] == xbegin->lo[0] &&
               xend->hi[0] == xbegin->hi[0])
          ++xend;
        const auto ybegin = std::lower_bound(
            xbegin, xend, point[1],
            [](const Cell &cell, double value) { return cell.hi[1] < value; });
        if (ybegin != xend) {
          auto yend = ybegin;
          while (yend != xend && yend->lo[1] == ybegin->lo[1] &&
                 yend->hi[1] == ybegin->hi[1])
            ++yend;
          for (auto cell = ybegin; cell != yend; ++cell)
            if (contains(*cell, point)) {
              cache = {this, &*cell, branch};
              return cache.cell;
            }
        }
      }
    }
    for (std::size_t index = range->second.first; index < range->second.second;
         ++index) {
      if (contains(cells_[index], point)) {
        cache = {this, &cells_[index], branch};
        return cache.cell;
      }
    }
    return nullptr;
  }

  bool load_table(std::istream &in, Kind expected) {
      char magic[8]{};
      in.read(magic, sizeof(magic));
      const std::string signature(magic, magic + 7);
      std::uint32_t version = 0, kind = 0, count = 0;
      in.read(reinterpret_cast<char *>(&version), sizeof(version));
      in.read(reinterpret_cast<char *>(&kind), sizeof(kind));
      in.read(reinterpret_cast<char *>(&dimensions_), sizeof(dimensions_));
      in.read(reinterpret_cast<char *>(&fields_), sizeof(fields_));
      in.read(reinterpret_cast<char *>(&count), sizeof(count));
      if (!in || signature != "CRATBL1" || version != 1 ||
          kind != static_cast<std::uint32_t>(expected) ||
          dimensions_ < 1 || dimensions_ > 3 || fields_ == 0)
        return false;
      cells_.clear();
      cells_.reserve(count);
      for (std::uint32_t index = 0; index < count; ++index) {
        Cell cell;
        in.read(reinterpret_cast<char *>(&cell.branch), sizeof(cell.branch));
        in.read(reinterpret_cast<char *>(cell.degree.data()),
                sizeof(std::uint16_t) * cell.degree.size());
        std::uint16_t reserved = 0;
        in.read(reinterpret_cast<char *>(&reserved), sizeof(reserved));
        in.read(reinterpret_cast<char *>(cell.lo.data()), sizeof(double) * 3);
        in.read(reinterpret_cast<char *>(cell.hi.data()), sizeof(double) * 3);
        std::uint64_t coefficient_count = 0;
        in.read(reinterpret_cast<char *>(&coefficient_count),
                sizeof(coefficient_count));
        std::uint64_t expected_count = fields_;
        for (std::uint32_t axis = 0; axis < dimensions_; ++axis)
          expected_count *= static_cast<std::uint64_t>(cell.degree[axis]) + 1;
        if (!in || coefficient_count != expected_count ||
            coefficient_count > 100000000ULL)
          return false;
        cell.coefficients.resize(static_cast<std::size_t>(coefficient_count));
        in.read(reinterpret_cast<char *>(cell.coefficients.data()),
                static_cast<std::streamsize>(sizeof(double) * coefficient_count));
        if (!in)
          return false;
        cells_.push_back(std::move(cell));
      }
      ranges_.clear();
      for (std::size_t index = 0; index < cells_.size();) {
        const std::uint32_t branch = cells_[index].branch;
        std::size_t end = index + 1;
        while (end < cells_.size() && cells_[end].branch == branch)
          ++end;
        Range range;
        range.first = index;
        range.second = end;
        range.lo.fill(std::numeric_limits<double>::infinity());
        range.hi.fill(-std::numeric_limits<double>::infinity());
        for (std::size_t cell_index = index; cell_index < end; ++cell_index)
          for (std::uint32_t axis = 0; axis < dimensions_; ++axis) {
            range.lo[axis] =
                std::min(range.lo[axis], cells_[cell_index].lo[axis]);
            range.hi[axis] =
                std::max(range.hi[axis], cells_[cell_index].hi[axis]);
          }
        ranges_[branch] = range;
        index = end;
      }
      ready_ = true;
      return true;
  }

  std::uint32_t dimensions_ = 0;
  std::uint32_t fields_ = 0;
  struct Range {
    std::size_t first = 0, second = 0;
    std::array<double, 3> lo{}, hi{};
  };
  std::vector<Cell> cells_;
  std::unordered_map<std::uint32_t, Range> ranges_;
  bool ready_ = false;
};

inline std::vector<RSCHACompat::path> candidate_paths(const char *name) {
  std::vector<RSCHACompat::path> paths;
#ifdef _WIN32
  wchar_t *directory = nullptr;
  std::size_t length = 0;
  if (_wdupenv_s(&directory, &length, L"RSCHA_CR_TABLE_DIR") == 0 && directory) {
    paths.emplace_back(RSCHACompat::path(directory) / name);
    std::free(directory);
  }
#else
  if (const char *directory = std::getenv("RSCHA_CR_TABLE_DIR"))
    paths.emplace_back(RSCHACompat::path(directory) / name);
#endif
  const RSCHACompat::path current = RSCHACompat::current_path();
  // Each deployable directory carries its own table, next to these headers.
  paths.emplace_back(current / name);
  paths.emplace_back(RSCHACompat::path(__FILE__).parent_path() / name);
  return paths;
}

inline const Atlas &atlas(const char *section, Kind kind) {
  struct Holder {
    Atlas data;
    std::once_flag once;
  };
  static Holder legendre_low, legendre_airy, legendre_outer,
      legendre_airy_high, legendre_outer_high, legendre_near, mehler,
      legendre_outer_patch, constant;
  Holder *holder = nullptr;
  const std::string name(section);
  if (name == "crl") holder = &legendre_low;
  else if (name == "crg") holder = &legendre_airy;
  else if (name == "cro") holder = &legendre_outer;
  else if (name == "crh") holder = &legendre_airy_high;
  else if (name == "crx") holder = &legendre_outer_high;
  else if (name == "crp") holder = &legendre_outer_patch;
  else if (name == "crn") holder = &legendre_near;
  else if (name == "crm") holder = &mehler;
  else if (name == "crc") holder = &constant;
  if (!holder)
    return legendre_low.data;
  std::call_once(holder->once, [holder, name, kind] {
    for (const auto &path : candidate_paths("cr.bin"))
      if (holder->data.load(path.string(), kind, name))
        break;
  });
  return holder->data;
}

inline bool legendre(double mu, double nu, double theta, double &p_tilde,
                     double &dtheta, double &d2theta) noexcept {
  if (mu < 0.0 || nu < mu)
    return false;
  std::array<double, 1> value{}, d1{}, d2{};
  if (!atlas("crl", Kind::legendre)
           .evaluate_into(0, {mu, nu - mu, theta}, value.data(), d1.data(),
                          d2.data(), 2))
    return false;
  const double st = std::sin(theta);
  const double ct = std::cos(theta);
  const double scale =
      std::exp(mu * std::log(st) - std::lgamma(mu + 1.0)) *
      std::sqrt((nu + 0.5) * st);
  const double a = (mu + 0.5) * ct / st;
  const double ap = -(mu + 0.5) / (st * st);
  p_tilde = scale * value[0];
  dtheta = scale * (d1[0] + a * value[0]);
  d2theta = scale *
            (d2[0] + 2.0 * a * d1[0] + (a * a + ap) * value[0]);
  return std::isfinite(p_tilde) && std::isfinite(dtheta) &&
         std::isfinite(d2theta);
}

inline bool legendre_degree(int m, double nu, double theta,
                            double &value_nu) noexcept {
  if (m < 0 || nu < m)
    return false;
  std::array<double, 1> value{}, derivative{};
  if (!atlas("crl", Kind::legendre)
           .evaluate_into(0, {static_cast<double>(m), nu - m, theta},
                          value.data(), derivative.data(), nullptr, 1))
    return false;
  const double st = std::sin(theta);
  const double longitude = m == 0 ? 1.0 : 2.0;
  const double factor = std::sqrt(longitude) *
                        std::exp(m * std::log(st) - std::lgamma(m + 1.0));
  value_nu = factor * derivative[0];
  return std::isfinite(value_nu);
}

inline const Atlas &legendre_global_atlas() {
  return atlas("crg", Kind::legendre);
}

inline bool legendre_airy(double inverse_carrier, double alpha,
                          double scaled_zeta, double &real_correction,
                          double &imag_correction) noexcept {
  std::array<double, 2> value{};
  if (!legendre_global_atlas().evaluate_into(
          0, {inverse_carrier, alpha, scaled_zeta}, value.data()))
    return false;
  real_correction = value[0];
  imag_correction = value[1];
  return std::isfinite(real_correction) && std::isfinite(imag_correction);
}

inline bool legendre_outer(double inverse_carrier, double alpha,
                           double log_action, double &real_correction,
                           double &imag_correction) noexcept {
  std::array<double, 2> value{};
  if (!atlas("cro", Kind::legendre)
           .evaluate_into(0, {inverse_carrier, alpha, log_action},
                          value.data()))
    return false;
  real_correction = value[0];
  imag_correction = value[1];
  return std::isfinite(real_correction) && std::isfinite(imag_correction);
}

inline bool legendre_airy_high(double inverse_carrier, double alpha,
                               double scaled_zeta,
                               double &real_correction,
                               double &imag_correction) noexcept {
  std::array<double, 2> value{};
  if (!atlas("crh", Kind::legendre)
           .evaluate_into(0, {inverse_carrier, alpha, scaled_zeta},
                          value.data()))
    return false;
  real_correction = value[0];
  imag_correction = value[1];
  return std::isfinite(real_correction) && std::isfinite(imag_correction);
}

inline bool legendre_outer_high(double inverse_carrier, double alpha,
                                double log_action,
                                double &real_correction,
                                double &imag_correction) noexcept {
  std::array<double, 2> value{};
  if (!atlas("crx", Kind::legendre)
           .evaluate_into(0, {inverse_carrier, alpha, log_action},
                          value.data()))
    return false;
  real_correction = value[0];
  imag_correction = value[1];
  return std::isfinite(real_correction) && std::isfinite(imag_correction);
}

inline bool legendre_outer_patch(double inverse_carrier, double alpha,
                                 double log_action,
                                 double &real_correction,
                                 double &imag_correction) noexcept {
  std::array<double, 2> value{};
  if (!atlas("crp", Kind::legendre)
           .evaluate_into(0, {inverse_carrier, alpha, log_action},
                          value.data()))
    return false;
  real_correction = value[0];
  imag_correction = value[1];
  return std::isfinite(real_correction) && std::isfinite(imag_correction);
}

inline bool legendre_near_evanescent_raw(double inverse_carrier, double gap,
                                         double log_action,
                                         double &log_correction,
                                         double &sampled_interpolation_estimate) noexcept {
  std::array<double, 2> value{0.0,
      std::numeric_limits<double>::infinity()};
  if (!atlas("crn", Kind::legendre)
           .evaluate_into(0, {inverse_carrier, gap, log_action},
                          value.data()))
    return false;
  log_correction = value[0];
  sampled_interpolation_estimate = value[1];
  return std::isfinite(log_correction);
}

inline bool legendre_near_enabled() noexcept {
  return atlas("crn", Kind::legendre).ready();
}

inline bool legendre_near_evanescent(double inverse_carrier, double gap,
                                     double log_action,
                                     double &log_correction) noexcept {
  // Empirical interpolation discrepancy only; not a total error bound.
  double sampled_interpolation_estimate = 0.0;
  return legendre_near_evanescent_raw(inverse_carrier, gap, log_action,
                                      log_correction, sampled_interpolation_estimate) &&
         std::isfinite(sampled_interpolation_estimate) && sampled_interpolation_estimate <= 4.0e-11;
}

inline double mehler_carrier(double order, double lambda, double t) noexcept {
  const double y = std::tanh(t);
  const double root_lambda = std::sqrt(lambda);
  if (order == 0.0)
    return root_lambda * std::asin(y);
  const double total = order * order + lambda;
  const double denominator = std::sqrt(total - lambda * y * y);
  const double argument = order * y / denominator;
  double one_plus_argument = 1.0 + argument;
  if (y < 0.0) {
    const double c = std::cosh(t);
    const double g = 1.0 / (c * c);
    one_plus_argument = total * g /
                        ((denominator - order * y) * denominator);
  }
  const double inverse_hyperbolic =
      0.5 * (std::log(one_plus_argument) - std::log1p(-argument));
  return root_lambda * std::asin(y * std::sqrt(lambda / total)) +
         order * inverse_hyperbolic;
}

inline bool mehler(int m, double tau, double theta, double &log_raw,
                   double &sensitivity, double &dlog_dtheta) noexcept {
  std::array<double, 2> value{}, d1{};
  const double t = std::log(std::tan(0.5 * theta));
  const Atlas &table = atlas("crm", Kind::mehler);
  if (!table.evaluate_into(0, {std::log1p(static_cast<double>(m)),
                               std::log1p(std::fabs(tau)), t},
                           value.data(), d1.data(), nullptr, 2))
    return false;
  const double lambda = tau * tau + 0.25;
  const double potential = static_cast<double>(m) * m +
                           lambda / (std::cosh(t) * std::cosh(t));
  log_raw = mehler_carrier(m, lambda, t) + value[0];
  sensitivity = std::exp(value[1]);
  dlog_dtheta = (std::sqrt(potential) + d1[0]) / std::sin(theta);
  if (m == 0) {
    std::array<double, 2> order_one{}, order_two{};
    const double transformed_tau = std::log1p(std::fabs(tau));
    if (!table.evaluate_into(0, {std::log(2.0), transformed_tau, t},
                             order_one.data(), nullptr, nullptr, 2) ||
        !table.evaluate_into(0, {std::log(3.0), transformed_tau, t},
                             order_two.data(), nullptr, nullptr, 2))
      return false;
    const double log_one = mehler_carrier(1, lambda, t) + order_one[0];
    const double log_two = mehler_carrier(2, lambda, t) + order_two[0] +
                           std::log((lambda + 2.0) / 8.0);
    const double cotangent = -std::sinh(t);
    if (cotangent > 0.0) {
      const double a = log_one + std::log(cotangent);
      const double maximum = std::max(a, log_two);
      log_raw = maximum + std::log(std::exp(a - maximum) +
                                   std::exp(log_two - maximum));
    } else {
      log_raw = log_two;
    }
  }
  return std::isfinite(log_raw) && sensitivity > 0.0 &&
         std::isfinite(dlog_dtheta);
}

inline bool constant_norm(int m, double theta0, double &log_norm) noexcept {
  if (m < 1)
    return false;
  std::array<double, 1> value{};
  const double t = std::log(std::tan(0.5 * theta0));
  if (!atlas("crc", Kind::constant)
           .evaluate_into(0, {std::log1p(static_cast<double>(m)), t, 0.0},
                          value.data()))
    return false;
  log_norm = value[0];
  return std::isfinite(log_norm);
}

} // namespace CRTable
