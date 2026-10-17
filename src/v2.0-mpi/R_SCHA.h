#pragma once
#include "Basis Function.h"
#include "CoordinateTransformation.h"
#include "MPIParallel.h"
#include "RSCHARoots.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
using namespace CoordinateTransformation;

namespace R_SCHA {

// Combined azimuthal terms have finite limits even when m/sin(theta)
// and cot(theta) individually diverge.  Components use the longitude's
// chosen tangent frame at the axis.
struct AxisAngularTerms { double azimuthal, cross, yy; };
inline AxisAngularTerms angular_terms(double m, double theta, double f,
                                      double d1, double d2, double lambda,
                                      double logf, int sign, double ratio) {
  if (theta == 0.0) return {
      m == 1.0 ? d1 : 0.0,
      m == 2.0 ? -d2 : 0.0,
      m == 0.0 ? -d2 : (m == 2.0 ? d2 : 0.0)};
  const double sine = sind(theta), cosine = cosd(theta);
  const auto product = [f, logf, sign](double coefficient) {
    if (coefficient == 0.0 || sign == 0)
      return 0.0;
    if (std::fabs(f) >= std::numeric_limits<double>::min())
      return coefficient * f;
    return sign * std::copysign(static_cast<double>(
        Mehler::ConicalDetail::safe_exp(logf + std::log(std::fabs(coefficient)))),
        coefficient);
  };
  const double azimuthal = product(m / sine);
  if (lambda == 0.0)
    return {azimuthal, -d2, d2};
  const double yy = d2 + product(lambda);
  if (m == 0.0)
    return {0.0, 0.0, yy};
  const double angle = theta * DegreeTrig::radians_per_degree;
  if (angle <= 0.01 &&
      (std::sqrt(std::fabs(lambda)) + m + 1.0) * angle <= 1.5) {
    const double z = std::pow(std::sin(angle / 2.0), 2);
    double term = 1.0, sum = 1.0, derivative = 0.0;
    bool converged = false;
    for (int k = 1; k <= 256; ++k) {
      const double coefficient = ((m+k) * (m+k-1.0) - lambda) /
                                 (double(k) * (m+k));
      const double previous = term * coefficient;
      term = previous * z;
      sum += term;
      const double increment = k * previous;
      derivative += increment;
      if (k >= 2 && std::fabs(increment) <=
          1.0e-18 * std::max(1.0, std::fabs(derivative))) {
        converged = true;
        break;
      }
    }
    if (!converged || sum == 0.0)
      throw std::runtime_error("Polar gradient series did not converge.");
    if (m > 1.0 && sine * sine == 0.0)
      return {azimuthal, -d2, yy};
    const double coefficient = m == 1.0 ? -0.5 * derivative / sum :
        m * ((1.0-m)*cosine/(sine*sine) - 0.5*derivative/sum);
    return {azimuthal, product(coefficient), yy};
  }
  return {azimuthal, product((m / sine) * (cosine / sine - ratio)), yy};
}

inline int checked_index(double value) {
  if (!std::isfinite(value) || value < 0.0 ||
      value > static_cast<double>(std::numeric_limits<int>::max()) ||
      std::floor(value) != value)
    throw std::invalid_argument("Indices and iteration counts must be nonnegative integers.");
  return static_cast<int>(value);
}

inline int coefficient_count(int internal, int external, int mehler, int constant) {
  const long double count = (internal + 1.0L) * (internal + 2.0L) +
      (external + 1.0L) * (external + 2.0L) +
      2.0L * mehler * (mehler + 1.0L) + 2.0L * (constant + 1.0L);
  if (internal < 0 || external < 0 || mehler < 0 || constant < 0 ||
      count > std::numeric_limits<int>::max())
    throw std::invalid_argument("Invalid or excessive model dimensions.");
  return static_cast<int>(count);
}

inline void check_geometry(double longitude, double latitude, double angle,
                           double lower, double upper) {
  if (!std::isfinite(longitude) || !std::isfinite(latitude) ||
      !std::isfinite(angle) || !std::isfinite(lower) || !std::isfinite(upper) ||
      longitude < -180.0 || longitude > 180.0 ||
      latitude < -90.0 || latitude > 90.0 || angle <= 0.0 || angle > 90.0 ||
      lower <= 0.0 || upper <= lower)
    throw std::invalid_argument("Invalid cap coordinates, angle, or radii.");
}


inline void get_nk(int K_Max, double Theta, std::vector<double> &nk) {
  RSCHARoots::compute(K_Max, Theta, nk, MPIParallel::rank(),
                      MPIParallel::size());
  MPIParallel::allreduce_sum(nk);
}

inline void R_SCHA_nk(int K_Max, double Theta, std::vector<double> &nk) {
  if (MPIParallel::rank() == 0) {
    cout << "Function: R_SCHA_nk: " << endl;
    cout << "computing non-integer roots with safeguarded continuation ... "
         << endl << endl;
    cout << "Computing the non-integer roots ... ";
  }

  std::fill(nk.begin(), nk.end(), 0.0);
  get_nk(K_Max, Theta, nk);
  if (std::any_of(nk.begin(), nk.end(),
                  [](double degree) { return !std::isfinite(degree); }))
    throw std::runtime_error("Failed to bracket an R-SCHA degree root");
  if (MPIParallel::rank() == 0)
    cout << "Complete! " << endl;
}

inline void
RSCHA_VectorKernel(double Radius, double Theta, double Theta0, double Lamda,
                   double Min_Radius, double Max_Radius, int Num_CoefIn,
                   int Num_CoefEx, int Num_CoefMe, int Num_CoefCo,
                   std::vector<double> &nk, std::vector<double> &m_Legendre,
                   std::vector<double> &p, std::vector<double> &m_Mehler,
                   std::vector<double> &m_Constant,
                   std::vector<double> &KernelBx, std::vector<double> &KernelBy,
                   std::vector<double> &KernelBz) {

  using namespace Legendre;
  using namespace Mehler;
  using namespace Constant;

  const double pi = (4.0 * std::atan(1.0));
  const double a = Min_Radius;
  const double S = log(Max_Radius / Min_Radius);
  const double rPow1 = std::pow(a / Radius, 1.5);
  const double logrInv = std::log(Radius / a);

  int i;
  int FBP, SBP;

  double RFunc, SinmLamda, CosmLamda, Bznk;
  double Value_Pnm, Value_dPnm;
  double Factor1, Factor2, FN, A, B, RFuncMeh1, RFuncMeh2;
  double Kpm, dKpm;
  double P0m, dP0m;

  for (i = 0; i < Num_CoefIn; ++i) {
    FBP = i;
    SBP = i + Num_CoefIn;

    RFunc = std::pow(a / Radius, nk[i] + 2.0);
    SinmLamda = sind(Lamda * m_Legendre[i]);
    CosmLamda = cosd(Lamda * m_Legendre[i]);
    Bznk = -(nk[i] + 1.0);

    const auto basis = Legendre::EvaluateLegacySchmidt(
        m_Legendre[i], nk[i], Theta * DegreeTrig::radians_per_degree);
    Value_Pnm = basis.status == 0 ? static_cast<double>(basis.value)
        : std::numeric_limits<double>::quiet_NaN();
    Value_dPnm = basis.status == 0 ? static_cast<double>(basis.dtheta)
        : std::numeric_limits<double>::quiet_NaN();
    const auto angular = angular_terms(m_Legendre[i], Theta, Value_Pnm,
        Value_dPnm, static_cast<double>(basis.d2theta), nk[i]*(nk[i]+1.0),
        static_cast<double>(basis.log_abs_value + basis.log_abs_value_lo), basis.sign, static_cast<double>(basis.dtheta_over_value));

    KernelBx[FBP] = RFunc * CosmLamda * Value_dPnm;
    KernelBx[SBP] = RFunc * SinmLamda * Value_dPnm;

    KernelBy[FBP] = RFunc * SinmLamda * angular.azimuthal;
    KernelBy[SBP] = RFunc * CosmLamda * angular.azimuthal * -1.0;

    KernelBz[FBP] = RFunc * CosmLamda * Value_Pnm * Bznk;
    KernelBz[SBP] = RFunc * SinmLamda * Value_Pnm * Bznk;
  }

  for (i = 0; i < Num_CoefEx; ++i) {
    FBP = i + 2 * Num_CoefIn;
    SBP = i + 2 * Num_CoefIn + Num_CoefEx;

    RFunc = std::pow(Radius / a, nk[i] - 1.0);
    SinmLamda = sind(Lamda * m_Legendre[i]);
    CosmLamda = cosd(Lamda * m_Legendre[i]);

    const auto basis = Legendre::EvaluateLegacySchmidt(
        m_Legendre[i], nk[i], Theta * DegreeTrig::radians_per_degree);
    Value_Pnm = basis.status == 0 ? static_cast<double>(basis.value)
        : std::numeric_limits<double>::quiet_NaN();
    Value_dPnm = basis.status == 0 ? static_cast<double>(basis.dtheta)
        : std::numeric_limits<double>::quiet_NaN();
    const auto angular = angular_terms(m_Legendre[i], Theta, Value_Pnm,
        Value_dPnm, static_cast<double>(basis.d2theta), nk[i]*(nk[i]+1.0),
        static_cast<double>(basis.log_abs_value + basis.log_abs_value_lo), basis.sign, static_cast<double>(basis.dtheta_over_value));

    KernelBx[FBP] = RFunc * CosmLamda * Value_dPnm;
    KernelBx[SBP] = RFunc * SinmLamda * Value_dPnm;

    KernelBy[FBP] = RFunc * SinmLamda * angular.azimuthal;
    KernelBy[SBP] = RFunc * CosmLamda * angular.azimuthal * -1.0;

    KernelBz[FBP] = RFunc * CosmLamda * Value_Pnm * nk[i];
    KernelBz[SBP] = RFunc * SinmLamda * Value_Pnm * nk[i];
  }

  for (i = 0; i < Num_CoefMe; ++i) {
    FBP = i + 2 * (Num_CoefIn + Num_CoefEx);
    SBP = i + 2 * (Num_CoefIn + Num_CoefEx) + Num_CoefMe;

    Factor1 = std::sqrt(2.0 / (S * (std::pow(2.0 * pi * p[i], 2.0) + S * S)));
    Factor2 = (pi / S) * p[i];
    FN = p[i] * pi / S;

    A = (2.0 * pi * p[i]) * std::cos((pi / S) * p[i] * logrInv) +
        S * std::sin((pi / S) * p[i] * logrInv);
    B = -(2.0 * pi * p[i]) * std::sin((pi / S) * p[i] * logrInv) +
        S * std::cos((pi / S) * p[i] * logrInv);

    RFuncMeh1 = Factor1 * A * rPow1;
    RFuncMeh2 = (-0.5 * Factor1 * A * rPow1) + (Factor1 * Factor2 * B * rPow1);

    SinmLamda = sind(Lamda * m_Mehler[i]);
    CosmLamda = cosd(Lamda * m_Mehler[i]);

    const auto basis = Mehler::Evaluate(
        m_Mehler[i], FN, Theta, Theta0);
    Kpm = basis.status == 0 ? basis.value : std::numeric_limits<double>::quiet_NaN();
    dKpm = basis.status == 0 ? basis.dtheta : std::numeric_limits<double>::quiet_NaN();
    const auto angular = angular_terms(m_Mehler[i], Theta, Kpm, dKpm,
        basis.d2theta, -(FN*FN+0.25), static_cast<double>(basis.log_abs_value), basis.sign, static_cast<double>(basis.dtheta_over_value));

    KernelBx[FBP] = RFuncMeh1 * CosmLamda * dKpm;
    KernelBx[SBP] = RFuncMeh1 * SinmLamda * dKpm;

    KernelBy[FBP] = RFuncMeh1 * SinmLamda * angular.azimuthal;
    KernelBy[SBP] = RFuncMeh1 * CosmLamda * angular.azimuthal * -1.0;

    KernelBz[FBP] = RFuncMeh2 * CosmLamda * Kpm;
    KernelBz[SBP] = RFuncMeh2 * SinmLamda * Kpm;
  }

  for (i = 0; i < Num_CoefCo; ++i) {
    FBP = i + 2 * (Num_CoefIn + Num_CoefEx + Num_CoefMe);
    SBP = i + 2 * (Num_CoefIn + Num_CoefEx + Num_CoefMe) + Num_CoefCo;

    RFunc = (1.0 / std::sqrt(std::exp(S) - 1.0)) * (a / Radius);
    SinmLamda = sind(Lamda * m_Constant[i]);
    CosmLamda = cosd(Lamda * m_Constant[i]);

    const auto basis = Constant::Evaluate(
        m_Constant[i], Theta, Theta0);
    P0m = basis.status == 0 ? basis.value : std::numeric_limits<double>::quiet_NaN();
    dP0m = basis.status == 0 ? basis.dtheta : std::numeric_limits<double>::quiet_NaN();
    const auto angular = angular_terms(m_Constant[i], Theta, P0m, dP0m,
        basis.d2theta, 0.0, static_cast<double>(basis.log_abs_value), basis.sign, static_cast<double>(basis.dtheta_over_value));

    KernelBx[FBP] = RFunc * CosmLamda * dP0m;
    KernelBx[SBP] = RFunc * SinmLamda * dP0m;

    KernelBy[FBP] = RFunc * SinmLamda * angular.azimuthal;
    KernelBy[SBP] = RFunc * CosmLamda * angular.azimuthal * -1.0;

    KernelBz[FBP] = 0.0;
    KernelBz[SBP] = 0.0;
  }
}

inline void RSCHA_TensorKernel(
    double Radius, double Theta, double Theta0, double Lamda,
    double Min_Radius, double Max_Radius, int Num_CoefIn, int Num_CoefEx,
    int Num_CoefMe, int Num_CoefCo, std::vector<double> &nk,
    std::vector<double> &m_Legendre, std::vector<double> &p,
    std::vector<double> &m_Mehler, std::vector<double> &m_Constant,
    std::vector<double> &KernelV, std::vector<double> &KernelBx,
    std::vector<double> &KernelBy, std::vector<double> &KernelBz,
    std::vector<double> &KernelTxx, std::vector<double> &KernelTxy,
    std::vector<double> &KernelTyy, std::vector<double> &KernelTzx,
    std::vector<double> &KernelTzy, std::vector<double> &KernelTzz) {

  using namespace Legendre;
  using namespace Mehler;
  using namespace Constant;

  const double pi = (4.0 * std::atan(1.0));

  const double a = Min_Radius;
  const double S = log(Max_Radius / Min_Radius);
  const double rInv = 1.0 / Radius;
  const double rPow1 = std::pow(a / Radius, 1.5);
  const double rPow2 = std::pow(a / Radius, 2.5);
  const double logrInv = std::log(Radius / a);
  const double rInvC = 1.0 / Radius;

  int i;
  int FBP, SBP;

  double RFunc, RFuncT, SinmLamda, CosmLamda, Bznk;
  double Value_Pnm, Value_dPnm, Value_ddPnm;
  double Factor1, Factor2, FN, A, B, RFuncMeh1, RFuncMeh2, RFuncMeh3;
  double Kpm, dKpm, ddKpm;
  double P0m, dP0m, ddP0m;

  for (i = 0; i < Num_CoefIn; ++i) {
    FBP = i;
    SBP = i + Num_CoefIn;

    RFunc = std::pow(a / Radius, nk[i] + 2.0);
    RFuncT = (1.0 / Radius) * RFunc;
    SinmLamda = sind(Lamda * m_Legendre[i]);
    CosmLamda = cosd(Lamda * m_Legendre[i]);
    Bznk = -(nk[i] + 1.0);

    const auto basis = Legendre::EvaluateLegacySchmidt(
        m_Legendre[i], nk[i], Theta * DegreeTrig::radians_per_degree);
    Value_Pnm = basis.status == 0 ? static_cast<double>(basis.value)
        : std::numeric_limits<double>::quiet_NaN();
    Value_dPnm = basis.status == 0 ? static_cast<double>(basis.dtheta)
        : std::numeric_limits<double>::quiet_NaN();
    Value_ddPnm = basis.status == 0 ? static_cast<double>(basis.d2theta)
        : std::numeric_limits<double>::quiet_NaN();
    const auto angular = angular_terms(m_Legendre[i], Theta, Value_Pnm,
        Value_dPnm, static_cast<double>(basis.d2theta), nk[i]*(nk[i]+1.0),
        static_cast<double>(basis.log_abs_value + basis.log_abs_value_lo), basis.sign, static_cast<double>(basis.dtheta_over_value));

    KernelV[FBP] = RFunc * Radius * CosmLamda * Value_Pnm;
    KernelV[SBP] = RFunc * Radius * SinmLamda * Value_Pnm;

    KernelBx[FBP] = RFunc * CosmLamda * Value_dPnm;
    KernelBx[SBP] = RFunc * SinmLamda * Value_dPnm;

    KernelBy[FBP] = RFunc * SinmLamda * angular.azimuthal;
    KernelBy[SBP] = RFunc * CosmLamda * angular.azimuthal * -1.0;

    KernelBz[FBP] = RFunc * CosmLamda * Value_Pnm * Bznk;
    KernelBz[SBP] = RFunc * SinmLamda * Value_Pnm * Bznk;

    KernelTxx[FBP] =
        RFuncT * CosmLamda * (-1.0 * Value_ddPnm - Bznk * Value_Pnm);
    KernelTxx[SBP] =
        RFuncT * SinmLamda * (-1.0 * Value_ddPnm - Bznk * Value_Pnm);

    KernelTxy[FBP] =
        RFuncT * SinmLamda * angular.cross;
    KernelTxy[SBP] = RFuncT * CosmLamda * angular.cross * -1.0;

    KernelTyy[FBP] =
        RFuncT * CosmLamda *
        (angular.yy - Bznk * Value_Pnm);
    KernelTyy[SBP] =
        RFuncT * SinmLamda *
        (angular.yy - Bznk * Value_Pnm);

    KernelTzx[FBP] = RFuncT * CosmLamda * ((-Bznk + 1.0) * Value_dPnm);
    KernelTzx[SBP] = RFuncT * SinmLamda * ((-Bznk + 1.0) * Value_dPnm);

    KernelTzy[FBP] =
        RFuncT * SinmLamda * (angular.azimuthal * (-Bznk + 1.0));
    KernelTzy[SBP] =
        RFuncT * CosmLamda * (angular.azimuthal * (-Bznk + 1.0)) * -1.0;

    KernelTzz[FBP] = RFuncT * CosmLamda * (Bznk * (-Bznk + 1.0)) * Value_Pnm;
    KernelTzz[SBP] = RFuncT * SinmLamda * (Bznk * (-Bznk + 1.0)) * Value_Pnm;
  }

  for (i = 0; i < Num_CoefEx; ++i) {
    FBP = i + 2 * Num_CoefIn;
    SBP = i + 2 * Num_CoefIn + Num_CoefEx;

    RFunc = std::pow(Radius / a, nk[i] - 1.0);
    RFuncT = (1.0 / Radius) * RFunc;
    SinmLamda = sind(Lamda * m_Legendre[i]);
    CosmLamda = cosd(Lamda * m_Legendre[i]);

    const auto basis = Legendre::EvaluateLegacySchmidt(
        m_Legendre[i], nk[i], Theta * DegreeTrig::radians_per_degree);
    Value_Pnm = basis.status == 0 ? static_cast<double>(basis.value)
        : std::numeric_limits<double>::quiet_NaN();
    Value_dPnm = basis.status == 0 ? static_cast<double>(basis.dtheta)
        : std::numeric_limits<double>::quiet_NaN();
    Value_ddPnm = basis.status == 0 ? static_cast<double>(basis.d2theta)
        : std::numeric_limits<double>::quiet_NaN();
    const auto angular = angular_terms(m_Legendre[i], Theta, Value_Pnm,
        Value_dPnm, static_cast<double>(basis.d2theta), nk[i]*(nk[i]+1.0),
        static_cast<double>(basis.log_abs_value + basis.log_abs_value_lo), basis.sign, static_cast<double>(basis.dtheta_over_value));

    KernelV[FBP] = RFunc * Radius * CosmLamda * Value_Pnm;
    KernelV[SBP] = RFunc * Radius * SinmLamda * Value_Pnm;

    KernelBx[FBP] = RFunc * CosmLamda * Value_dPnm;
    KernelBx[SBP] = RFunc * SinmLamda * Value_dPnm;

    KernelBy[FBP] = RFunc * SinmLamda * angular.azimuthal;
    KernelBy[SBP] = RFunc * CosmLamda * angular.azimuthal * -1.0;

    KernelBz[FBP] = RFunc * CosmLamda * Value_Pnm * nk[i];
    KernelBz[SBP] = RFunc * SinmLamda * Value_Pnm * nk[i];

    KernelTxx[FBP] =
        RFuncT * CosmLamda * (-1.0 * Value_ddPnm - nk[i] * Value_Pnm);
    KernelTxx[SBP] =
        RFuncT * SinmLamda * (-1.0 * Value_ddPnm - nk[i] * Value_Pnm);

    KernelTxy[FBP] =
        RFuncT * SinmLamda * angular.cross;
    KernelTxy[SBP] = RFuncT * CosmLamda * angular.cross * -1.0;

    KernelTyy[FBP] =
        RFuncT * CosmLamda *
        (angular.yy - nk[i] * Value_Pnm);
    KernelTyy[SBP] =
        RFuncT * SinmLamda *
        (angular.yy - nk[i] * Value_Pnm);

    KernelTzx[FBP] = RFuncT * CosmLamda * (-(nk[i] - 1.0) * Value_dPnm);
    KernelTzx[SBP] = RFuncT * SinmLamda * (-(nk[i] - 1.0) * Value_dPnm);

    KernelTzy[FBP] =
        RFuncT * SinmLamda * (angular.azimuthal * (1.0 - nk[i]));
    KernelTzy[SBP] =
        RFuncT * CosmLamda * (angular.azimuthal * (1.0 - nk[i])) * -1.0;

    KernelTzz[FBP] = RFuncT * CosmLamda * (-nk[i] * (nk[i] - 1.0)) * Value_Pnm;
    KernelTzz[SBP] = RFuncT * SinmLamda * (-nk[i] * (nk[i] - 1.0)) * Value_Pnm;
  }

  for (i = 0; i < Num_CoefMe; ++i) {
    FBP = i + 2 * (Num_CoefIn + Num_CoefEx);
    SBP = i + 2 * (Num_CoefIn + Num_CoefEx) + Num_CoefMe;

    Factor1 = std::sqrt(2.0 / (S * (std::pow(2.0 * pi * p[i], 2.0) + S * S)));
    Factor2 = (pi / S) * p[i];
    FN = p[i] * pi / S;

    A = (2.0 * pi * p[i]) * std::cos((pi / S) * p[i] * logrInv) +
        S * std::sin((pi / S) * p[i] * logrInv);
    B = -(2.0 * pi * p[i]) * std::sin((pi / S) * p[i] * logrInv) +
        S * std::cos((pi / S) * p[i] * logrInv);

    RFuncMeh1 = Factor1 * A * rPow1;

    RFuncMeh2 = -0.5 * Factor1 * A * rPow1 + Factor1 * Factor2 * B * rPow1;

    RFuncMeh3 = (1.0 / a) * (0.75 - Factor2 * Factor2) * Factor1 * rPow2 * A -
                2.0 * (1.0 / a) * Factor1 * Factor2 * rPow2 * B;

    SinmLamda = sind(Lamda * m_Mehler[i]);
    CosmLamda = cosd(Lamda * m_Mehler[i]);

    const auto basis = Mehler::Evaluate(
        m_Mehler[i], FN, Theta, Theta0);
    Kpm = basis.status == 0 ? basis.value : std::numeric_limits<double>::quiet_NaN();
    dKpm = basis.status == 0 ? basis.dtheta : std::numeric_limits<double>::quiet_NaN();
    ddKpm = basis.status == 0 ? basis.d2theta : std::numeric_limits<double>::quiet_NaN();
    const auto angular = angular_terms(m_Mehler[i], Theta, Kpm, dKpm,
        basis.d2theta, -(FN*FN+0.25), static_cast<double>(basis.log_abs_value), basis.sign, static_cast<double>(basis.dtheta_over_value));

    KernelV[FBP] = RFuncMeh1 * Radius * CosmLamda * Kpm;
    KernelV[SBP] = RFuncMeh1 * Radius * SinmLamda * Kpm;

    KernelBx[FBP] = RFuncMeh1 * CosmLamda * dKpm;
    KernelBx[SBP] = RFuncMeh1 * SinmLamda * dKpm;

    KernelBy[FBP] = RFuncMeh1 * SinmLamda * angular.azimuthal;
    KernelBy[SBP] = RFuncMeh1 * CosmLamda * angular.azimuthal * -1.0;

    KernelBz[FBP] = RFuncMeh2 * CosmLamda * Kpm;
    KernelBz[SBP] = RFuncMeh2 * SinmLamda * Kpm;

    KernelTxx[FBP] =
        CosmLamda * rInv * (RFuncMeh1 * ddKpm + RFuncMeh2 * Kpm) * -1.0;
    KernelTxx[SBP] =
        SinmLamda * rInv * (RFuncMeh1 * ddKpm + RFuncMeh2 * Kpm) * -1.0;

    KernelTxy[FBP] =
        RFuncMeh1 * SinmLamda * rInv * angular.cross;
    KernelTxy[SBP] = RFuncMeh1 * CosmLamda * rInv * angular.cross * -1.0;

    KernelTyy[FBP] = RFuncMeh1 * CosmLamda * rInv *
                         angular.yy -
                     RFuncMeh2 * CosmLamda * rInv * Kpm;
    KernelTyy[SBP] = RFuncMeh1 * SinmLamda * rInv *
                         angular.yy -
                     RFuncMeh2 * SinmLamda * rInv * Kpm;

    KernelTzx[FBP] = CosmLamda * rInv * dKpm * (RFuncMeh1 - RFuncMeh2);
    KernelTzx[SBP] = SinmLamda * rInv * dKpm * (RFuncMeh1 - RFuncMeh2);

    KernelTzy[FBP] =
        SinmLamda * rInv * angular.azimuthal * (RFuncMeh1 - RFuncMeh2);
    KernelTzy[SBP] =
        CosmLamda * rInv * angular.azimuthal * (RFuncMeh1 - RFuncMeh2) * -1.0;

    KernelTzz[FBP] = -RFuncMeh3 * CosmLamda * Kpm;
    KernelTzz[SBP] = -RFuncMeh3 * SinmLamda * Kpm;
  }

  for (i = 0; i < Num_CoefCo; ++i) {
    FBP = i + 2 * (Num_CoefIn + Num_CoefEx + Num_CoefMe);
    SBP = i + 2 * (Num_CoefIn + Num_CoefEx + Num_CoefMe) + Num_CoefCo;

    RFunc = (1.0 / std::sqrt(std::exp(S) - 1.0)) * (a / Radius);
    SinmLamda = sind(Lamda * m_Constant[i]);
    CosmLamda = cosd(Lamda * m_Constant[i]);

    const auto basis = Constant::Evaluate(
        m_Constant[i], Theta, Theta0);
    P0m = basis.status == 0 ? basis.value : std::numeric_limits<double>::quiet_NaN();
    dP0m = basis.status == 0 ? basis.dtheta : std::numeric_limits<double>::quiet_NaN();
    ddP0m = basis.status == 0 ? basis.d2theta : std::numeric_limits<double>::quiet_NaN();
    const auto angular = angular_terms(m_Constant[i], Theta, P0m, dP0m,
        basis.d2theta, 0.0, static_cast<double>(basis.log_abs_value), basis.sign, static_cast<double>(basis.dtheta_over_value));

    KernelV[FBP] = RFunc * Radius * CosmLamda * P0m;
    KernelV[SBP] = RFunc * Radius * SinmLamda * P0m;

    KernelBx[FBP] = RFunc * CosmLamda * dP0m;
    KernelBx[SBP] = RFunc * SinmLamda * dP0m;

    KernelBy[FBP] = RFunc * SinmLamda * angular.azimuthal;
    KernelBy[SBP] = RFunc * CosmLamda * angular.azimuthal * -1.0;

    KernelBz[FBP] = 0.0;
    KernelBz[SBP] = 0.0;

    KernelTxx[FBP] = RFunc * CosmLamda * rInvC * ddP0m * -1.0;
    KernelTxx[SBP] = RFunc * SinmLamda * rInvC * ddP0m * -1.0;

    KernelTxy[FBP] =
        RFunc * SinmLamda * rInvC * angular.cross;
    KernelTxy[SBP] =
        RFunc * CosmLamda * rInvC * angular.cross * -1.0;

    KernelTyy[FBP] = RFunc * CosmLamda * rInvC *
                     angular.yy;
    KernelTyy[SBP] = RFunc * SinmLamda * rInvC *
                     angular.yy;

    KernelTzx[FBP] = RFunc * CosmLamda * rInvC * dP0m;
    KernelTzx[SBP] = RFunc * SinmLamda * rInvC * dP0m;

    KernelTzy[FBP] = RFunc * SinmLamda * rInvC * angular.azimuthal;
    KernelTzy[SBP] = RFunc * CosmLamda * rInvC * angular.azimuthal * -1.0;

    KernelTzz[FBP] = 0.0;
    KernelTzz[SBP] = 0.0;
  }
}

inline int readmatrix(const std::string &filename,
                      std::vector<std::vector<double>> &Matrix) {
  std::ifstream fp_File(filename);
  if (!fp_File) {
    std::cout << "Fail to open file: " << filename << "\n";
    return -1;
  }

  auto normalize = [](std::string &s) {
    for (char &ch : s) {
      if (ch == ',' || ch == ';' || ch == '\t' || ch == '\r')
        ch = ' ';
    }
  };

  auto parse_row = [&](const std::string &line_in,
                       std::vector<double> &row_out) -> void {
    row_out.clear();
    std::string line = line_in;
    normalize(line);
    std::istringstream iss(line);
    std::string token;
    while (iss >> token) {
      std::size_t consumed = 0;
      const double v = std::stod(token, &consumed);
      if (consumed != token.size() || !std::isfinite(v))
        throw std::invalid_argument("Non-finite value in input file: " + filename);
      row_out.push_back(v);
    }
    if (!iss.eof())
      throw std::invalid_argument("Invalid numeric input in file: " + filename);
  };

  Matrix.clear();

  std::string line;
  std::vector<double> row;
  size_t ncol = 0;

  while (std::getline(fp_File, line)) {
    parse_row(line, row);
    if (!row.empty()) {
      ncol = row.size();
      Matrix.push_back(row);
      break;
    }
  }

  if (ncol == 0) {

    return -1;
  }

  while (std::getline(fp_File, line)) {
    parse_row(line, row);
    if (row.empty())
      continue;

    if (row.size() != ncol) {
      std::cout << "Error: inconsistent column count in file: " << filename
                << "\n";
      std::cout << "Expected " << ncol << " columns, but got " << row.size()
                << ".\n";
      return -1;
    }

    Matrix.push_back(row);
  }

  return fp_File.bad() || Matrix.empty() ? -1 : 1;
}

inline int writematrix(const string &filename,
                       const vector<vector<double>> &Matrix) {
  ofstream fp(filename);
  if (!fp) {
    cout << "Fail to open file: " << filename << "\n";
    return -1;
  }
  fp << fixed << setprecision(15);
  for (const auto &row : Matrix) {
    for (size_t j = 0; j < row.size(); ++j) {
      if (j)
        fp << ',';
      fp << row[j];
    }
    fp << '\n';
  }
  fp.close();
  return fp ? 1 : -1;
}

inline int ReadModellingParameters(

    const string &ModellingParametersFile,
    const string &InversionParametersFile, const string &InitialModelFile,
    const string &ModelWeightFile, const string &PriorModelFile,
    const bool IncludeTMI, const bool IncludeVector, const bool IncludeDiffTMI,
    const bool IncludeDiffVector, const string &DataFile_TMI,
    const string &DataFile_Vector, const string &DataFile_DiffTMI,
    const string &DataFile_DiffVector, const string &DataWeightFile_TMI,
    const string &DataWeightFile_Vector, const string &DataWeightFile_DiffTMI,
    const string &DataWeightFile_DiffVector, const string &SigmaFile_TMI,
    const string &SigmaFile_Vector, const string &SigmaFile_DiffTMI,
    const string &SigmaFile_DiffVector,

    double &Lon0, double &Lat0, double &Theta0, int &K_Max_Internal,
    int &K_Max_External, int &P_Max, int &M_Max, double &Max_Radius,
    double &Min_Radius,

    double &epsilon, double &damping, int &Max_Iteration, int &RobustSwitch,
    double &RobustConstant,

    std::vector<double> &LongitudeTMI, std::vector<double> &LatitudeTMI,
    std::vector<double> &RadiusTMI, std::vector<double> &IncTMI,
    std::vector<double> &DecTMI, std::vector<double> &TMI,
    std::vector<double> &DataWeight_TMI, std::vector<double> &Sigma_TMI,

    std::vector<double> &LongitudeVector, std::vector<double> &LatitudeVector,
    std::vector<double> &RadiusVector, std::vector<double> &Bx,
    std::vector<double> &By, std::vector<double> &Bz,
    std::vector<double> &DataWeight_Bx, std::vector<double> &DataWeight_By,
    std::vector<double> &DataWeight_Bz, std::vector<double> &Sigma_Bx,
    std::vector<double> &Sigma_By, std::vector<double> &Sigma_Bz,

    std::vector<double> &LongitudeDiffTMI1,
    std::vector<double> &LatitudeDiffTMI1, std::vector<double> &RadiusDiffTMI1,
    std::vector<double> &IncDiffTMI1, std::vector<double> &DecDiffTMI1,
    std::vector<double> &DiffTMI1, std::vector<double> &LongitudeDiffTMI2,
    std::vector<double> &LatitudeDiffTMI2, std::vector<double> &RadiusDiffTMI2,
    std::vector<double> &IncDiffTMI2, std::vector<double> &DecDiffTMI2,
    std::vector<double> &DiffTMI2, std::vector<double> &DataWeight_DiffTMI,
    std::vector<double> &Sigma_DiffTMI,

    std::vector<double> &LongitudeDiffVector1,
    std::vector<double> &LatitudeDiffVector1,
    std::vector<double> &RadiusDiffVector1, std::vector<double> &DiffBx1,
    std::vector<double> &DiffBy1, std::vector<double> &DiffBz1,
    std::vector<double> &LongitudeDiffVector2,
    std::vector<double> &LatitudeDiffVector2,
    std::vector<double> &RadiusDiffVector2, std::vector<double> &DiffBx2,
    std::vector<double> &DiffBy2, std::vector<double> &DiffBz2,
    std::vector<double> &DataWeight_DiffBx,
    std::vector<double> &DataWeight_DiffBy,
    std::vector<double> &DataWeight_DiffBz, std::vector<double> &Sigma_DiffBx,
    std::vector<double> &Sigma_DiffBy, std::vector<double> &Sigma_DiffBz,
    std::vector<double> &m0, std::vector<double> &Weight_Model,
    std::vector<double> &Prior_Model) {

  cout << "\nReading the modelling parameters and data ..." << endl << endl;

  vector<vector<double>> MatrixTemp;
  int NumRows, NumCols;
  int Flag = 0;

  {

    Flag = readmatrix(ModellingParametersFile, MatrixTemp);

    if (Flag != 1) {
      cout << "An error occurred while reading the modelling parameters file."
           << endl;
      return -1;
    }
    NumRows = MatrixTemp.size();
    NumCols = MatrixTemp[0].size();

    if ((NumRows != 9) || (NumCols != 1)) {
      cout << "An error occurred while reading the modelling parameters file."
           << endl;
      cout << "Please check the matrix dimensions in the input text file. "
              "(Expected: 9 * 1)"
           << endl;
      cout << ModellingParametersFile << endl;
      return -1;
    }

    Lon0 = MatrixTemp[0][0];
    Lat0 = MatrixTemp[1][0];
    Theta0 = MatrixTemp[2][0];
    K_Max_Internal = checked_index(MatrixTemp[3][0]);
    K_Max_External = checked_index(MatrixTemp[4][0]);
    P_Max = checked_index(MatrixTemp[5][0]);
    M_Max = checked_index(MatrixTemp[6][0]);
    Max_Radius = MatrixTemp[7][0];
    Min_Radius = MatrixTemp[8][0];
    check_geometry(Lon0, Lat0, Theta0, Min_Radius, Max_Radius);
    coefficient_count(K_Max_Internal, K_Max_External, P_Max, M_Max);

    MatrixTemp.resize(0);

    cout << "Modelling parameters loaded from:" << "\n     "
         << ModellingParametersFile << endl;
  }

  {

    Flag = readmatrix(InversionParametersFile, MatrixTemp);

    if (Flag != 1) {
      cout << "An error occurred while reading the inversion parameters file."
           << endl;
      return -2;
    }

    NumRows = MatrixTemp.size();
    NumCols = MatrixTemp[0].size();

    if ((NumRows != 5) || (NumCols != 1)) {
      cout << "An error occurred while reading the inversion parameters file."
           << endl;
      cout << "Please check the matrix dimensions in the input text file. "
              "(Expected: 5 * 1)"
           << endl;
      cout << InversionParametersFile << endl;
      return -2;
    }

    epsilon = MatrixTemp[0][0];
    damping = MatrixTemp[1][0];
    Max_Iteration = checked_index(MatrixTemp[2][0]);
    RobustSwitch = checked_index(MatrixTemp[3][0]);
    RobustConstant = MatrixTemp[4][0];

    MatrixTemp.resize(0);

    cout << "Inversion parameters loaded from:" << "\n     "
         << InversionParametersFile << endl;
  }

  {
    if (IncludeTMI) {

      Flag = readmatrix(DataFile_TMI, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data file of TMI data."
             << endl;
        return -3;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if (NumCols != 6) {
        cout << "An error occurred while reading the data file of TMI data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file. "
                "(Expected: 6 columns)"
             << endl;
        cout << DataFile_TMI << endl;
        return -3;
      }

      LongitudeTMI.resize(NumRows);
      LatitudeTMI.resize(NumRows);
      RadiusTMI.resize(NumRows);
      TMI.resize(NumRows);
      IncTMI.resize(NumRows);
      DecTMI.resize(NumRows);

      for (int i = 0; i < NumRows; ++i) {
        LongitudeTMI[i] = MatrixTemp[i][0];
        LatitudeTMI[i] = MatrixTemp[i][1];
        RadiusTMI[i] = MatrixTemp[i][2];
        TMI[i] = MatrixTemp[i][3];
        IncTMI[i] = MatrixTemp[i][4];
        DecTMI[i] = MatrixTemp[i][5];
      }

      MatrixTemp.resize(0);

      Flag = readmatrix(DataWeightFile_TMI, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data weight file of TMI "
                "data."
             << endl;
        return -4;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if ((static_cast<std::size_t>(NumRows) != LongitudeTMI.size()) || (NumCols != 1)) {
        cout << "An error occurred while reading the data weight file of TMI "
                "data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file."
             << endl;
        cout << DataWeightFile_TMI << endl;
        return -4;
      }

      DataWeight_TMI.resize(NumRows);
      for (int i = 0; i < NumRows; ++i) {
        DataWeight_TMI[i] = MatrixTemp[i][0];
      }

      MatrixTemp.resize(0);

      Flag = readmatrix(SigmaFile_TMI, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data sigma file of TMI "
                "data."
             << endl;
        return -4;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if ((static_cast<std::size_t>(NumRows) != LongitudeTMI.size()) || (NumCols != 1)) {
        cout << "An error occurred while reading the data sigma file of TMI "
                "data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file."
             << endl;
        cout << SigmaFile_TMI << endl;
        return -4;
      }

      Sigma_TMI.resize(NumRows);
      for (int i = 0; i < NumRows; ++i) {
        Sigma_TMI[i] = MatrixTemp[i][0];
      }

      MatrixTemp.resize(0);

      cout << NumRows << " TMI data points loaded from: " << "\n     "
           << DataFile_TMI << endl;
      cout << "Corresponding data weights have been loaded from: " << "\n     "
           << DataWeightFile_TMI << endl;
      cout << "Corresponding data sigma have been loaded from: " << "\n     "
           << SigmaFile_TMI << endl;
    }
  }

  {
    if (IncludeVector) {

      Flag = readmatrix(DataFile_Vector, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data file of Vector data."
             << endl;
        return -5;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if (NumCols != 6) {
        cout << "An error occurred while reading the data file of Vector data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file. "
                "(Expected: 6 columns)"
             << endl;
        cout << DataFile_Vector << endl;
        return -5;
      }

      LongitudeVector.resize(NumRows);
      LatitudeVector.resize(NumRows);
      RadiusVector.resize(NumRows);
      Bx.resize(NumRows);
      By.resize(NumRows);
      Bz.resize(NumRows);

      for (int i = 0; i < NumRows; ++i) {
        LongitudeVector[i] = MatrixTemp[i][0];
        LatitudeVector[i] = MatrixTemp[i][1];
        RadiusVector[i] = MatrixTemp[i][2];
        Bx[i] = MatrixTemp[i][3];
        By[i] = MatrixTemp[i][4];
        Bz[i] = MatrixTemp[i][5];
      }

      MatrixTemp.resize(0);

      Flag = readmatrix(DataWeightFile_Vector, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data weight file of "
                "Vector data."
             << endl;
        return -6;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if ((static_cast<std::size_t>(NumRows) != LongitudeVector.size()) || (NumCols != 3)) {
        cout << "An error occurred while reading the data weight file of "
                "Vector data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file."
             << endl;
        cout << DataWeightFile_Vector << endl;
        return -6;
      }

      DataWeight_Bx.resize(NumRows);
      DataWeight_By.resize(NumRows);
      DataWeight_Bz.resize(NumRows);
      for (int i = 0; i < NumRows; ++i) {
        DataWeight_Bx[i] = MatrixTemp[i][0];
        DataWeight_By[i] = MatrixTemp[i][1];
        DataWeight_Bz[i] = MatrixTemp[i][2];
      }

      MatrixTemp.resize(0);

      Flag = readmatrix(SigmaFile_Vector, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data sigma file of Vector "
                "data."
             << endl;
        return -6;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if ((static_cast<std::size_t>(NumRows) != LongitudeVector.size()) || (NumCols != 3)) {
        cout << "An error occurred while reading the data sigma file of Vector "
                "data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file."
             << endl;
        cout << SigmaFile_Vector << endl;
        return -6;
      }

      Sigma_Bx.resize(NumRows);
      Sigma_By.resize(NumRows);
      Sigma_Bz.resize(NumRows);
      for (int i = 0; i < NumRows; ++i) {
        Sigma_Bx[i] = MatrixTemp[i][0];
        Sigma_By[i] = MatrixTemp[i][1];
        Sigma_Bz[i] = MatrixTemp[i][2];
      }

      MatrixTemp.resize(0);

      cout << NumRows << " vector data points loaded from: " << "\n     "
           << DataFile_Vector << endl;
      cout << "Corresponding data weights have been loaded from: " << "\n     "
           << DataWeightFile_Vector << endl;
      cout << "Corresponding data sigma have been loaded from: " << "\n     "
           << SigmaFile_Vector << endl;
    }
  }

  {
    if (IncludeDiffTMI) {

      Flag = readmatrix(DataFile_DiffTMI, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data file of along track "
                "differenced TMI data."
             << endl;
        return -7;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if (NumCols != 12) {
        cout << "An error occurred while reading the data file of along track "
                "differenced TMI data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file. "
                "(Expected: 12 columns)"
             << endl;
        cout << DataFile_DiffTMI << endl;
        return -7;
      }

      LongitudeDiffTMI1.resize(NumRows);
      LatitudeDiffTMI1.resize(NumRows);
      RadiusDiffTMI1.resize(NumRows);
      DiffTMI1.resize(NumRows);
      IncDiffTMI1.resize(NumRows);
      DecDiffTMI1.resize(NumRows);
      LongitudeDiffTMI2.resize(NumRows);
      LatitudeDiffTMI2.resize(NumRows);
      RadiusDiffTMI2.resize(NumRows);
      DiffTMI2.resize(NumRows);
      IncDiffTMI2.resize(NumRows);
      DecDiffTMI2.resize(NumRows);

      for (int i = 0; i < NumRows; ++i) {
        LongitudeDiffTMI1[i] = MatrixTemp[i][0];
        LatitudeDiffTMI1[i] = MatrixTemp[i][1];
        RadiusDiffTMI1[i] = MatrixTemp[i][2];
        DiffTMI1[i] = MatrixTemp[i][3];
        IncDiffTMI1[i] = MatrixTemp[i][4];
        DecDiffTMI1[i] = MatrixTemp[i][5];
        LongitudeDiffTMI2[i] = MatrixTemp[i][6];
        LatitudeDiffTMI2[i] = MatrixTemp[i][7];
        RadiusDiffTMI2[i] = MatrixTemp[i][8];
        DiffTMI2[i] = MatrixTemp[i][9];
        IncDiffTMI2[i] = MatrixTemp[i][10];
        DecDiffTMI2[i] = MatrixTemp[i][11];
      }

      MatrixTemp.resize(0);

      Flag = readmatrix(DataWeightFile_DiffTMI, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data weight file of along "
                "track differenced TMI data."
             << endl;
        return -8;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if ((static_cast<std::size_t>(NumRows) != LongitudeDiffTMI1.size()) || (NumCols != 1)) {
        cout << "An error occurred while reading the data weight file of along "
                "track differenced TMI data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file."
             << endl;
        cout << DataWeightFile_DiffTMI << endl;
        return -8;
      }

      DataWeight_DiffTMI.resize(NumRows);
      for (int i = 0; i < NumRows; ++i) {
        DataWeight_DiffTMI[i] = MatrixTemp[i][0];
      }

      MatrixTemp.resize(0);

      Flag = readmatrix(SigmaFile_DiffTMI, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data sigma file of along "
                "track differenced TMI data."
             << endl;
        return -8;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if ((static_cast<std::size_t>(NumRows) != LongitudeDiffTMI1.size()) || (NumCols != 1)) {
        cout << "An error occurred while reading the data sigma file of along "
                "track differenced TMI data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file."
             << endl;
        cout << SigmaFile_DiffTMI << endl;
        return -8;
      }

      Sigma_DiffTMI.resize(NumRows);
      for (int i = 0; i < NumRows; ++i) {
        Sigma_DiffTMI[i] = MatrixTemp[i][0];
      }

      MatrixTemp.resize(0);

      cout << NumRows
           << " along track differenced TMI data points loaded from: "
           << "\n     " << DataFile_DiffTMI << endl;
      cout << "Corresponding data weights have been loaded from: " << "\n     "
           << DataWeightFile_DiffTMI << endl;
      cout << "Corresponding data sigma have been loaded from: " << "\n     "
           << SigmaFile_DiffTMI << endl;
    }
  }

  {
    if (IncludeDiffVector) {

      Flag = readmatrix(DataFile_DiffVector, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data file of along track "
                "differenced Vector data."
             << endl;
        return -7;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if (NumCols != 12) {
        cout << "An error occurred while reading the data file of along track "
                "differenced Vector data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file. "
                "(Expected: 12 columns)"
             << endl;
        cout << DataFile_DiffVector << endl;
        return -7;
      }

      LongitudeDiffVector1.resize(NumRows);
      LatitudeDiffVector1.resize(NumRows);
      RadiusDiffVector1.resize(NumRows);
      DiffBx1.resize(NumRows);
      DiffBy1.resize(NumRows);
      DiffBz1.resize(NumRows);
      LongitudeDiffVector2.resize(NumRows);
      LatitudeDiffVector2.resize(NumRows);
      RadiusDiffVector2.resize(NumRows);
      DiffBx2.resize(NumRows);
      DiffBy2.resize(NumRows);
      DiffBz2.resize(NumRows);

      for (int i = 0; i < NumRows; ++i) {
        LongitudeDiffVector1[i] = MatrixTemp[i][0];
        LatitudeDiffVector1[i] = MatrixTemp[i][1];
        RadiusDiffVector1[i] = MatrixTemp[i][2];
        DiffBx1[i] = MatrixTemp[i][3];
        DiffBy1[i] = MatrixTemp[i][4];
        DiffBz1[i] = MatrixTemp[i][5];
        LongitudeDiffVector2[i] = MatrixTemp[i][6];
        LatitudeDiffVector2[i] = MatrixTemp[i][7];
        RadiusDiffVector2[i] = MatrixTemp[i][8];
        DiffBx2[i] = MatrixTemp[i][9];
        DiffBy2[i] = MatrixTemp[i][10];
        DiffBz2[i] = MatrixTemp[i][11];
      }

      MatrixTemp.resize(0);

      Flag = readmatrix(DataWeightFile_DiffVector, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data weight file of along "
                "track differenced Vector data."
             << endl;
        return -8;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if ((static_cast<std::size_t>(NumRows) != LongitudeDiffVector1.size()) || (NumCols != 3)) {
        cout << "An error occurred while reading the data weight file of along "
                "track differenced Vector data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file."
             << endl;
        cout << DataWeightFile_DiffVector << endl;
        return -8;
      }

      DataWeight_DiffBx.resize(NumRows);
      DataWeight_DiffBy.resize(NumRows);
      DataWeight_DiffBz.resize(NumRows);
      for (int i = 0; i < NumRows; ++i) {
        DataWeight_DiffBx[i] = MatrixTemp[i][0];
        DataWeight_DiffBy[i] = MatrixTemp[i][1];
        DataWeight_DiffBz[i] = MatrixTemp[i][2];
      }

      MatrixTemp.resize(0);

      Flag = readmatrix(SigmaFile_DiffVector, MatrixTemp);

      if (Flag != 1) {
        cout << "An error occurred while reading the data sigma file of along "
                "track differenced Vector data."
             << endl;
        return -8;
      }

      NumRows = MatrixTemp.size();
      NumCols = MatrixTemp[0].size();

      if ((static_cast<std::size_t>(NumRows) != LongitudeDiffVector1.size()) || (NumCols != 3)) {
        cout << "An error occurred while reading the data sigma file of along "
                "track differenced Vector data."
             << endl;
        cout << "Please check the matrix dimensions in the input text file."
             << endl;
        cout << SigmaFile_DiffVector << endl;
        return -8;
      }

      Sigma_DiffBx.resize(NumRows);
      Sigma_DiffBy.resize(NumRows);
      Sigma_DiffBz.resize(NumRows);
      for (int i = 0; i < NumRows; ++i) {
        Sigma_DiffBx[i] = MatrixTemp[i][0];
        Sigma_DiffBy[i] = MatrixTemp[i][1];
        Sigma_DiffBz[i] = MatrixTemp[i][2];
      }

      MatrixTemp.resize(0);

      cout << NumRows
           << " along track differenced Vector data points loaded from: "
           << "\n     " << DataFile_DiffVector << endl;
      cout << "Corresponding data weights have been loaded from: " << "\n     "
           << DataWeightFile_DiffVector << endl;
      cout << "Corresponding data sigma have been loaded from: " << "\n     "
           << SigmaFile_DiffVector << endl;
    }
  }

  {

    Flag = readmatrix(InitialModelFile, MatrixTemp);

    if (Flag != 1) {
      cout << "An error occurred while reading the initial model file." << endl;
      return -9;
    }

    NumRows = MatrixTemp.size();
    NumCols = MatrixTemp[0].size();

    if (NumCols != 1) {
      cout << "An error occurred while reading the initial model file." << endl;
      cout << "Please check the matrix dimensions in the input text file. "
              "(Expected: 1 column)"
           << endl;
      cout << InitialModelFile << endl;
      return -9;
    }

    m0.resize(NumRows);

    for (int i = 0; i < NumRows; ++i) {
      m0[i] = MatrixTemp[i][0];
    }
    cout << "The initial model for the inversion have been loaded from: "
         << "\n     " << InitialModelFile << endl;
  }

  {

    Flag = readmatrix(ModelWeightFile, MatrixTemp);

    if (Flag != 1) {
      cout << "An error occurred while reading the initial model file." << endl;
      return -10;
    }

    NumRows = MatrixTemp.size();
    NumCols = MatrixTemp[0].size();

    if (NumCols != 1) {
      cout << "An error occurred while reading the initial model file." << endl;
      cout << "Please check the matrix dimensions in the input text file. "
              "(Expected: 1 column)"
           << endl;
      cout << ModelWeightFile << endl;
      return -10;
    }

    Weight_Model.resize(NumRows);

    for (int i = 0; i < NumRows; ++i) {
      Weight_Model[i] = MatrixTemp[i][0];
    }
    cout << "The model weight for the inversion have been loaded from: "
         << "\n     " << ModelWeightFile << endl;
  }

  {

    Flag = readmatrix(PriorModelFile, MatrixTemp);

    if (Flag != 1) {
      cout << "An error occurred while reading the initial model file." << endl;
      return -10;
    }

    NumRows = MatrixTemp.size();
    NumCols = MatrixTemp[0].size();

    if (NumCols != 1) {
      cout << "An error occurred while reading the initial model file." << endl;
      cout << "Please check the matrix dimensions in the input text file. "
              "(Expected: 1 column)"
           << endl;
      cout << PriorModelFile << endl;
      return -10;
    }

    Prior_Model.resize(NumRows);

    for (int i = 0; i < NumRows; ++i) {
      Prior_Model[i] = MatrixTemp[i][0];
    }
    cout << "The prior model for the inversion have been loaded from: "
         << "\n     " << PriorModelFile << endl;
  }

  cout << endl
       << "The modelling parameters and data have been loaded, no errors "
          "found.\n\n"
       << endl;
  return 1;
}

int ParametersCheck(

    const bool IncludeTMI, const bool IncludeVector, const bool IncludeDiffTMI,
    const bool IncludeDiffVector,

    double Lon0, double Lat0, double Theta0, int K_Max_Internal,
    int K_Max_External, int P_Max, int M_Max, double Max_Radius,
    double Min_Radius,

    double epsilon, double damping, int Max_Iteration, int RobustSwitch,
    double RobustConstant,

    std::vector<double> &LongitudeTMI, std::vector<double> &LatitudeTMI,
    std::vector<double> &RadiusTMI, std::vector<double> &IncTMI,
    std::vector<double> &DecTMI, std::vector<double> &TMI,
    std::vector<double> &DataWeight_TMI, std::vector<double> &Sigma_TMI,

    std::vector<double> &LongitudeVector, std::vector<double> &LatitudeVector,
    std::vector<double> &RadiusVector, std::vector<double> &Bx,
    std::vector<double> &By, std::vector<double> &Bz,
    std::vector<double> &DataWeight_Bx, std::vector<double> &DataWeight_By,
    std::vector<double> &DataWeight_Bz, std::vector<double> &Sigma_Bx,
    std::vector<double> &Sigma_By, std::vector<double> &Sigma_Bz,

    std::vector<double> &LongitudeDiffTMI1,
    std::vector<double> &LatitudeDiffTMI1, std::vector<double> &RadiusDiffTMI1,
    std::vector<double> &IncDiffTMI1, std::vector<double> &DecDiffTMI1,
    std::vector<double> &DiffTMI1, std::vector<double> &LongitudeDiffTMI2,
    std::vector<double> &LatitudeDiffTMI2, std::vector<double> &RadiusDiffTMI2,
    std::vector<double> &IncDiffTMI2, std::vector<double> &DecDiffTMI2,
    std::vector<double> &DiffTMI2, std::vector<double> &DataWeight_DiffTMI,
    std::vector<double> &Sigma_DiffTMI,

    std::vector<double> &LongitudeDiffVector1,
    std::vector<double> &LatitudeDiffVector1,
    std::vector<double> &RadiusDiffVector1, std::vector<double> &DiffBx1,
    std::vector<double> &DiffBy1, std::vector<double> &DiffBz1,
    std::vector<double> &LongitudeDiffVector2,
    std::vector<double> &LatitudeDiffVector2,
    std::vector<double> &RadiusDiffVector2, std::vector<double> &DiffBx2,
    std::vector<double> &DiffBy2, std::vector<double> &DiffBz2,
    std::vector<double> &DataWeight_DiffBx,
    std::vector<double> &DataWeight_DiffBy,
    std::vector<double> &DataWeight_DiffBz, std::vector<double> &Sigma_DiffBx,
    std::vector<double> &Sigma_DiffBy, std::vector<double> &Sigma_DiffBz) {
  check_geometry(Lon0, Lat0, Theta0, Min_Radius, Max_Radius);
  for (const auto *sigma : {&Sigma_TMI, &Sigma_Bx, &Sigma_By, &Sigma_Bz,
                            &Sigma_DiffTMI, &Sigma_DiffBx, &Sigma_DiffBy, &Sigma_DiffBz})
    for (double value : *sigma)
      if (!std::isfinite(value) || value <= 0.0)
        throw std::invalid_argument("Sigma values must be finite and positive.");
  for (const auto *weights : {&DataWeight_TMI, &DataWeight_Bx, &DataWeight_By,
                              &DataWeight_Bz, &DataWeight_DiffTMI, &DataWeight_DiffBx,
                              &DataWeight_DiffBy, &DataWeight_DiffBz})
    for (double value : *weights)
      if (!std::isfinite(value) || value < 0.0)
        throw std::invalid_argument("Data weights must be finite and nonnegative.");
  cout << "Checking the modelling parameters and data ...\n" << endl;

  if ((Lon0 < -180.0) || (Lon0 > 180.0)) {
    cout << "The longitude range of the centre point of the spherical cap is "
            "set incorrectly."
         << endl;
    cout << "The correct range is [-180.0 180.0]\n" << endl;
    return -11;
  }
  if ((Lat0 < -90.0) || (Lat0 > 90.0)) {
    cout << "The latitude range of the centre point of the spherical cap is "
            "set incorrectly."
         << endl;
    cout << "The correct range is [-90.0 90.0]\n" << endl;
    return -11;
  }
  if ((Theta0 <= 0.0) || (Theta0 > 90.0)) {
    cout << "The half-aperture angle of the spherical cap must be greater than "
            "zero."
         << endl;
    cout << "The correct range is (0.0 90.0]\n" << endl;
    return -11;
  }
  if ((K_Max_Internal < 0) || (K_Max_External < 0) || (P_Max < 0) ||
      (M_Max < 0)) {
    cout << "The maximum truncation indices of the basis series are set "
            "incorrectly."
         << endl;
    cout << "The correct range for K_Max_Internal, K_Max_External, P_Max and "
            "M_Max is: [0 +Inf]\n"
         << endl;
    return -11;
  }
  if (Max_Radius <= Min_Radius) {
    cout << "The radius of the upper cap must be greater than the lower cap.\n"
         << endl;
    return -11;
  }
  if (!std::isfinite(epsilon) || epsilon < 0.0) {
    cout << "The minimum tolerance for the optimization process is set "
            "incorrectly."
         << endl;
    cout << "The correct range is [0.0 +Inf]\n" << endl;
    return -11;
  }
  if (!std::isfinite(damping) || damping < 0.0) {
    cout
        << "The damping factor for the optimization process is set incorrectly."
        << endl;
    cout << "The correct range is [0.0 +Inf]\n" << endl;
    return -11;
  }
  if (Max_Iteration <= 0) {
    cout << "The maximum number of iterations for the optimization process is "
            "set incorrectly."
         << endl;
    cout << "The correct range is [1 +Inf]\n" << endl;
    return -11;
  }
  if ((RobustSwitch != 0) && (RobustSwitch != 1) && (RobustSwitch != 2)) {
    cout << "The Robust switch is set incorrectly." << endl;
    cout << "The correct choice is either 0(none), 1(Huber) or 2(Tukey).\n"
         << endl;
    return -11;
  }
  if (!std::isfinite(RobustConstant) || RobustConstant < 0 ||
      (RobustSwitch != 0 && RobustConstant == 0)) {
    cout << "The Robust constant is set incorrectly." << endl;
    cout << "The correct range is [0 +Inf]\n" << endl;
    return -11;
  }

  int EraseCount;
  if (IncludeTMI) {

    const int NumTMI = LongitudeTMI.size();
    if ((LatitudeTMI.size() != static_cast<std::size_t>(NumTMI)) || (RadiusTMI.size() != static_cast<std::size_t>(NumTMI)) ||
        (IncTMI.size() != static_cast<std::size_t>(NumTMI)) || (DecTMI.size() != static_cast<std::size_t>(NumTMI)) ||
        (TMI.size() != static_cast<std::size_t>(NumTMI)) || (DataWeight_TMI.size() != static_cast<std::size_t>(NumTMI)) ||
        (Sigma_TMI.size() != static_cast<std::size_t>(NumTMI))) {

      cout << "\nThe dimensions of TMI data (coordinates, data, data weight or "
              "data sigma) are inconsistent."
           << endl;
      cout << "Possible incorrect parameters:" << endl;
      cout << "	LongitudeTMI LatitudeTMI RadiusTMI IncTMI DecTMI TMI "
              "DataWeight_TMI Sigma_TMI\n"
           << endl;
      return -12;
    }

    EraseCount = 0;
    std::vector<double> LamdaTMI(NumTMI), ThetaTMI(NumTMI);
    SCHA_CT(LongitudeTMI, LatitudeTMI, Lon0, Lat0, 1, LamdaTMI, ThetaTMI);

    for (int i = 0; i < NumTMI; ++i) {

      if ((ThetaTMI[i - EraseCount] >= Theta0) ||
          (ThetaTMI[i - EraseCount] < 1.0e-3) ||
          (std::fabs(LatitudeTMI[i - EraseCount] - 90.0) < 1.0e-3)) {
        LongitudeTMI.erase(LongitudeTMI.begin() + i - EraseCount);
        LatitudeTMI.erase(LatitudeTMI.begin() + i - EraseCount);
        LamdaTMI.erase(LamdaTMI.begin() + i - EraseCount);
        ThetaTMI.erase(ThetaTMI.begin() + i - EraseCount);
        RadiusTMI.erase(RadiusTMI.begin() + i - EraseCount);
        TMI.erase(TMI.begin() + i - EraseCount);
        IncTMI.erase(IncTMI.begin() + i - EraseCount);
        DecTMI.erase(DecTMI.begin() + i - EraseCount);
        DataWeight_TMI.erase(DataWeight_TMI.begin() + i - EraseCount);
        Sigma_TMI.erase(Sigma_TMI.begin() + i - EraseCount);
        EraseCount++;
        continue;
      }
      if ((RadiusTMI[i - EraseCount] < Min_Radius) ||
          (RadiusTMI[i - EraseCount] > Max_Radius)) {
        LongitudeTMI.erase(LongitudeTMI.begin() + i - EraseCount);
        LatitudeTMI.erase(LatitudeTMI.begin() + i - EraseCount);
        RadiusTMI.erase(RadiusTMI.begin() + i - EraseCount);
        LamdaTMI.erase(LamdaTMI.begin() + i - EraseCount);
        ThetaTMI.erase(ThetaTMI.begin() + i - EraseCount);
        TMI.erase(TMI.begin() + i - EraseCount);
        IncTMI.erase(IncTMI.begin() + i - EraseCount);
        DecTMI.erase(DecTMI.begin() + i - EraseCount);
        DataWeight_TMI.erase(DataWeight_TMI.begin() + i - EraseCount);
        Sigma_TMI.erase(Sigma_TMI.begin() + i - EraseCount);
        EraseCount++;
        continue;
      }
    }

    LamdaTMI.resize(0);
    ThetaTMI.resize(0);

    if (EraseCount > 0) {
      cout << EraseCount
           << " TMI data points were outside the modelling area or " << endl;
      cout << "    too close to the centre of the modeling area or" << endl;
      cout << "    too close to the polar axis and have been deleted.\n"
           << endl;
    }
  }

  if (IncludeVector) {

    const int NumVector = LongitudeVector.size();
    if ((LatitudeVector.size() != static_cast<std::size_t>(NumVector)) ||
        (RadiusVector.size() != static_cast<std::size_t>(NumVector)) || (Bx.size() != static_cast<std::size_t>(NumVector)) ||
        (By.size() != static_cast<std::size_t>(NumVector)) || (Bz.size() != static_cast<std::size_t>(NumVector)) ||
        (DataWeight_Bx.size() != static_cast<std::size_t>(NumVector)) ||
        (DataWeight_By.size() != static_cast<std::size_t>(NumVector)) ||
        (DataWeight_Bz.size() != static_cast<std::size_t>(NumVector)) || (Sigma_Bx.size() != static_cast<std::size_t>(NumVector)) ||
        (Sigma_By.size() != static_cast<std::size_t>(NumVector)) || (Sigma_Bz.size() != static_cast<std::size_t>(NumVector))) {

      cout << "\nThe dimensions of Vector data (coordinates, data, or data "
              "weight) are inconsistent."
           << endl;
      cout << "Possible incorrect parameters:" << endl;
      cout << "	LongitudeVector LatitudeVector RadiusVector Bx By Bz" << endl;
      cout << "	DataWeight_Bx DataWeight_By DataWeight_Bz Sigma_Bx Sigma_By "
              "Sigma_Bz\n"
           << endl;
      return -12;
    }

    EraseCount = 0;
    std::vector<double> LamdaVector(NumVector), ThetaVector(NumVector);
    SCHA_CT(LongitudeVector, LatitudeVector, Lon0, Lat0, 1, LamdaVector,
            ThetaVector);

    for (int i = 0; i < NumVector; ++i) {
      if ((ThetaVector[i - EraseCount] >= Theta0) ||
          (ThetaVector[i - EraseCount] < 1.0e-3)) {
        LongitudeVector.erase(LongitudeVector.begin() + i - EraseCount);
        LatitudeVector.erase(LatitudeVector.begin() + i - EraseCount);
        RadiusVector.erase(RadiusVector.begin() + i - EraseCount);
        LamdaVector.erase(LamdaVector.begin() + i - EraseCount);
        ThetaVector.erase(ThetaVector.begin() + i - EraseCount);
        Bx.erase(Bx.begin() + i - EraseCount);
        By.erase(By.begin() + i - EraseCount);
        Bz.erase(Bz.begin() + i - EraseCount);
        DataWeight_Bx.erase(DataWeight_Bx.begin() + i - EraseCount);
        DataWeight_By.erase(DataWeight_By.begin() + i - EraseCount);
        DataWeight_Bz.erase(DataWeight_Bz.begin() + i - EraseCount);
        Sigma_Bx.erase(Sigma_Bx.begin() + i - EraseCount);
        Sigma_By.erase(Sigma_By.begin() + i - EraseCount);
        Sigma_Bz.erase(Sigma_Bz.begin() + i - EraseCount);
        EraseCount++;
        continue;
      }
      if ((RadiusVector[i - EraseCount] < Min_Radius) ||
          (RadiusVector[i - EraseCount] > Max_Radius)) {
        LongitudeVector.erase(LongitudeVector.begin() + i - EraseCount);
        LatitudeVector.erase(LatitudeVector.begin() + i - EraseCount);
        RadiusVector.erase(RadiusVector.begin() + i - EraseCount);
        LamdaVector.erase(LamdaVector.begin() + i - EraseCount);
        ThetaVector.erase(ThetaVector.begin() + i - EraseCount);
        Bx.erase(Bx.begin() + i - EraseCount);
        By.erase(By.begin() + i - EraseCount);
        Bz.erase(Bz.begin() + i - EraseCount);
        DataWeight_Bx.erase(DataWeight_Bx.begin() + i - EraseCount);
        DataWeight_By.erase(DataWeight_By.begin() + i - EraseCount);
        DataWeight_Bz.erase(DataWeight_Bz.begin() + i - EraseCount);
        Sigma_Bx.erase(Sigma_Bx.begin() + i - EraseCount);
        Sigma_By.erase(Sigma_By.begin() + i - EraseCount);
        Sigma_Bz.erase(Sigma_Bz.begin() + i - EraseCount);
        EraseCount++;
        continue;
      }
    }

    LamdaVector.resize(0);
    ThetaVector.resize(0);

    if (EraseCount > 0) {
      cout << EraseCount
           << " Vector data points were outside the modelling area or " << endl;
      cout << "    too close to the centre of the modeling area and have been "
              "deleted.\n"
           << endl;
    }
  }

  if (IncludeDiffTMI) {

    const int NumDiffTMI = LongitudeDiffTMI1.size();
    if ((LatitudeDiffTMI1.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (RadiusDiffTMI1.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (LongitudeDiffTMI2.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (LatitudeDiffTMI2.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (RadiusDiffTMI2.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (IncDiffTMI1.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (DecDiffTMI1.size() != static_cast<std::size_t>(NumDiffTMI)) || (DiffTMI1.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (IncDiffTMI2.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (DecDiffTMI2.size() != static_cast<std::size_t>(NumDiffTMI)) || (DiffTMI2.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (DataWeight_DiffTMI.size() != static_cast<std::size_t>(NumDiffTMI)) ||
        (Sigma_DiffTMI.size() != static_cast<std::size_t>(NumDiffTMI))) {

      cout << "\nThe dimensions of along track differenced TMI data "
              "(coordinates, data, or data weight) are inconsistent."
           << endl;
      cout << "Possible incorrect parameters:" << endl;
      cout << "	LongitudeDiffTMI1 LatitudeDiffTMI1 RadiusDiffTMI1" << endl;
      cout << "	IncDiffTMI1 DecDiffTMI1 DiffTMI1" << endl;
      cout << "	LongitudeDiffTMI2 LatitudeDiffTMI2 RadiusDiffTMI2" << endl;
      cout << "	IncDiffTMI2 DecDiffTMI2 DiffTMI2" << endl;
      cout << "	DataWeight_DiffTMI Sigma_DiffTMI\n" << endl;
      return -12;
    }

    EraseCount = 0;
    std::vector<double> LamdaDiffTMI1(NumDiffTMI), ThetaDiffTMI1(NumDiffTMI);
    std::vector<double> LamdaDiffTMI2(NumDiffTMI), ThetaDiffTMI2(NumDiffTMI);
    SCHA_CT(LongitudeDiffTMI1, LatitudeDiffTMI1, Lon0, Lat0, 1, LamdaDiffTMI1,
            ThetaDiffTMI1);
    SCHA_CT(LongitudeDiffTMI2, LatitudeDiffTMI2, Lon0, Lat0, 1, LamdaDiffTMI2,
            ThetaDiffTMI2);

    for (int i = 0; i < NumDiffTMI; ++i) {
      if ((ThetaDiffTMI1[i - EraseCount] >= Theta0) ||
          (ThetaDiffTMI1[i - EraseCount] < 1.0e-3) ||
          (ThetaDiffTMI2[i - EraseCount] >= Theta0) ||
          (ThetaDiffTMI2[i - EraseCount] < 1.0e-3) ||
          (std::fabs(LatitudeDiffTMI1[i - EraseCount] - 90.0) < 1.0e-3) ||
          (std::fabs(LatitudeDiffTMI2[i - EraseCount] - 90.0) < 1.0e-3)) {
        LongitudeDiffTMI1.erase(LongitudeDiffTMI1.begin() + i - EraseCount);
        LatitudeDiffTMI1.erase(LatitudeDiffTMI1.begin() + i - EraseCount);
        RadiusDiffTMI1.erase(RadiusDiffTMI1.begin() + i - EraseCount);
        LamdaDiffTMI1.erase(LamdaDiffTMI1.begin() + i - EraseCount);
        ThetaDiffTMI1.erase(ThetaDiffTMI1.begin() + i - EraseCount);
        DiffTMI1.erase(DiffTMI1.begin() + i - EraseCount);
        IncDiffTMI1.erase(IncDiffTMI1.begin() + i - EraseCount);
        DecDiffTMI1.erase(DecDiffTMI1.begin() + i - EraseCount);
        LongitudeDiffTMI2.erase(LongitudeDiffTMI2.begin() + i - EraseCount);
        LatitudeDiffTMI2.erase(LatitudeDiffTMI2.begin() + i - EraseCount);
        RadiusDiffTMI2.erase(RadiusDiffTMI2.begin() + i - EraseCount);
        LamdaDiffTMI2.erase(LamdaDiffTMI2.begin() + i - EraseCount);
        ThetaDiffTMI2.erase(ThetaDiffTMI2.begin() + i - EraseCount);
        DiffTMI2.erase(DiffTMI2.begin() + i - EraseCount);
        IncDiffTMI2.erase(IncDiffTMI2.begin() + i - EraseCount);
        DecDiffTMI2.erase(DecDiffTMI2.begin() + i - EraseCount);
        DataWeight_DiffTMI.erase(DataWeight_DiffTMI.begin() + i - EraseCount);
        Sigma_DiffTMI.erase(Sigma_DiffTMI.begin() + i - EraseCount);
        EraseCount++;
        continue;
      }
      if ((RadiusDiffTMI1[i - EraseCount] < Min_Radius) ||
          (RadiusDiffTMI1[i - EraseCount] > Max_Radius) ||
          (RadiusDiffTMI2[i - EraseCount] < Min_Radius) ||
          (RadiusDiffTMI2[i - EraseCount] > Max_Radius)) {
        LongitudeDiffTMI1.erase(LongitudeDiffTMI1.begin() + i - EraseCount);
        LatitudeDiffTMI1.erase(LatitudeDiffTMI1.begin() + i - EraseCount);
        RadiusDiffTMI1.erase(RadiusDiffTMI1.begin() + i - EraseCount);
        LamdaDiffTMI1.erase(LamdaDiffTMI1.begin() + i - EraseCount);
        ThetaDiffTMI1.erase(ThetaDiffTMI1.begin() + i - EraseCount);
        DiffTMI1.erase(DiffTMI1.begin() + i - EraseCount);
        IncDiffTMI1.erase(IncDiffTMI1.begin() + i - EraseCount);
        DecDiffTMI1.erase(DecDiffTMI1.begin() + i - EraseCount);
        LongitudeDiffTMI2.erase(LongitudeDiffTMI2.begin() + i - EraseCount);
        LatitudeDiffTMI2.erase(LatitudeDiffTMI2.begin() + i - EraseCount);
        RadiusDiffTMI2.erase(RadiusDiffTMI2.begin() + i - EraseCount);
        LamdaDiffTMI2.erase(LamdaDiffTMI2.begin() + i - EraseCount);
        ThetaDiffTMI2.erase(ThetaDiffTMI2.begin() + i - EraseCount);
        DiffTMI2.erase(DiffTMI2.begin() + i - EraseCount);
        IncDiffTMI2.erase(IncDiffTMI2.begin() + i - EraseCount);
        DecDiffTMI2.erase(DecDiffTMI2.begin() + i - EraseCount);
        DataWeight_DiffTMI.erase(DataWeight_DiffTMI.begin() + i - EraseCount);
        Sigma_DiffTMI.erase(Sigma_DiffTMI.begin() + i - EraseCount);
        EraseCount++;
        continue;
      }
    }

    LamdaDiffTMI1.resize(0);
    ThetaDiffTMI1.resize(0);
    LamdaDiffTMI2.resize(0);
    ThetaDiffTMI2.resize(0);

    if (EraseCount > 0) {
      cout << EraseCount
           << " along track differenced TMI data points were outside the "
              "modelling area or "
           << endl;
      cout << "    too close to the centre of the modeling area or" << endl;
      cout << "    too close to the polar axis and have been deleted.\n"
           << endl;
    }
  }

  if (IncludeDiffVector) {

    const int NumDiffVector = LongitudeDiffVector1.size();
    if ((LatitudeDiffVector1.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (RadiusDiffVector1.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (LongitudeDiffVector2.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (LatitudeDiffVector2.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (RadiusDiffVector2.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (DiffBx1.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (DiffBy1.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (DiffBz1.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (DiffBx2.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (DiffBy2.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (DiffBz2.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (DataWeight_DiffBx.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (DataWeight_DiffBy.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (DataWeight_DiffBz.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (Sigma_DiffBx.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (Sigma_DiffBy.size() != static_cast<std::size_t>(NumDiffVector)) ||
        (Sigma_DiffBz.size() != static_cast<std::size_t>(NumDiffVector))) {

      cout << "\nThe dimensions of along track differenced vector data "
              "(coordinates, data, or data weight) are inconsistent."
           << endl;
      cout << "Possible incorrect parameters:" << endl;
      cout << "	LongitudeDiffVector1 LatitudeDiffVector1 RadiusDiffVector1"
           << endl;
      cout << "	DiffBx1 DiffBy1 DiffBz1" << endl;
      cout << "	LongitudeDiffVector2 LatitudeDiffVector2 RadiusDiffVector2"
           << endl;
      cout << "	DiffBx2 DiffBy2 DiffBz2" << endl;
      cout << "	DataWeight_DiffBx DataWeight_DiffBy DataWeight_DiffBz" << endl;
      cout << "	Sigma_DiffBx Sigma_DiffBy Sigma_DiffBz\n" << endl;
      return -12;
    }

    EraseCount = 0;
    std::vector<double> LamdaDiffVector1(NumDiffVector),
        ThetaDiffVector1(NumDiffVector);
    std::vector<double> LamdaDiffVector2(NumDiffVector),
        ThetaDiffVector2(NumDiffVector);
    SCHA_CT(LongitudeDiffVector1, LatitudeDiffVector1, Lon0, Lat0, 1,
            LamdaDiffVector1, ThetaDiffVector1);
    SCHA_CT(LongitudeDiffVector2, LatitudeDiffVector2, Lon0, Lat0, 1,
            LamdaDiffVector2, ThetaDiffVector2);

    for (int i = 0; i < NumDiffVector; ++i) {
      if ((ThetaDiffVector1[i - EraseCount] >= Theta0) ||
          (ThetaDiffVector1[i - EraseCount] < 1.0e-3) ||
          (ThetaDiffVector2[i - EraseCount] >= Theta0) ||
          (ThetaDiffVector2[i - EraseCount] < 1.0e-3)) {
        LongitudeDiffVector1.erase(LongitudeDiffVector1.begin() + i -
                                   EraseCount);
        LatitudeDiffVector1.erase(LatitudeDiffVector1.begin() + i - EraseCount);
        RadiusDiffVector1.erase(RadiusDiffVector1.begin() + i - EraseCount);
        LamdaDiffVector1.erase(LamdaDiffVector1.begin() + i - EraseCount);
        ThetaDiffVector1.erase(ThetaDiffVector1.begin() + i - EraseCount);
        DiffBx1.erase(DiffBx1.begin() + i - EraseCount);
        DiffBy1.erase(DiffBy1.begin() + i - EraseCount);
        DiffBz1.erase(DiffBz1.begin() + i - EraseCount);
        LongitudeDiffVector2.erase(LongitudeDiffVector2.begin() + i -
                                   EraseCount);
        LatitudeDiffVector2.erase(LatitudeDiffVector2.begin() + i - EraseCount);
        RadiusDiffVector2.erase(RadiusDiffVector2.begin() + i - EraseCount);
        LamdaDiffVector2.erase(LamdaDiffVector2.begin() + i - EraseCount);
        ThetaDiffVector2.erase(ThetaDiffVector2.begin() + i - EraseCount);
        DiffBx2.erase(DiffBx2.begin() + i - EraseCount);
        DiffBy2.erase(DiffBy2.begin() + i - EraseCount);
        DiffBz2.erase(DiffBz2.begin() + i - EraseCount);
        DataWeight_DiffBx.erase(DataWeight_DiffBx.begin() + i - EraseCount);
        DataWeight_DiffBy.erase(DataWeight_DiffBy.begin() + i - EraseCount);
        DataWeight_DiffBz.erase(DataWeight_DiffBz.begin() + i - EraseCount);
        Sigma_DiffBx.erase(Sigma_DiffBx.begin() + i - EraseCount);
        Sigma_DiffBy.erase(Sigma_DiffBy.begin() + i - EraseCount);
        Sigma_DiffBz.erase(Sigma_DiffBz.begin() + i - EraseCount);
        EraseCount++;
        continue;
      }
      if ((RadiusDiffVector1[i - EraseCount] < Min_Radius) ||
          (RadiusDiffVector1[i - EraseCount] > Max_Radius) ||
          (RadiusDiffVector2[i - EraseCount] < Min_Radius) ||
          (RadiusDiffVector2[i - EraseCount] > Max_Radius)) {
        LongitudeDiffVector1.erase(LongitudeDiffVector1.begin() + i -
                                   EraseCount);
        LatitudeDiffVector1.erase(LatitudeDiffVector1.begin() + i - EraseCount);
        RadiusDiffVector1.erase(RadiusDiffVector1.begin() + i - EraseCount);
        LamdaDiffVector1.erase(LamdaDiffVector1.begin() + i - EraseCount);
        ThetaDiffVector1.erase(ThetaDiffVector1.begin() + i - EraseCount);
        DiffBx1.erase(DiffBx1.begin() + i - EraseCount);
        DiffBy1.erase(DiffBy1.begin() + i - EraseCount);
        DiffBz1.erase(DiffBz1.begin() + i - EraseCount);
        LongitudeDiffVector2.erase(LongitudeDiffVector2.begin() + i -
                                   EraseCount);
        LatitudeDiffVector2.erase(LatitudeDiffVector2.begin() + i - EraseCount);
        RadiusDiffVector2.erase(RadiusDiffVector2.begin() + i - EraseCount);
        LamdaDiffVector2.erase(LamdaDiffVector2.begin() + i - EraseCount);
        ThetaDiffVector2.erase(ThetaDiffVector2.begin() + i - EraseCount);
        DiffBx2.erase(DiffBx2.begin() + i - EraseCount);
        DiffBy2.erase(DiffBy2.begin() + i - EraseCount);
        DiffBz2.erase(DiffBz2.begin() + i - EraseCount);
        DataWeight_DiffBx.erase(DataWeight_DiffBx.begin() + i - EraseCount);
        DataWeight_DiffBy.erase(DataWeight_DiffBy.begin() + i - EraseCount);
        DataWeight_DiffBz.erase(DataWeight_DiffBz.begin() + i - EraseCount);
        Sigma_DiffBx.erase(Sigma_DiffBx.begin() + i - EraseCount);
        Sigma_DiffBy.erase(Sigma_DiffBy.begin() + i - EraseCount);
        Sigma_DiffBz.erase(Sigma_DiffBz.begin() + i - EraseCount);
        EraseCount++;
        continue;
      }
    }

    LamdaDiffVector1.resize(0);
    ThetaDiffVector1.resize(0);
    LamdaDiffVector2.resize(0);
    ThetaDiffVector2.resize(0);

    if (EraseCount > 0) {
      cout << EraseCount
           << " along track differenced vector data points were outside the "
              "modelling area or "
           << endl;
      cout << "    too close to the centre of the modeling area and have been "
              "deleted.\n"
           << endl;
    }
  }

  cout << "Parameter check completed, no errors found.\n" << endl;
  return 1;
}

inline bool SaveModellingResults(

    const string &ReconsTMIFile, const string &ReconsVectorFile,
    const string &ReconsDiffTMIFile, const string &ReconsDiffVectorFile,

    const string &ModellingResultsFile, const string &ModelSpectrumFile,

    const bool IncludeTMI, const bool IncludeVector, const bool IncludeDiffTMI,
    const bool IncludeDiffVector,

    double Lon0, double Lat0, double Theta0, int K_Max_Internal,
    int K_Max_External, int P_Max, int M_Max, double Max_Radius,
    double Min_Radius,

    std::vector<double> LongitudeTMI, std::vector<double> LatitudeTMI,
    std::vector<double> RadiusTMI, std::vector<double> TMI,
    std::vector<double> ReconsTMI,

    std::vector<double> LongitudeVector, std::vector<double> LatitudeVector,
    std::vector<double> RadiusVector, std::vector<double> Bx,
    std::vector<double> By, std::vector<double> Bz,
    std::vector<double> ReconsBx, std::vector<double> ReconsBy,
    std::vector<double> ReconsBz,

    std::vector<double> LongitudeDiffTMI1, std::vector<double> LatitudeDiffTMI1,
    std::vector<double> RadiusDiffTMI1, std::vector<double> DiffTMI1,
    std::vector<double> ReconsDiffTMI1, std::vector<double> LongitudeDiffTMI2,
    std::vector<double> LatitudeDiffTMI2, std::vector<double> RadiusDiffTMI2,
    std::vector<double> DiffTMI2, std::vector<double> ReconsDiffTMI2,

    std::vector<double> LongitudeDiffVector1,
    std::vector<double> LatitudeDiffVector1,
    std::vector<double> RadiusDiffVector1, std::vector<double> DiffBx1,
    std::vector<double> DiffBy1, std::vector<double> DiffBz1,
    std::vector<double> ReconsDiffBx1, std::vector<double> ReconsDiffBy1,
    std::vector<double> ReconsDiffBz1, std::vector<double> LongitudeDiffVector2,
    std::vector<double> LatitudeDiffVector2,
    std::vector<double> RadiusDiffVector2, std::vector<double> DiffBx2,
    std::vector<double> DiffBy2, std::vector<double> DiffBz2,
    std::vector<double> ReconsDiffBx2, std::vector<double> ReconsDiffBy2,
    std::vector<double> ReconsDiffBz2, std::vector<double> Coefficients,
    std::vector<double> n, std::vector<double> Spectrum) {

  const int NumTMI = LongitudeTMI.size();
  const int NumVector = LongitudeVector.size();
  const int NumDiffTMI = LongitudeDiffTMI1.size();
  const int NumDiffVector = LongitudeDiffVector1.size();
  int Flag = 0;

  cout << "\nSaveing the modelling results and coefficients ..." << endl
       << endl;

  if (IncludeTMI) {
    vector<vector<double>> MatrixSave(NumTMI, vector<double>(6));

    for (int i = 0; i < NumTMI; ++i) {
      MatrixSave[i][0] = LongitudeTMI[i];
      MatrixSave[i][1] = LatitudeTMI[i];
      MatrixSave[i][2] = RadiusTMI[i];
      MatrixSave[i][3] = TMI[i];
      MatrixSave[i][4] = ReconsTMI[i];
      MatrixSave[i][5] = TMI[i] - ReconsTMI[i];
    }
    Flag = writematrix(ReconsTMIFile, MatrixSave);
    MatrixSave.resize(0);

    if (Flag == 1) {
      cout << "Reconstructed TMI data saved in:" << "\n     " << ReconsTMIFile
           << endl;
    } else {
      cout << "An error occurred while saving the file." << endl;
      return false;
    }
  }

  if (IncludeVector) {
    vector<vector<double>> MatrixSave(NumVector, vector<double>(12));

    for (int i = 0; i < NumVector; ++i) {
      MatrixSave[i][0] = LongitudeVector[i];
      MatrixSave[i][1] = LatitudeVector[i];
      MatrixSave[i][2] = RadiusVector[i];
      MatrixSave[i][3] = Bx[i];
      MatrixSave[i][4] = By[i];
      MatrixSave[i][5] = Bz[i];
      MatrixSave[i][6] = ReconsBx[i];
      MatrixSave[i][7] = ReconsBy[i];
      MatrixSave[i][8] = ReconsBz[i];
      MatrixSave[i][9] = Bx[i] - ReconsBx[i];
      MatrixSave[i][10] = By[i] - ReconsBy[i];
      MatrixSave[i][11] = Bz[i] - ReconsBz[i];
    }
    Flag = writematrix(ReconsVectorFile, MatrixSave);
    MatrixSave.resize(0);

    if (Flag == 1) {
      cout << "Reconstructed vector data saved in:" << "\n     "
           << ReconsVectorFile << endl;
    } else {
      cout << "An error occurred while saving the file." << endl;
      return false;
    }
  }

  if (IncludeDiffTMI) {
    vector<vector<double>> MatrixSave(NumDiffTMI, vector<double>(12));

    for (int i = 0; i < NumDiffTMI; ++i) {
      MatrixSave[i][0] = LongitudeDiffTMI1[i];
      MatrixSave[i][1] = LatitudeDiffTMI1[i];
      MatrixSave[i][2] = RadiusDiffTMI1[i];
      MatrixSave[i][3] = DiffTMI1[i];
      MatrixSave[i][4] = ReconsDiffTMI1[i];
      MatrixSave[i][5] = DiffTMI1[i] - ReconsDiffTMI1[i];
      MatrixSave[i][6] = LongitudeDiffTMI2[i];
      MatrixSave[i][7] = LatitudeDiffTMI2[i];
      MatrixSave[i][8] = RadiusDiffTMI2[i];
      MatrixSave[i][9] = DiffTMI2[i];
      MatrixSave[i][10] = ReconsDiffTMI2[i];
      MatrixSave[i][11] = DiffTMI2[i] - ReconsDiffTMI2[i];
    }
    Flag = writematrix(ReconsDiffTMIFile, MatrixSave);
    MatrixSave.resize(0);

    if (Flag == 1) {
      cout << "Reconstructed along track differenced TMI data saved in:"
           << "\n     " << ReconsDiffTMIFile << endl;
    } else {
      cout << "An error occurred while saving the file." << endl;
      return false;
    }
  }

  if (IncludeDiffVector) {
    vector<vector<double>> MatrixSave(NumDiffVector, vector<double>(24));

    for (int i = 0; i < NumDiffVector; ++i) {
      MatrixSave[i][0] = LongitudeDiffVector1[i];
      MatrixSave[i][1] = LatitudeDiffVector1[i];
      MatrixSave[i][2] = RadiusDiffVector1[i];
      MatrixSave[i][3] = DiffBx1[i];
      MatrixSave[i][4] = DiffBy1[i];
      MatrixSave[i][5] = DiffBz1[i];
      MatrixSave[i][6] = ReconsDiffBx1[i];
      MatrixSave[i][7] = ReconsDiffBy1[i];
      MatrixSave[i][8] = ReconsDiffBz1[i];
      MatrixSave[i][9] = DiffBx1[i] - ReconsDiffBx1[i];
      MatrixSave[i][10] = DiffBy1[i] - ReconsDiffBy1[i];
      MatrixSave[i][11] = DiffBz1[i] - ReconsDiffBz1[i];
      MatrixSave[i][12] = LongitudeDiffVector2[i];
      MatrixSave[i][13] = LatitudeDiffVector2[i];
      MatrixSave[i][14] = RadiusDiffVector2[i];
      MatrixSave[i][15] = DiffBx2[i];
      MatrixSave[i][16] = DiffBy2[i];
      MatrixSave[i][17] = DiffBz2[i];
      MatrixSave[i][18] = ReconsDiffBx2[i];
      MatrixSave[i][19] = ReconsDiffBy2[i];
      MatrixSave[i][20] = ReconsDiffBz2[i];
      MatrixSave[i][21] = DiffBx2[i] - ReconsDiffBx2[i];
      MatrixSave[i][22] = DiffBy2[i] - ReconsDiffBy2[i];
      MatrixSave[i][23] = DiffBz2[i] - ReconsDiffBz2[i];
    }
    Flag = writematrix(ReconsDiffVectorFile, MatrixSave);
    MatrixSave.resize(0);

    if (Flag == 1) {
      cout << "Reconstructed along track differenced vector data saved in:"
           << "\n     " << ReconsDiffVectorFile << endl;
    } else {
      cout << "An error occurred while saving the file." << endl;
      return false;
    }
  }

  ofstream fp;
  fp.open(ModellingResultsFile, ios::out);
  if (!fp) {
    cout << "Fail to open file: " << ModellingResultsFile << "\n";
    cout << "An error occurred while saving the file." << endl;
      return false;
  } else {
    fp << fixed << setprecision(15);

    fp << Lon0 << endl;
    fp << Lat0 << endl;
    fp << Theta0 << endl;
    fp << K_Max_Internal << endl;
    fp << K_Max_External << endl;
    fp << P_Max << endl;
    fp << M_Max << endl;
    fp << Max_Radius << endl;
    fp << Min_Radius << endl;

    for (size_t i = 0; i < Coefficients.size(); ++i) {
      fp << Coefficients[i] << endl;
    }
    fp.close();
    if (!fp)
      return false;

    cout << "Modelling parameters and RSCHA coefficients saved in:" << "\n     "
         << ModellingResultsFile << endl;
  }

  fp.open(ModelSpectrumFile, ios::out);
  if (!fp) {
    cout << "Fail to open file: " << ModelSpectrumFile << "\n";
    cout << "An error occurred while saving the file." << endl;
      return false;
  } else {
    fp << fixed << setprecision(15);

    const int NumSpectrum = n.size();

    for (int i = 0; i < NumSpectrum; ++i) {
      fp << n[i] << ",";
      fp << Spectrum[i] << endl;
    }
    fp.close();
    if (!fp)
      return false;

    cout << "RSCHA spectrum of the built model saved in:" << "\n     "
         << ModelSpectrumFile << endl;
  }
  return true;
}

struct CoefPacks {
  std::vector<std::vector<double>> Gikm, Hikm, Gekm, Hekm;
  std::vector<std::vector<double>> Gpm, Hpm;
  std::vector<double> G0m, H0m;
};

inline CoefPacks Sort_Coefficient(const std::vector<double> &Coefficient,
                                  int K_Max_Internal, int K_Max_External,
                                  int P_Max, int M_Max) {
  CoefPacks out;
  out.Gikm.assign(K_Max_Internal + 1,
                  std::vector<double>(K_Max_Internal + 1, 0.0));
  out.Hikm = out.Gikm;
  out.Gekm.assign(K_Max_External + 1,
                  std::vector<double>(K_Max_External + 1, 0.0));
  out.Hekm = out.Gekm;
  out.Gpm.assign(P_Max + 1, std::vector<double>(P_Max + 1, 0.0));
  out.Hpm = out.Gpm;
  out.G0m.assign(M_Max + 1, 0.0);
  out.H0m.assign(M_Max + 1, 0.0);

  size_t Num = 0;

  for (int k = 0; k <= K_Max_Internal; ++k)
    for (int m = 0; m <= k; ++m)
      out.Gikm[k][m] = Coefficient.at(Num++);

  for (int k = 0; k <= K_Max_Internal; ++k)
    for (int m = 0; m <= k; ++m)
      out.Hikm[k][m] = Coefficient.at(Num++);

  for (int k = 0; k <= K_Max_External; ++k)
    for (int m = 0; m <= k; ++m)
      out.Gekm[k][m] = Coefficient.at(Num++);

  for (int k = 0; k <= K_Max_External; ++k)
    for (int m = 0; m <= k; ++m)
      out.Hekm[k][m] = Coefficient.at(Num++);

  for (int p = 1; p <= P_Max; ++p)
    for (int m = 0; m <= P_Max; ++m)
      out.Gpm[p][m] = Coefficient.at(Num++);

  for (int p = 1; p <= P_Max; ++p)
    for (int m = 0; m <= P_Max; ++m)
      out.Hpm[p][m] = Coefficient.at(Num++);

  for (int m = 0; m <= M_Max; ++m)
    out.G0m[m] = Coefficient.at(Num++);

  for (int m = 0; m <= M_Max; ++m)
    out.H0m[m] = Coefficient.at(Num++);

  return out;
}

inline double Pnmdn(double m, double nk, double Theta0) {
  const double theta0 = Theta0 * (4.0 * std::atan(1.0)) / 180.0;
  return Legendre::BoundaryDegreeDerivativeLegacySchmidt(m, nk, theta0);
}

inline double Norm_Legendre(double m, double nk, double Theta0) {
  const double theta0 = Theta0 * (4.0 * std::atan(1.0)) / 180.0;
  return Legendre::CapNormLegacySchmidt(m, nk, theta0);
}

inline void RSCHA_Spectrum(vector<double> Coefficient, int K_Max_Internal,
                           int K_Max_External, int P_Max, int M_Max,
                           vector<double> nk, double Theta0, double a,
                           double Radius, vector<double> &n,
                           vector<double> &Spectrum) {
  if (MPIParallel::rank() == 0)
    cout << "\nEstimating the RSCHA spectrum of the built model ...     ";

  const double dn = 180.0 / Theta0;

  auto packs = Sort_Coefficient(Coefficient, K_Max_Internal, K_Max_External,
                                P_Max, M_Max);

  const int K_Max = std::max(K_Max_Internal, K_Max_External);

  std::vector<std::vector<double>> Gikm(K_Max + 1,
                                        std::vector<double>(K_Max + 1, 0.0));
  std::vector<std::vector<double>> Hikm = Gikm;
  std::vector<std::vector<double>> Gekm = Gikm;
  std::vector<std::vector<double>> Hekm = Gikm;

  for (int k = 0; k <= K_Max_Internal; ++k) {
    for (int m = 0; m <= k; ++m) {
      Gikm[k][m] = packs.Gikm[k][m];
      Hikm[k][m] = packs.Hikm[k][m];
    }
  }

  for (int k = 0; k <= K_Max_External; ++k) {
    for (int m = 0; m <= k; ++m) {
      Gekm[k][m] = packs.Gekm[k][m];
      Hekm[k][m] = packs.Hekm[k][m];
    }
  }

  int Position = K_Max * (1 + K_Max) / 2 + 0;
  double Max_nk = nk.at(Position);

  const int Num_n = static_cast<int>(std::ceil(Max_nk / dn));
  n.resize(Num_n);
  for (int i = 0; i < Num_n; ++i)
    n[i] = dn * (i + 1) - dn / 2.0;

  Spectrum.resize(Num_n, 0.0);

  for (int k = 0; k <= K_Max; ++k) {
    for (int m = 0; m <= k; ++m) {

      Position = k * (1 + k) / 2 + m;
      if (Position % MPIParallel::size() != MPIParallel::rank())
        continue;
      const double nk_Local = nk.at(Position);
      int Location = static_cast<int>(std::floor(nk_Local / dn));
      if (Location < 0)
        Location = 0;
      if (Location >= Num_n)
        Location = Num_n - 1;

      const double Factor = (m == 0) ? (1.0 / (1.0 - cosd(Theta0)))
                                     : (1.0 / (2.0 * (1.0 - cosd(Theta0))));

      const double Norm_Pnm = Norm_Legendre((double)m, nk_Local, Theta0);

      const double term_int =
          (nk_Local + 1.0) * (2.0 * nk_Local + 1.0) *
          std::pow(a / Radius, 2.0 * nk_Local + 4.0) *
          (Gikm[k][m] * Gikm[k][m] + Hikm[k][m] * Hikm[k][m]);

      const double term_ext =
          nk_Local * (2.0 * nk_Local + 1.0) *
          std::pow(Radius / a, 2.0 * nk_Local - 2.0) *
          (Gekm[k][m] * Gekm[k][m] + Hekm[k][m] * Hekm[k][m]);

      Spectrum[Location] += Factor * (term_int + term_ext) * Norm_Pnm;
    }
  }

  MPIParallel::allreduce_sum(Spectrum);

  const double scale = Theta0 / 180.0;
  for (double &v : Spectrum)
    v *= scale;

  if (MPIParallel::rank() == 0)
    cout << "Completed!\n" << endl;
}

}
