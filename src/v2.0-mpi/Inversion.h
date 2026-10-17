#pragma once
#include "MPIParallel.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <stdexcept>
#include <limits>

namespace Inversion {

enum class SolveStatus { converged, data_fit, iteration_limit };

class DenseMatrix {
public:
  DenseMatrix() = default;
  DenseMatrix(int rows, int columns) : columns_(columns) {
    if (rows < 0 || columns <= 0)
      throw std::invalid_argument("Invalid matrix dimensions.");
    values_.assign(static_cast<std::size_t>(rows) * columns, 0.0);
  }

  double *operator[](int row) {
    return values_.data() + static_cast<std::size_t>(row) * columns_;
  }
  const double *operator[](int row) const {
    return values_.data() + static_cast<std::size_t>(row) * columns_;
  }
  int rows() const {
    return columns_ == 0 ? 0 : static_cast<int>(values_.size() / columns_);
  }
  int columns() const { return columns_; }
  std::vector<double> &values() { return values_; }
  const std::vector<double> &values() const { return values_; }
private:
  int columns_ = 0;
  std::vector<double> values_;
};

inline SolveStatus PreCond_BICGSTAB(
    int NumData, int NumCoefficients, const DenseMatrix &A,
    const std::vector<int> &GlobalRows,
    std::vector<double> &m_predict, const std::vector<double> &d,
    double epsilon, double damping, int Max_Iteration,
    const std::vector<double> &Weight_Data,
    const std::vector<double> &Weight_Model,
    const std::vector<double> &Prior_Model, const std::vector<double> &Sigma,
    int RobustSwitch, double RobustConstant, std::vector<double> &m,
    std::vector<double> &Misfit, std::vector<double> &d_Predict) {
  if (NumData <= 0 || NumCoefficients <= 0 || Max_Iteration <= 0)
    throw std::invalid_argument("Inversion requires data, coefficients, and positive iteration count.");
  const std::size_t columns = static_cast<std::size_t>(NumCoefficients);
  if (m_predict.size() != columns || Weight_Model.size() != columns ||
      Prior_Model.size() != columns || m.size() != columns ||
      Misfit.size() < static_cast<std::size_t>(Max_Iteration) ||
      d_Predict.size() != static_cast<std::size_t>(NumData) ||
      d.size() != Sigma.size() || d.size() != Weight_Data.size() ||
      A.columns() != NumCoefficients ||
      d.size() != static_cast<std::size_t>(A.rows()) ||
      GlobalRows.size() != d.size())
    throw std::invalid_argument("Inconsistent inversion vector dimensions.");
  for (int row : GlobalRows)
    if (row < 0 || row >= NumData)
      throw std::invalid_argument("Invalid global data row.");
  for (std::size_t i = 0; i < d.size(); ++i)
    if (!std::isfinite(d[i]) || !std::isfinite(Sigma[i]) || Sigma[i] <= 0.0 ||
        !std::isfinite(Weight_Data[i]) || Weight_Data[i] < 0.0)
      throw std::invalid_argument("Data must be finite, sigma positive, and weights nonnegative.");
  for (std::size_t j = 0; j < columns; ++j)
    if (!std::isfinite(m_predict[j]) || !std::isfinite(Prior_Model[j]) ||
        !std::isfinite(Weight_Model[j]) || Weight_Model[j] < 0.0)
      throw std::invalid_argument("Invalid initial model, prior model, or model weight.");
  if (!std::isfinite(epsilon) || epsilon < 0.0 ||
      !std::isfinite(damping) || damping < 0.0 ||
      RobustSwitch < 0 || RobustSwitch > 2 ||
      !std::isfinite(RobustConstant) || RobustConstant < 0.0 ||
      (RobustSwitch != 0 && RobustConstant == 0.0))
    throw std::invalid_argument("Invalid inversion controls.");
  const double normal_tolerance = 64.0 * std::numeric_limits<double>::epsilon();
  const bool is_root = MPIParallel::rank() == 0;
  const int NumLocal = A.rows();
  std::vector<double> TempVecRows(NumLocal), TempVecCols(columns);
  std::vector<double> DampingFactor1(columns), DampingFactor2(columns);
  std::vector<double> F(columns), newd(columns), r0(columns), r_new(columns);
  std::vector<double> pj(columns), m_old(columns), Apj(columns), sj(columns), Asj(columns);
  std::vector<double> ReW(NumLocal), LocalPredict(NumLocal);
  MPIParallel::phase_begin("inversion partition validation");
  double started = MPIParallel::time();
  MPIParallel::IndexedGather prediction_gather(NumData, GlobalRows);
  MPIParallel::phase_complete("inversion partition validation", started);
  if (is_root) {
    std::printf("\nFunction PreCond_BICGSTAB:\n");
    std::printf("Using Fortran-style MPI Gatherv / rank-order sum / Bcast\n");
    std::printf("Fixed weights: continuous BiCGSTAB; robust weights: IRLS with Krylov restart\n");
    std::fflush(stdout);
  }
  for (std::size_t j = 0; j < columns; ++j) {
    DampingFactor1[j] = damping * Weight_Model[j] * Weight_Model[j];
    DampingFactor2[j] = DampingFactor1[j] * Prior_Model[j];
  }
  const auto dot = [](const std::vector<double> &a, const std::vector<double> &b) {
    long double result = 0.0L;
    for (std::size_t j = 0; j < a.size(); ++j)
      result += static_cast<long double>(a[j]) * b[j];
    return result;
  };
  const auto maximum = [](const std::vector<double> &a) {
    double result = 0.0;
    for (double value : a) {
      if (!std::isfinite(value)) throw std::runtime_error("Non-finite inversion vector.");
      result = std::max(result, std::abs(value));
    }
    return result;
  };
  const auto predict = [&]() {
    for (int i = 0; i < NumLocal; ++i) {
      double sum = 0.0;
      for (std::size_t j = 0; j < columns; ++j) sum += A[i][j] * m_predict[j];
      LocalPredict[i] = sum;
    }
  };
  const auto misfit = [&]() {
    MPIParallel::phase_begin("misfit gather and broadcast");
    const double begin = MPIParallel::time();
    double squared = 0.0;
    for (int i = 0; i < NumLocal; ++i) {
      const double error = d[i] - LocalPredict[i];
      squared += error * error;
    }
    MPIParallel::allreduce_sum(squared);
    const double result = std::sqrt(squared / NumData);
    if (!std::isfinite(result)) throw std::runtime_error("Non-finite data misfit.");
    MPIParallel::phase_complete("misfit gather and broadcast", begin);
    return result;
  };
  const auto update_weights = [&]() {
    for (int i = 0; i < NumLocal; ++i) {
      const double error = std::abs(d[i] - LocalPredict[i]);
      const double cutoff = RobustConstant * Sigma[i];
      double factor = 1.0;
      if (RobustSwitch == 1 && error > cutoff) factor = std::sqrt(cutoff / error);
      if (RobustSwitch == 2)
        factor = error < cutoff ? 1.0 - std::pow(error / cutoff, 2.0) : 0.0;
      ReW[i] = Weight_Data[i] * factor / Sigma[i];
      if (!std::isfinite(ReW[i])) throw std::runtime_error("Non-finite robust data weight.");
    }
  };
  // Apply the symmetric preconditioned normal operator without forming a normal matrix.
  const auto normal = [&](const std::vector<double> &x, std::vector<double> &result) {
    MPIParallel::phase_begin("normal operator local products");
    double begin = MPIParallel::time();
    for (int i = 0; i < NumLocal; ++i) {
      double sum = 0.0;
      for (std::size_t j = 0; j < columns; ++j) sum += A[i][j] * F[j] * x[j];
      TempVecRows[i] = ReW[i] * sum;
    }
    std::fill(TempVecCols.begin(), TempVecCols.end(), 0.0);
    for (int i = 0; i < NumLocal; ++i) {
      const double multiplier = ReW[i] * TempVecRows[i];
      const double *row = A[i];
      for (std::size_t j = 0; j < columns; ++j) TempVecCols[j] += row[j] * multiplier;
    }
    MPIParallel::phase_complete("normal operator local products", begin);
    MPIParallel::phase_begin("normal operator coefficient gather and broadcast");
    begin = MPIParallel::time();
    MPIParallel::allreduce_sum(TempVecCols);
    for (std::size_t j = 0; j < columns; ++j)
      result[j] = F[j] * (TempVecCols[j] + DampingFactor1[j] * F[j] * x[j]);
    maximum(result);
    MPIParallel::phase_complete("normal operator coefficient gather and broadcast", begin);
  };

  MPIParallel::phase_begin("initial model local prediction");
  started = MPIParallel::time();
  predict();
  MPIParallel::phase_complete("initial model local prediction", started);
  const double initial_misfit = misfit();
  Misfit[0] = initial_misfit;
  m = m_predict;
  if (is_root) {
    std::printf("The Misfit of The Initial Model is %.17g\n", initial_misfit);
    std::fflush(stdout);
  }
  if (damping == 0.0 && initial_misfit < epsilon) {
    prediction_gather.collect(LocalPredict, d_Predict);
    if (is_root) { std::puts("Initial model meets data fitting tolerance."); std::fflush(stdout); }
    return SolveStatus::data_fit;
  }
  update_weights();
  SolveStatus termination = SolveStatus::iteration_limit;
  double residual_scale = 1.0, reference_scale = 0.0, residual_tolerance = normal_tolerance;
  bool restart = true;
  for (int iteration = 1; iteration <= Max_Iteration; ++iteration) {
    if (is_root) {
      std::printf("Iteration Progress Activated, No. %d\n", iteration);
      std::fflush(stdout);
    }
    const bool rebuild = iteration == 1 || RobustSwitch != 0;
    if (rebuild) {
      MPIParallel::phase_begin("preconditioner and RHS local products");
      started = MPIParallel::time();
      std::fill(F.begin(), F.end(), 0.0);
      std::fill(newd.begin(), newd.end(), 0.0);
      for (int i = 0; i < NumLocal; ++i) {
        const double *row = A[i];
        const double weight2 = ReW[i] * ReW[i], weighted_data = weight2 * d[i];
        for (std::size_t j = 0; j < columns; ++j) {
          F[j] += row[j] * row[j] * weight2;
          newd[j] += row[j] * weighted_data;
        }
      }
      MPIParallel::phase_complete("preconditioner and RHS local products", started);
      MPIParallel::phase_begin("preconditioner and RHS gather and broadcast");
      started = MPIParallel::time();
      MPIParallel::allreduce_sum(F);
      MPIParallel::allreduce_sum(newd);
      for (std::size_t j = 0; j < columns; ++j) {
        const double diagonal = F[j] + DampingFactor1[j];
        if (!std::isfinite(diagonal) || !std::isfinite(newd[j]))
          throw std::runtime_error("Non-finite preconditioned system.");
        if (diagonal <= 0.0) { F[j] = 0.0; newd[j] = 0.0; }
        else { F[j] = 1.0 / std::sqrt(diagonal); newd[j] = (newd[j]+DampingFactor2[j])*F[j]; }
      }
      MPIParallel::phase_complete("preconditioner and RHS gather and broadcast", started);
      restart = true;
    }
    if (restart) {
      for (std::size_t j = 0; j < columns; ++j)
        m_old[j] = F[j] != 0.0 ? m_predict[j] / F[j] : m_predict[j];
      normal(m_old, Apj);
      reference_scale = 0.0;
      for (std::size_t j = 0; j < columns; ++j) {
        r_new[j] = newd[j] - Apj[j];
        reference_scale = std::max(reference_scale, std::max(std::abs(newd[j]), std::abs(Apj[j])));
      }
      residual_scale = maximum(r_new);
      if (residual_scale <= normal_tolerance * reference_scale) {
        termination = SolveStatus::converged;
        break;
      }
      residual_tolerance = normal_tolerance * (reference_scale / residual_scale);
      for (std::size_t j = 0; j < columns; ++j) pj[j] = r0[j] = r_new[j] /= residual_scale;
      restart = false;
    }
    normal(pj, Apj);
    const long double rho = dot(r_new, r0), denominator = dot(Apj, r0);
    const long double rho_bound = std::sqrt(dot(r_new, r_new)) * std::sqrt(dot(r0, r0));
    const long double alpha_bound = std::sqrt(dot(Apj, Apj)) * std::sqrt(dot(r0, r0));
    if (!std::isfinite(rho) || !std::isfinite(denominator) ||
        std::abs(rho) <= normal_tolerance * rho_bound ||
        std::abs(denominator) <= normal_tolerance * alpha_bound)
      throw std::runtime_error("BiCGSTAB alpha or shadow-residual breakdown.");
    const double alpha = static_cast<double>(rho / denominator);
    if (!std::isfinite(alpha)) throw std::runtime_error("Non-finite BiCGSTAB alpha.");
    for (std::size_t j = 0; j < columns; ++j) sj[j] = r_new[j] - alpha * Apj[j];
    const bool short_step = maximum(sj) <= residual_tolerance;
    double omega = 0.0;
    if (!short_step) {
      normal(sj, Asj);
      const long double norm2 = dot(Asj, Asj);
      if (!(norm2 > 0.0L) || !std::isfinite(norm2))
        throw std::runtime_error("BiCGSTAB omega breakdown.");
      omega = static_cast<double>(dot(Asj, sj) / norm2);
      if (!std::isfinite(omega) || omega == 0.0)
        throw std::runtime_error("Non-finite or zero BiCGSTAB omega.");
    }
    for (std::size_t j = 0; j < columns; ++j) {
      m_old[j] += residual_scale * (alpha * pj[j] + omega * sj[j]);
      if (F[j] != 0.0) m_predict[j] = F[j] * m_old[j];
      m[j] = m_predict[j];
      r_new[j] = short_step ? sj[j] : sj[j] - omega * Asj[j];
    }
    maximum(m_predict);
    MPIParallel::phase_begin("updated model local prediction");
    started = MPIParallel::time();
    predict();
    MPIParallel::phase_complete("updated model local prediction", started);
    Misfit[iteration-1] = misfit();
    if (is_root) {
      std::printf("The Misfit of Current Iteration is: %.17g\n", Misfit[iteration-1]);
      std::fflush(stdout);
    }
    if (damping == 0.0 && Misfit[iteration-1] < epsilon) {
      termination = SolveStatus::data_fit;
      break;
    }
    if (RobustSwitch != 0) {
      // Changing weights changes the normal operator, so restart Krylov after each IRLS update.
      update_weights();
      restart = true;
      continue;
    }
    if (short_step || maximum(r_new) <= residual_tolerance) {
      // Verify a recursively small residual against the actual weighted normal equations.
      normal(m_old, Apj);
      for (std::size_t j = 0; j < columns; ++j) TempVecCols[j] = newd[j] - Apj[j];
      if (maximum(TempVecCols) <= normal_tolerance * std::max(reference_scale, maximum(Apj))) {
        termination = SolveStatus::converged;
        break;
      }
      restart = true;
      continue;
    }
    const double beta = static_cast<double>((dot(r_new, r0) / rho) * (alpha / omega));
    if (!std::isfinite(beta)) throw std::runtime_error("Non-finite BiCGSTAB beta.");
    for (std::size_t j = 0; j < columns; ++j)
      pj[j] = r_new[j] + beta * (pj[j] - omega * Apj[j]);
  }
  MPIParallel::phase_begin("final prediction gather and broadcast");
  started = MPIParallel::time();
  prediction_gather.collect(LocalPredict, d_Predict);
  MPIParallel::phase_complete("final prediction gather and broadcast", started);
  if (is_root) {
    std::printf("Termination: %s\n",
        termination == SolveStatus::converged ? "normal equations converged" :
        termination == SolveStatus::data_fit ? "data fitting tolerance" : "iteration limit");
    std::puts("Function PrecondBiCGSTAB: Step out!\n");
    std::fflush(stdout);
  }
  return termination;
}
}
