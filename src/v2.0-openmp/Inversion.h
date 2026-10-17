#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <omp.h>
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
private:
  int columns_ = 0;
  std::vector<double> values_;
};

inline SolveStatus PreCond_BICGSTAB(
    int NumData, int NumCoefficients, const DenseMatrix &A,
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
      A.rows() != NumData || d.size() != static_cast<std::size_t>(NumData))
    throw std::invalid_argument("Inconsistent inversion vector dimensions.");
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

  std::vector<double> TempVecRows(NumData, 0.0);
  std::vector<double> TempVecCols(NumCoefficients, 0.0);

  std::vector<double> DampingFactor1(NumCoefficients, 0.0);
  std::vector<double> DampingFactor2(NumCoefficients, 0.0);
  std::vector<double> F(NumCoefficients, 0.0);
  std::vector<double> newd(NumCoefficients, 0.0);

  std::vector<double> r0(NumCoefficients, 0.0);
  std::vector<double> pj(NumCoefficients, 0.0);
  std::vector<double> m_old(NumCoefficients, 0.0);
  std::vector<double> m_new(NumCoefficients, 0.0);
  std::vector<double> Apj(NumCoefficients, 0.0);
  std::vector<double> sj(NumCoefficients, 0.0);
  std::vector<double> Asj(NumCoefficients, 0.0);

  std::vector<double> ReW(NumData, 0.0);

  std::printf("\nFunction PreCond_BICGSTAB:\n");
  std::printf(
      "Using Preconditioned BiCGSTAB (IRLS outer, 1 step per outer iter)\n");

#pragma omp parallel for if (NumCoefficients > 256)
  for (int j = 0; j < NumCoefficients; ++j) {
    DampingFactor1[j] = damping * Weight_Model[j] * Weight_Model[j];
    DampingFactor2[j] = DampingFactor1[j] * Prior_Model[j];
  }

#pragma omp parallel for if (NumData > 64)
  for (int i = 0; i < NumData; ++i) {
    double s = 0.0;
    for (int j = 0; j < NumCoefficients; ++j)
      s += A[i][j] * m_predict[j];
    d_Predict[i] = s;
  }

  double mis0 = 0.0;
#pragma omp parallel for reduction(+ : mis0) if (NumData > 256)
  for (int i = 0; i < NumData; ++i) {
    const double e = d[i] - d_Predict[i];
    mis0 += e * e;
  }
  Misfit[0] = std::sqrt(mis0 / std::max(1, NumData));

  if (damping == 0.0 && Misfit[0] < epsilon) {
    std::puts("\nThe Initial Model Meets the Fitting Requirements");
    std::puts("Iterations Ended without Beginning");
    std::printf(
        "The Iteration Terminated When the Fitting Difference is %f\n\n",
        Misfit[0]);
    m = m_predict;
    std::puts("Function PrecondBiCGSTAB: Step out!\n");
    return SolveStatus::data_fit;
  }
  std::printf("The Misfit of The Initial Model is %lf\n\n", Misfit[0]);

  if (RobustSwitch == 1) {
#pragma omp parallel for if (NumData > 64)
    for (int i = 0; i < NumData; ++i) {
      const double absres = std::abs(d[i] - d_Predict[i]);
      if (absres > RobustConstant * Sigma[i]) {
        ReW[i] = Weight_Data[i] * (1.0 / Sigma[i]) *
                 std::sqrt(RobustConstant * Sigma[i] / absres);
      } else {
        ReW[i] = Weight_Data[i] * (1.0 / Sigma[i]);
      }
    }
  } else if (RobustSwitch == 2) {
#pragma omp parallel for if (NumData > 64)
    for (int i = 0; i < NumData; ++i) {
      const double absres = std::abs(d[i] - d_Predict[i]);
      if (absres < RobustConstant * Sigma[i]) {
        ReW[i] = Weight_Data[i] * (1.0 / Sigma[i]) *
                 (1.0 - pow(absres / (RobustConstant * Sigma[i]), 2.0));
      } else {
        ReW[i] = 0.0;
      }
    }
  } else {
#pragma omp parallel for if (NumData > 64)
    for (int i = 0; i < NumData; ++i) {
      ReW[i] = Weight_Data[i] / Sigma[i];
    }
  }

  SolveStatus termination = SolveStatus::iteration_limit;
  int Num_Iteration = 0;
  for (Num_Iteration = 1; Num_Iteration <= Max_Iteration; ++Num_Iteration) {
    std::printf("Iteration Progress Activated, No. %d\n", Num_Iteration);

#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j) {
      double Fj = 0.0, nd = 0.0;
      for (int i = 0; i < NumData; ++i) {
        const double aij = A[i][j];
        const double rw2 = ReW[i] * ReW[i];
        Fj += aij * aij * rw2;
        nd += aij * rw2 * d[i];
      }
      Fj += DampingFactor1[j];
      if (Fj <= 0.0) {
        F[j] = 0.0;
        newd[j] = 0.0;
      } else {
        F[j] = 1.0 / std::sqrt(Fj);
        newd[j] = (nd + DampingFactor2[j]) * F[j];
      }
    }

#pragma omp parallel for if (NumData > 64)
    for (int i = 0; i < NumData; ++i) {
      double s = 0.0;
      const double w = ReW[i];
      for (int j = 0; j < NumCoefficients; ++j)
        s += A[i][j] * w * m_predict[j];
      TempVecRows[i] = s;
    }
#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j) {
      double s = 0.0;
      for (int i = 0; i < NumData; ++i)
        s += A[i][j] * ReW[i] * TempVecRows[i];
      TempVecCols[j] = s;
    }
#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j) {
      const double col = TempVecCols[j] + DampingFactor1[j] * m_predict[j];
      r0[j] = newd[j] - F[j] * col;
      pj[j] = r0[j];
      m_old[j] = (F[j] != 0.0) ? (m_predict[j] / F[j]) : m_predict[j];
    }

    double residual_scale = 0.0, reference_scale = 0.0;
    for (int j = 0; j < NumCoefficients; ++j) {
      if (!std::isfinite(r0[j]) || !std::isfinite(newd[j]) || !std::isfinite(F[j]))
        throw std::runtime_error("Non-finite preconditioned system.");
      residual_scale = std::max(residual_scale, std::abs(r0[j]));
      reference_scale = std::max(reference_scale,
          std::max(std::abs(newd[j]), std::abs(newd[j]-r0[j])));
    }
    if (residual_scale <= normal_tolerance * reference_scale) {
      m = m_predict;
      termination = SolveStatus::converged;
      break;
    }
    for (int j = 0; j < NumCoefficients; ++j)
      pj[j] = r0[j] = r0[j] / residual_scale;

#pragma omp parallel for if (NumData > 64)
    for (int i = 0; i < NumData; ++i) {
      double s = 0.0;
      const double w = ReW[i];
      for (int j = 0; j < NumCoefficients; ++j)
        s += A[i][j] * F[j] * pj[j];
      TempVecRows[i] = w * s;
    }
#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j) {
      double s = 0.0;
      for (int i = 0; i < NumData; ++i)
        s += A[i][j] * ReW[i] * TempVecRows[i];
      TempVecCols[j] = s;
    }
#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j) {
      Apj[j] = F[j] * (TempVecCols[j] + DampingFactor1[j] * F[j] * pj[j]);
    }

    double dot_r0_Apj = 0.0, dot_r0_r0 = 0.0, dot_Apj_Apj = 0.0;
#pragma omp parallel for reduction(+ : dot_r0_Apj,                             \
                                       dot_r0_r0, dot_Apj_Apj) if (NumCoefficients > 256)
    for (int j = 0; j < NumCoefficients; ++j) {
      dot_r0_Apj += Apj[j] * r0[j];
      dot_r0_r0 += r0[j] * r0[j];
      dot_Apj_Apj += Apj[j] * Apj[j];
    }
    if (!std::isfinite(dot_r0_Apj) || !(dot_r0_Apj >
        normal_tolerance * std::sqrt(dot_r0_r0) * std::sqrt(dot_Apj_Apj)))
      throw std::runtime_error("BiCGSTAB alpha breakdown.");
    const double aj = dot_r0_r0 / dot_r0_Apj;

#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j)
      sj[j] = r0[j] - aj * Apj[j];

#pragma omp parallel for if (NumData > 64)
    for (int i = 0; i < NumData; ++i) {
      double s = 0.0;
      const double w = ReW[i];
      for (int j = 0; j < NumCoefficients; ++j)
        s += A[i][j] * F[j] * sj[j];
      TempVecRows[i] = w * s;
    }
#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j) {
      double s = 0.0;
      for (int i = 0; i < NumData; ++i)
        s += A[i][j] * ReW[i] * TempVecRows[i];
      TempVecCols[j] = s;
    }
#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j) {
      Asj[j] = F[j] * (TempVecCols[j] + DampingFactor1[j] * F[j] * sj[j]);
    }

    double s_scale = 0.0, as_scale = 0.0;
    for (int j = 0; j < NumCoefficients; ++j) {
      s_scale = std::max(s_scale, std::abs(sj[j]));
      as_scale = std::max(as_scale, std::abs(Asj[j]));
    }
    double wj = 0.0;
    if (s_scale != 0.0) {
      if (!(as_scale > 0.0) || !std::isfinite(as_scale))
        throw std::runtime_error("BiCGSTAB omega breakdown.");
      double numerator = 0.0, denominator = 0.0;
      for (int j = 0; j < NumCoefficients; ++j) {
        const double a = Asj[j] / as_scale, b = sj[j] / s_scale;
        numerator += a * b;
        denominator += a * a;
      }
      wj = (numerator / denominator) * (s_scale / as_scale);
      if (!std::isfinite(wj))
        throw std::runtime_error("Non-finite BiCGSTAB step.");
    }

#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j) {
      m_new[j] = m_old[j] + residual_scale * (aj * pj[j] + wj * sj[j]);
    }

#pragma omp parallel for if (NumCoefficients > 64)
    for (int j = 0; j < NumCoefficients; ++j) {
      m_predict[j] = F[j] != 0.0 ? F[j] * m_new[j] : m_predict[j];
      m[j] = m_predict[j];
    }

#pragma omp parallel for if (NumData > 64)
    for (int i = 0; i < NumData; ++i) {
      double s = 0.0;
      for (int j = 0; j < NumCoefficients; ++j)
        s += A[i][j] * m_predict[j];
      d_Predict[i] = s;
    }

    double cur = 0.0;
#pragma omp parallel for reduction(+ : cur) if (NumData > 256)
    for (int i = 0; i < NumData; ++i) {
      const double e = d[i] - d_Predict[i];
      cur += e * e;
    }
    Misfit[Num_Iteration - 1] = std::sqrt(cur / std::max(1, NumData));
    if (!std::isfinite(Misfit[Num_Iteration - 1]) ||
        !std::all_of(m_predict.begin(), m_predict.end(), [](double x){return std::isfinite(x);}) ||
        !std::all_of(d_Predict.begin(), d_Predict.end(), [](double x){return std::isfinite(x);}))
      throw std::runtime_error("Non-finite inversion iterate.");

    std::printf("The Misfit of Current Iteration is: %f\n\n",
                Misfit[Num_Iteration - 1]);
    if (damping == 0.0 && Misfit[Num_Iteration - 1] < epsilon) {
      termination = SolveStatus::data_fit;
      break;
    }

    if (RobustSwitch == 1) {
#pragma omp parallel for if (NumData > 64)
      for (int i = 0; i < NumData; ++i) {
        const double absres = std::abs(d[i] - d_Predict[i]);
        if (absres > RobustConstant * Sigma[i]) {
          ReW[i] = Weight_Data[i] * (1.0 / Sigma[i]) *
                   std::sqrt(RobustConstant * Sigma[i] / absres);
        } else {
          ReW[i] = Weight_Data[i] * (1.0 / Sigma[i]);
        }
      }
    } else if (RobustSwitch == 2) {
#pragma omp parallel for if (NumData > 64)
      for (int i = 0; i < NumData; ++i) {
        const double absres = std::abs(d[i] - d_Predict[i]);
        if (absres < RobustConstant * Sigma[i]) {
          ReW[i] = Weight_Data[i] * (1.0 / Sigma[i]) *
                   (1.0 - pow(absres / (RobustConstant * Sigma[i]), 2.0));
        } else {
          ReW[i] = 0.0;
        }
      }
    }
  }

  std::puts("Function PrecondBiCGSTAB: Step out!\n");
  std::printf("Termination: %s\n",
      termination == SolveStatus::converged ? "normal equations converged" :
      termination == SolveStatus::data_fit ? "data fitting tolerance" : "iteration limit");
  return termination;
}

}
