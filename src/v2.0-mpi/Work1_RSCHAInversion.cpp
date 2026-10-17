#include "Basis Function.h"
#include "RSCHAConsole.h"
#include "CoordinateTransformation.h"
#include "Inversion.h"
#include "MPIParallel.h"
#include "R_SCHA.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace CoordinateTransformation;
using namespace R_SCHA;
using namespace Inversion;

namespace {
std::string input_failure_detail(int code, bool loading) {
  const char *where = "modelling parameters or data dimensions";
  if (loading) {
    switch (code) {
    case -1: where = "Parameters_Modelling.dat"; break;
    case -2: where = "Parameters_Inversion.dat"; break;
    case -3: where = "Data_TMI.dat (expected 6 columns)"; break;
    case -4: where = "Parameters_DataWeightTMI.dat or Parameters_SigmaTMI.dat (expected 1 column)"; break;
    case -5: where = "Data_Vector.dat (expected 6 columns); check IncludeVector"; break;
    case -6: where = "Parameters_DataWeightVector.dat or Parameters_SigmaVector.dat (expected 3 columns)"; break;
    case -7: where = "satellite endpoint-difference data files (expected 12 columns)"; break;
    case -8: where = "satellite endpoint-difference weight or Sigma files"; break;
    case -9: where = "Parameters_InitialModel.dat"; break;
    case -10: where = "Parameters_ModelWeight.dat or Parameters_PriorModel.dat"; break;
    }
  }
  return std::string(loading ? "Input loading failed" : "Input validation failed") +
      " (code " + std::to_string(code) + "): " + where;
}
}

int main(int argc, char **argv) {

  MPIParallel::Environment mpi(argc, argv);
  try {
  if (!mpi.is_root())
    std::cout.setstate(std::ios_base::failbit);

  if (argc == 2 && std::string(argv[1]) == "--mpi-check") {
    MPIParallel::communication_check();
    return 0;
  }
  RSCHAConsole::banner("R-SCHA model inversion", "MPI");
  RSCHAConsole::section(1, "Input files and validation");

  const std::string ModellingParametersFile = "Parameters_Modelling.dat";
  const std::string InversionParametersFile = "Parameters_Inversion.dat";

  const bool IncludeTMI = true;
  const bool IncludeVector = false;
  const bool IncludeDiffTMI = true;
  const bool IncludeDiffVector = true;

  const std::string DataFile_TMI = "Data_TMI.dat";
  const std::string DataFile_Vector = "Data_Vector.dat";
  const std::string DataFile_DiffTMI = "Data_ALongTrackDiffTMI.dat";
  const std::string DataFile_DiffVector = "Data_ALongTrackDiffVector.dat";

  const std::string DataWeightFile_TMI = "Parameters_DataWeightTMI.dat";
  const std::string DataWeightFile_Vector = "Parameters_DataWeightVector.dat";
  const std::string DataWeightFile_DiffTMI =
      "Parameters_DataWeightALongTrackDiffTMI.dat";
  const std::string DataWeightFile_DiffVector =
      "Parameters_DataWeightALongTrackDiffVector.dat";

  const std::string SigmaFile_TMI = "Parameters_SigmaTMI.dat";
  const std::string SigmaFile_Vector = "Parameters_SigmaVector.dat";
  const std::string SigmaFile_DiffTMI = "Parameters_SigmaALongTrackDiffTMI.dat";
  const std::string SigmaFile_DiffVector =
      "Parameters_SigmaALongTrackDiffVector.dat";

  const std::string InitialModelFile = "Parameters_InitialModel.dat";
  const std::string ModelWeightFile = "Parameters_ModelWeight.dat";
  const std::string PriorModelFile = "Parameters_PriorModel.dat";

  const std::string ReconsTMIFile = "Recons_TMIData.dat";
  const std::string ReconsVectorFile = "Recons_VectorData.dat";
  const std::string ReconsDiffTMIFile = "Recons_ALongTrackDiffTMIData.dat";
  const std::string ReconsDiffVectorFile =
      "Recons_ALongTrackDiffVectorData.dat";

  const std::string ModellingResultsFile = "Results_RSCHAModel.dat";
  const std::string ModelSpectrumFile = "Results_RSCHASpectrum.dat";

  double Lon0, Lat0, Theta0;
  int K_Max_Internal, K_Max_External, P_Max, M_Max;
  double Max_Radius, Min_Radius;

  double epsilon, damping;
  int Max_Iteration;
  int RobustSwitch;
  double RobustConstant;

  std::vector<double> LongitudeTMI, LatitudeTMI, RadiusTMI, IncTMI, DecTMI, TMI;

  std::vector<double> LongitudeVector, LatitudeVector, RadiusVector, Bx, By, Bz;

  std::vector<double> LongitudeDiffTMI1, LatitudeDiffTMI1, RadiusDiffTMI1,
      IncDiffTMI1, DecDiffTMI1, DiffTMI1;
  std::vector<double> LongitudeDiffTMI2, LatitudeDiffTMI2, RadiusDiffTMI2,
      IncDiffTMI2, DecDiffTMI2, DiffTMI2;

  std::vector<double> LongitudeDiffVector1, LatitudeDiffVector1,
      RadiusDiffVector1, DiffBx1, DiffBy1, DiffBz1;
  std::vector<double> LongitudeDiffVector2, LatitudeDiffVector2,
      RadiusDiffVector2, DiffBx2, DiffBy2, DiffBz2;

  std::vector<double> DataWeight_TMI, DataWeight_Bx, DataWeight_By,
      DataWeight_Bz, DataWeight_DiffTMI, DataWeight_DiffBx, DataWeight_DiffBy,
      DataWeight_DiffBz;

  std::vector<double> Sigma_TMI, Sigma_Bx, Sigma_By, Sigma_Bz, Sigma_DiffTMI,
      Sigma_DiffBx, Sigma_DiffBy, Sigma_DiffBz;

  std::vector<double> m0;

  std::vector<double> Weight_Model;

  std::vector<double> Prior_Model;

  int Flag = 1;

  {
    Flag = ReadModellingParameters(
        ModellingParametersFile, InversionParametersFile, InitialModelFile,
        ModelWeightFile, PriorModelFile, IncludeTMI, IncludeVector,
        IncludeDiffTMI, IncludeDiffVector, DataFile_TMI, DataFile_Vector,
        DataFile_DiffTMI, DataFile_DiffVector, DataWeightFile_TMI,
        DataWeightFile_Vector, DataWeightFile_DiffTMI,
        DataWeightFile_DiffVector, SigmaFile_TMI, SigmaFile_Vector,
        SigmaFile_DiffTMI, SigmaFile_DiffVector, Lon0, Lat0, Theta0,
        K_Max_Internal, K_Max_External, P_Max, M_Max, Max_Radius, Min_Radius,
        epsilon, damping, Max_Iteration, RobustSwitch, RobustConstant,
        LongitudeTMI, LatitudeTMI, RadiusTMI, IncTMI, DecTMI, TMI,
        DataWeight_TMI, Sigma_TMI, LongitudeVector, LatitudeVector,
        RadiusVector, Bx, By, Bz, DataWeight_Bx, DataWeight_By, DataWeight_Bz,
        Sigma_Bx, Sigma_By, Sigma_Bz, LongitudeDiffTMI1, LatitudeDiffTMI1,
        RadiusDiffTMI1, IncDiffTMI1, DecDiffTMI1, DiffTMI1, LongitudeDiffTMI2,
        LatitudeDiffTMI2, RadiusDiffTMI2, IncDiffTMI2, DecDiffTMI2, DiffTMI2,
        DataWeight_DiffTMI, Sigma_DiffTMI, LongitudeDiffVector1,
        LatitudeDiffVector1, RadiusDiffVector1, DiffBx1, DiffBy1, DiffBz1,
        LongitudeDiffVector2, LatitudeDiffVector2, RadiusDiffVector2, DiffBx2,
        DiffBy2, DiffBz2, DataWeight_DiffBx, DataWeight_DiffBy,
        DataWeight_DiffBz, Sigma_DiffBx, Sigma_DiffBy, Sigma_DiffBz, m0,
        Weight_Model, Prior_Model);

    if (Flag <= 0) {
      std::cout << "\nModelling parameter settings and data format are "
                   "incorrect, calculation failed, exiting ..."
                << std::endl;
      std::cout << "Error code: " << Flag << std::endl;
      throw std::runtime_error(input_failure_detail(Flag, true));
    }
  }

  {
    Flag = ParametersCheck(
        IncludeTMI, IncludeVector, IncludeDiffTMI, IncludeDiffVector, Lon0,
        Lat0, Theta0, K_Max_Internal, K_Max_External, P_Max, M_Max, Max_Radius,
        Min_Radius, epsilon, damping, Max_Iteration, RobustSwitch,
        RobustConstant, LongitudeTMI, LatitudeTMI, RadiusTMI, IncTMI, DecTMI,
        TMI, DataWeight_TMI, Sigma_TMI, LongitudeVector, LatitudeVector,
        RadiusVector, Bx, By, Bz, DataWeight_Bx, DataWeight_By, DataWeight_Bz,
        Sigma_Bx, Sigma_By, Sigma_Bz, LongitudeDiffTMI1, LatitudeDiffTMI1,
        RadiusDiffTMI1, IncDiffTMI1, DecDiffTMI1, DiffTMI1, LongitudeDiffTMI2,
        LatitudeDiffTMI2, RadiusDiffTMI2, IncDiffTMI2, DecDiffTMI2, DiffTMI2,
        DataWeight_DiffTMI, Sigma_DiffTMI, LongitudeDiffVector1,
        LatitudeDiffVector1, RadiusDiffVector1, DiffBx1, DiffBy1, DiffBz1,
        LongitudeDiffVector2, LatitudeDiffVector2, RadiusDiffVector2, DiffBx2,
        DiffBy2, DiffBz2, DataWeight_DiffBx, DataWeight_DiffBy,
        DataWeight_DiffBz, Sigma_DiffBx, Sigma_DiffBy, Sigma_DiffBz);

    if (Flag <= 0) {
      std::cout << "\nModelling parameter settings and data format checking "
                   "failed, calculation failed, exiting ..."
                << std::endl;
      std::cout << "Error code: " << Flag << std::endl;
      throw std::runtime_error(input_failure_detail(Flag, false));
    }
  }

  const int NumTMI = LongitudeTMI.size();
  const int NumVector = LongitudeVector.size();
  const int NumDiffTMI = LongitudeDiffTMI1.size();
  const int NumDiffVector = LongitudeDiffVector1.size();
  const int NumData = NumTMI + 3 * NumVector + NumDiffTMI + 3 * NumDiffVector;

  const int Num_CoefIn = (K_Max_Internal + 1) * (K_Max_Internal + 2) / 2;
  const int Num_CoefEx = (K_Max_External + 1) * (K_Max_External + 2) / 2;
  const int Num_CoefMe = P_Max * (P_Max + 1);
  const int Num_CoefCo = M_Max + 1;
  const int Num_Coef = coefficient_count(K_Max_Internal, K_Max_External, P_Max, M_Max);
  if (NumData <= 0 || m0.size() != static_cast<std::size_t>(Num_Coef) ||
      Weight_Model.size() != m0.size() || Prior_Model.size() != m0.size())
    throw std::invalid_argument("Empty data or inconsistent model vector dimensions.");

  // Return before constructing any design matrix or computing basis functions.
  if (argc == 2 && std::string(argv[1]) == "--validate-only") {
    std::cout << "Input validation passed: " << NumData << " active data rows, "
              << Num_Coef << " model coefficients. No inversion performed." << std::endl;
    return 0;
  }


  RSCHAConsole::section(2, "Model, data and solver settings");
  RSCHAConsole::model(Lon0, Lat0, Theta0, K_Max_Internal, K_Max_External,
                      P_Max, M_Max, Min_Radius, Max_Radius);
  RSCHAConsole::row("TMI records", NumTMI);
  RSCHAConsole::row("Vector records", NumVector);
  RSCHAConsole::row("Differenced TMI pairs", NumDiffTMI);
  RSCHAConsole::row("Differenced vector pairs", NumDiffVector);
  RSCHAConsole::row("Active scalar rows", NumData);
  RSCHAConsole::row("Model coefficients", Num_Coef);
  RSCHAConsole::row("Maximum iterations", Max_Iteration);
  RSCHAConsole::row("Fitting tolerance", epsilon);
  RSCHAConsole::row("Damping", damping);
  RSCHAConsole::row("Robust switch / constant", std::to_string(RobustSwitch)+" / "+std::to_string(RobustConstant));
  RSCHAConsole::row("OpenMP threads per process", omp_get_max_threads());
  RSCHAConsole::section(3, "Coordinates and regional basis");

  std::vector<double> LamdaTMI, ThetaTMI;
  std::vector<double> sinalpha_TMI, cosalpha_TMI;
  std::vector<double> LamdaVector, ThetaVector;
  std::vector<double> Bx_SCCS, By_SCCS, Bz_SCCS;
  std::vector<double> LamdaDiffTMI1, ThetaDiffTMI1;
  std::vector<double> LamdaDiffTMI2, ThetaDiffTMI2;
  std::vector<double> sinalpha_DiffTMI1, cosalpha_DiffTMI1;
  std::vector<double> sinalpha_DiffTMI2, cosalpha_DiffTMI2;
  std::vector<double> LamdaDiffVector1, ThetaDiffVector1;
  std::vector<double> LamdaDiffVector2, ThetaDiffVector2;
  std::vector<double> DiffBx1_SCCS, DiffBy1_SCCS, DiffBz1_SCCS;
  std::vector<double> DiffBx2_SCCS, DiffBy2_SCCS, DiffBz2_SCCS;

  if (IncludeTMI) {

    LamdaTMI.resize(NumTMI);
    ThetaTMI.resize(NumTMI);
    sinalpha_TMI.resize(NumTMI);
    cosalpha_TMI.resize(NumTMI);

    SCHA_CT(LongitudeTMI, LatitudeTMI, Lon0, Lat0, 1, LamdaTMI, ThetaTMI);

    for (int i = 0; i < NumTMI; ++i) {
      sinalpha_TMI[i] =
          sind(90.0 - Lat0) * sind(LongitudeTMI[i] - Lon0) / sind(ThetaTMI[i]);
      cosalpha_TMI[i] = (cosd(90.0 - Lat0) -
                         cosd(ThetaTMI[i]) * cosd(90.0 - LatitudeTMI[i])) /
                        (sind(ThetaTMI[i]) * sind(90.0 - LatitudeTMI[i]));
    }
  }

  if (IncludeVector) {

    Bx_SCCS.resize(NumVector);
    By_SCCS.resize(NumVector);
    Bz_SCCS.resize(NumVector);
    LamdaVector.resize(NumVector);
    ThetaVector.resize(NumVector);

    SCHA_CT(LongitudeVector, LatitudeVector, Lon0, Lat0, 1, LamdaVector,
            ThetaVector);
    SCHA_3CompsTrans(LongitudeVector, LatitudeVector, Lon0, Lat0, Bx, By, Bz, 1,
                     Bx_SCCS, By_SCCS, Bz_SCCS);
  }

  if (IncludeDiffTMI) {

    LamdaDiffTMI1.resize(NumDiffTMI);
    ThetaDiffTMI1.resize(NumDiffTMI);
    LamdaDiffTMI2.resize(NumDiffTMI);
    ThetaDiffTMI2.resize(NumDiffTMI);
    sinalpha_DiffTMI1.resize(NumDiffTMI);
    cosalpha_DiffTMI1.resize(NumDiffTMI);
    sinalpha_DiffTMI2.resize(NumDiffTMI);
    cosalpha_DiffTMI2.resize(NumDiffTMI);
    SCHA_CT(LongitudeDiffTMI1, LatitudeDiffTMI1, Lon0, Lat0, 1, LamdaDiffTMI1,
            ThetaDiffTMI1);
    SCHA_CT(LongitudeDiffTMI2, LatitudeDiffTMI2, Lon0, Lat0, 1, LamdaDiffTMI2,
            ThetaDiffTMI2);

    for (int i = 0; i < NumDiffTMI; ++i) {
      sinalpha_DiffTMI1[i] = sind(90.0 - Lat0) *
                             sind(LongitudeDiffTMI1[i] - Lon0) /
                             sind(ThetaDiffTMI1[i]);
      cosalpha_DiffTMI1[i] =
          (cosd(90.0 - Lat0) -
           cosd(ThetaDiffTMI1[i]) * cosd(90.0 - LatitudeDiffTMI1[i])) /
          (sind(ThetaDiffTMI1[i]) * sind(90.0 - LatitudeDiffTMI1[i]));
    }
    for (int i = 0; i < NumDiffTMI; ++i) {
      sinalpha_DiffTMI2[i] = sind(90.0 - Lat0) *
                             sind(LongitudeDiffTMI2[i] - Lon0) /
                             sind(ThetaDiffTMI2[i]);
      cosalpha_DiffTMI2[i] =
          (cosd(90.0 - Lat0) -
           cosd(ThetaDiffTMI2[i]) * cosd(90.0 - LatitudeDiffTMI2[i])) /
          (sind(ThetaDiffTMI2[i]) * sind(90.0 - LatitudeDiffTMI2[i]));
    }
  }

  if (IncludeDiffVector) {

    LamdaDiffVector1.resize(NumDiffVector);
    ThetaDiffVector1.resize(NumDiffVector);
    LamdaDiffVector2.resize(NumDiffVector);
    ThetaDiffVector2.resize(NumDiffVector);
    DiffBx1_SCCS.resize(NumDiffVector);
    DiffBy1_SCCS.resize(NumDiffVector);
    DiffBz1_SCCS.resize(NumDiffVector);
    DiffBx2_SCCS.resize(NumDiffVector);
    DiffBy2_SCCS.resize(NumDiffVector);
    DiffBz2_SCCS.resize(NumDiffVector);

    SCHA_CT(LongitudeDiffVector1, LatitudeDiffVector1, Lon0, Lat0, 1,
            LamdaDiffVector1, ThetaDiffVector1);
    SCHA_CT(LongitudeDiffVector2, LatitudeDiffVector2, Lon0, Lat0, 1,
            LamdaDiffVector2, ThetaDiffVector2);
    SCHA_3CompsTrans(LongitudeDiffVector1, LatitudeDiffVector1, Lon0, Lat0,
                     DiffBx1, DiffBy1, DiffBz1, 1, DiffBx1_SCCS, DiffBy1_SCCS,
                     DiffBz1_SCCS);
    SCHA_3CompsTrans(LongitudeDiffVector2, LatitudeDiffVector2, Lon0, Lat0,
                     DiffBx2, DiffBy2, DiffBz2, 1, DiffBx2_SCCS, DiffBy2_SCCS,
                     DiffBz2_SCCS);
  }

  int Count;
  std::vector<double> m_Legendre, nk;
  if (K_Max_Internal >= K_Max_External) {

    m_Legendre.resize(Num_CoefIn);
    nk.resize(Num_CoefIn);

    Count = 0;
    for (int i = 0; i <= K_Max_Internal; ++i) {
      for (int j = 0; j <= i; ++j) {
        m_Legendre[Count] = j;
        Count++;
      }
    }

    R_SCHA_nk(K_Max_Internal, Theta0, nk);
    std::cout << "\n";
  } else {

    m_Legendre.resize(Num_CoefEx);
    nk.resize(Num_CoefEx);

    Count = 0;
    for (int i = 0; i <= K_Max_External; ++i) {
      for (int j = 0; j <= i; ++j) {
        m_Legendre[Count] = j;
        Count++;
      }
    }

    R_SCHA_nk(K_Max_External, Theta0, nk);
    std::cout << "\n";
  }

  for (size_t i = 0; i < nk.size(); ++i)
    Legendre::PrepareLegacySchmidtMode(m_Legendre[i], nk[i]);

  std::vector<double> p(Num_CoefMe), m_Mehler(Num_CoefMe),
      m_Constant(Num_CoefCo);
  Count = 0;
  for (int i = 1; i <= P_Max; i++) {
    for (int j = 0; j <= P_Max; j++) {
      p[Count] = i;
      m_Mehler[Count] = j;
      Count++;
    }
  }

  for (int i = 0; i <= M_Max; i++) {
    m_Constant[i] = (double)i;
  }

  // Match the Fortran architecture: contiguous ownership within each data family.
  const std::vector<int> OwnedTMI = MPIParallel::owned_indices(NumTMI, mpi.rank(), mpi.size());
  const std::vector<int> OwnedVector = MPIParallel::owned_indices(NumVector, mpi.rank(), mpi.size());
  const std::vector<int> OwnedDiffTMI = MPIParallel::owned_indices(NumDiffTMI, mpi.rank(), mpi.size());
  const std::vector<int> OwnedDiffVector = MPIParallel::owned_indices(NumDiffVector, mpi.rank(), mpi.size());

  const int LocalTMIBase = 0;
  const int LocalVectorBase = static_cast<int>(OwnedTMI.size());
  const int LocalDiffTMIBase =
      LocalVectorBase + 3 * static_cast<int>(OwnedVector.size());
  const int LocalDiffVectorBase =
      LocalDiffTMIBase + static_cast<int>(OwnedDiffTMI.size());
  const int NumLocal =
      LocalDiffVectorBase + 3 * static_cast<int>(OwnedDiffVector.size());

  std::vector<int> GlobalRows;
  GlobalRows.reserve(NumLocal);
  for (int i : OwnedTMI)
    GlobalRows.push_back(i);
  for (int component = 0; component < 3; ++component)
    for (int i : OwnedVector)
      GlobalRows.push_back(NumTMI + component * NumVector + i);
  for (int i : OwnedDiffTMI)
    GlobalRows.push_back(NumTMI + 3 * NumVector + i);
  for (int component = 0; component < 3; ++component)
    for (int i : OwnedDiffVector)
      GlobalRows.push_back(NumTMI + 3 * NumVector + NumDiffTMI +
                           component * NumDiffVector + i);

  RSCHAConsole::section(4, "Kernel assembly");
  Inversion::DenseMatrix Matrix(NumLocal, Num_Coef);
  Inversion::DenseMatrix MatrixDiffTMI1(OwnedDiffTMI.size(), Num_Coef);
  Inversion::DenseMatrix MatrixDiffTMI2(OwnedDiffTMI.size(), Num_Coef);
  Inversion::DenseMatrix MatrixDiffVector1(3 * OwnedDiffVector.size(), Num_Coef);
  Inversion::DenseMatrix MatrixDiffVector2(3 * OwnedDiffVector.size(), Num_Coef);

  int dataidx, coefidx;
  double kernel_started = 0.0;
  std::vector<double> TempBx(Num_Coef), TempBy(Num_Coef);
  std::vector<double> KernelBx(Num_Coef), KernelBy(Num_Coef),
      KernelBz(Num_Coef);

  if (IncludeTMI) {
    MPIParallel::phase_begin("TMI kernels");
    kernel_started = MPIParallel::time();

    RSCHAConsole::kernel("TMI");

    for (int local = 0; local < static_cast<int>(OwnedTMI.size()); ++local) {
      if (mpi.is_root() && local % std::max(1, static_cast<int>(OwnedTMI.size())/20) == 0) {
        std::printf("[Kernel progress] TMI kernels: rank 0 local %d/%zu; elapsed=%.3f s\n",
            local, OwnedTMI.size(), MPIParallel::time()-kernel_started);
        std::fflush(stdout);
      }
      dataidx = OwnedTMI[local];

      RSCHA_VectorKernel(RadiusTMI[dataidx], ThetaTMI[dataidx], Theta0,
                         LamdaTMI[dataidx], Min_Radius, Max_Radius, Num_CoefIn,
                         Num_CoefEx, Num_CoefMe, Num_CoefCo, nk, m_Legendre, p,
                         m_Mehler, m_Constant, TempBx, TempBy, KernelBz);

      for (coefidx = 0; coefidx < Num_Coef; ++coefidx) {

        KernelBx[coefidx] = TempBx[coefidx] * cosalpha_TMI[dataidx] +
                            TempBy[coefidx] * sinalpha_TMI[dataidx];
        KernelBy[coefidx] = TempBy[coefidx] * cosalpha_TMI[dataidx] -
                            TempBx[coefidx] * sinalpha_TMI[dataidx];

        Matrix[LocalTMIBase + local][coefidx] =
            KernelBx[coefidx] * cosd(IncTMI[dataidx]) * cosd(DecTMI[dataidx]) +
            KernelBy[coefidx] * cosd(IncTMI[dataidx]) * sind(DecTMI[dataidx]) +
            KernelBz[coefidx] * sind(IncTMI[dataidx]);
      }
    }

    MPIParallel::phase_complete("TMI kernels", kernel_started);
    RSCHAConsole::done("Kernel block complete");
  }

  if (IncludeVector) {
    MPIParallel::phase_begin("vector kernels");
    kernel_started = MPIParallel::time();

    RSCHAConsole::kernel("Vector");

    for (int local = 0; local < static_cast<int>(OwnedVector.size()); ++local) {
      if (mpi.is_root() && local % std::max(1, static_cast<int>(OwnedVector.size())/20) == 0) {
        std::printf("[Kernel progress] vector kernels: rank 0 local %d/%zu; elapsed=%.3f s\n",
            local, OwnedVector.size(), MPIParallel::time()-kernel_started);
        std::fflush(stdout);
      }
      dataidx = OwnedVector[local];

      RSCHA_VectorKernel(RadiusVector[dataidx], ThetaVector[dataidx], Theta0,
                         LamdaVector[dataidx], Min_Radius, Max_Radius, Num_CoefIn,
                         Num_CoefEx, Num_CoefMe, Num_CoefCo, nk, m_Legendre, p,
                         m_Mehler, m_Constant, KernelBx, KernelBy, KernelBz);

      for (coefidx = 0; coefidx < Num_Coef; ++coefidx) {

        Matrix[LocalVectorBase + local][coefidx] = KernelBx[coefidx];
        Matrix[LocalVectorBase + local + OwnedVector.size()][coefidx] = KernelBy[coefidx];
        Matrix[LocalVectorBase + local + 2 * OwnedVector.size()][coefidx] = KernelBz[coefidx];
      }
    }

    MPIParallel::phase_complete("vector kernels", kernel_started);
    RSCHAConsole::done("Kernel block complete");
  }

  if (IncludeDiffTMI) {
    MPIParallel::phase_begin("differenced TMI kernels");
    kernel_started = MPIParallel::time();

    RSCHAConsole::kernel("Along-track differenced TMI");

    for (int local = 0; local < static_cast<int>(OwnedDiffTMI.size()); ++local) {
      if (mpi.is_root() && local % std::max(1, static_cast<int>(OwnedDiffTMI.size())/20) == 0) {
        std::printf("[Kernel progress] differenced TMI kernels: rank 0 local %d/%zu; elapsed=%.3f s\n",
            local, OwnedDiffTMI.size(), MPIParallel::time()-kernel_started);
        std::fflush(stdout);
      }
      dataidx = OwnedDiffTMI[local];

      RSCHA_VectorKernel(RadiusDiffTMI1[dataidx], ThetaDiffTMI1[dataidx],
                         Theta0, LamdaDiffTMI1[dataidx], Min_Radius, Max_Radius,
                         Num_CoefIn, Num_CoefEx, Num_CoefMe, Num_CoefCo, nk,
                         m_Legendre, p, m_Mehler, m_Constant, TempBx, TempBy,
                         KernelBz);

      for (coefidx = 0; coefidx < Num_Coef; ++coefidx) {

        KernelBx[coefidx] = TempBx[coefidx] * cosalpha_DiffTMI1[dataidx] +
                            TempBy[coefidx] * sinalpha_DiffTMI1[dataidx];
        KernelBy[coefidx] = TempBy[coefidx] * cosalpha_DiffTMI1[dataidx] -
                            TempBx[coefidx] * sinalpha_DiffTMI1[dataidx];

        MatrixDiffTMI1[local][coefidx] =
            KernelBx[coefidx] * cosd(IncDiffTMI1[dataidx]) *
                cosd(DecDiffTMI1[dataidx]) +
            KernelBy[coefidx] * cosd(IncDiffTMI1[dataidx]) *
                sind(DecDiffTMI1[dataidx]) +
            KernelBz[coefidx] * sind(IncDiffTMI1[dataidx]);
      }

      RSCHA_VectorKernel(RadiusDiffTMI2[dataidx], ThetaDiffTMI2[dataidx],
                         Theta0, LamdaDiffTMI2[dataidx], Min_Radius, Max_Radius,
                         Num_CoefIn, Num_CoefEx, Num_CoefMe, Num_CoefCo, nk,
                         m_Legendre, p, m_Mehler, m_Constant, TempBx, TempBy,
                         KernelBz);

      for (coefidx = 0; coefidx < Num_Coef; ++coefidx) {

        KernelBx[coefidx] = TempBx[coefidx] * cosalpha_DiffTMI2[dataidx] +
                            TempBy[coefidx] * sinalpha_DiffTMI2[dataidx];
        KernelBy[coefidx] = TempBy[coefidx] * cosalpha_DiffTMI2[dataidx] -
                            TempBx[coefidx] * sinalpha_DiffTMI2[dataidx];

        MatrixDiffTMI2[local][coefidx] =
            KernelBx[coefidx] * cosd(IncDiffTMI2[dataidx]) *
                cosd(DecDiffTMI2[dataidx]) +
            KernelBy[coefidx] * cosd(IncDiffTMI2[dataidx]) *
                sind(DecDiffTMI2[dataidx]) +
            KernelBz[coefidx] * sind(IncDiffTMI2[dataidx]);

        Matrix[LocalDiffTMIBase + local][coefidx] =
            MatrixDiffTMI1[local][coefidx] - MatrixDiffTMI2[local][coefidx];
      }
    }

    MPIParallel::phase_complete("differenced TMI kernels", kernel_started);
    RSCHAConsole::done("Kernel block complete");
  }

  if (IncludeDiffVector) {
    MPIParallel::phase_begin("differenced vector kernels");
    kernel_started = MPIParallel::time();

    RSCHAConsole::kernel("Along-track differenced vector");

    for (int local = 0; local < static_cast<int>(OwnedDiffVector.size()); ++local) {
      if (mpi.is_root() && local % std::max(1, static_cast<int>(OwnedDiffVector.size())/20) == 0) {
        std::printf("[Kernel progress] differenced vector kernels: rank 0 local %d/%zu; elapsed=%.3f s\n",
            local, OwnedDiffVector.size(), MPIParallel::time()-kernel_started);
        std::fflush(stdout);
      }
      dataidx = OwnedDiffVector[local];

      RSCHA_VectorKernel(RadiusDiffVector1[dataidx], ThetaDiffVector1[dataidx],
                         Theta0, LamdaDiffVector1[dataidx], Min_Radius, Max_Radius,
                         Num_CoefIn, Num_CoefEx, Num_CoefMe, Num_CoefCo, nk,
                         m_Legendre, p, m_Mehler, m_Constant, KernelBx,
                         KernelBy, KernelBz);

      for (coefidx = 0; coefidx < Num_Coef; ++coefidx) {

        MatrixDiffVector1[local][coefidx] = KernelBx[coefidx];
        MatrixDiffVector1[local + OwnedDiffVector.size()][coefidx] = KernelBy[coefidx];
        MatrixDiffVector1[local + 2 * OwnedDiffVector.size()][coefidx] =
            KernelBz[coefidx];
      }

      RSCHA_VectorKernel(RadiusDiffVector2[dataidx], ThetaDiffVector2[dataidx],
                         Theta0, LamdaDiffVector2[dataidx], Min_Radius, Max_Radius,
                         Num_CoefIn, Num_CoefEx, Num_CoefMe, Num_CoefCo, nk,
                         m_Legendre, p, m_Mehler, m_Constant, KernelBx,
                         KernelBy, KernelBz);

      for (coefidx = 0; coefidx < Num_Coef; ++coefidx) {

        MatrixDiffVector2[local][coefidx] = KernelBx[coefidx];
        MatrixDiffVector2[local + OwnedDiffVector.size()][coefidx] = KernelBy[coefidx];
        MatrixDiffVector2[local + 2 * OwnedDiffVector.size()][coefidx] =
            KernelBz[coefidx];

        Matrix[LocalDiffVectorBase + local][coefidx] =
            MatrixDiffVector1[local][coefidx] - MatrixDiffVector2[local][coefidx];
        Matrix[LocalDiffVectorBase + local + OwnedDiffVector.size()][coefidx] =
            MatrixDiffVector1[local + OwnedDiffVector.size()][coefidx] -
            MatrixDiffVector2[local + OwnedDiffVector.size()][coefidx];
        Matrix[LocalDiffVectorBase + local + 2 * OwnedDiffVector.size()][coefidx] =
            MatrixDiffVector1[local + 2 * OwnedDiffVector.size()][coefidx] -
            MatrixDiffVector2[local + 2 * OwnedDiffVector.size()][coefidx];
      }
    }
    MPIParallel::phase_complete("differenced vector kernels", kernel_started);
    RSCHAConsole::done("Kernel block complete");
  }

  RSCHAConsole::row("Root-owned scalar rows ready", NumLocal);
  std::fflush(stderr);
  std::vector<double> d(NumData), Weight_Data(NumData), Sigma(NumData);

  Count = 0;
  for (int i = 0; i < NumTMI; ++i) {
    d[Count] = TMI[i];
    Weight_Data[Count] = DataWeight_TMI[i];
    Sigma[Count] = Sigma_TMI[i];
    ++Count;
  }
  for (int i = 0; i < NumVector; ++i) {
    d[Count] = Bx_SCCS[i];
    Weight_Data[Count] = DataWeight_Bx[i];
    Sigma[Count] = Sigma_Bx[i];
    ++Count;
  }
  for (int i = 0; i < NumVector; ++i) {
    d[Count] = By_SCCS[i];
    Weight_Data[Count] = DataWeight_By[i];
    Sigma[Count] = Sigma_By[i];
    ++Count;
  }
  for (int i = 0; i < NumVector; ++i) {
    d[Count] = Bz_SCCS[i];
    Weight_Data[Count] = DataWeight_Bz[i];
    Sigma[Count] = Sigma_Bz[i];
    ++Count;
  }
  for (int i = 0; i < NumDiffTMI; ++i) {
    d[Count] = DiffTMI1[i] - DiffTMI2[i];
    Weight_Data[Count] = DataWeight_DiffTMI[i];
    Sigma[Count] = Sigma_DiffTMI[i];
    ++Count;
  }
  for (int i = 0; i < NumDiffVector; ++i) {
    d[Count] = DiffBx1_SCCS[i] - DiffBx2_SCCS[i];
    Weight_Data[Count] = DataWeight_DiffBx[i];
    Sigma[Count] = Sigma_DiffBx[i];
    ++Count;
  }
  for (int i = 0; i < NumDiffVector; ++i) {
    d[Count] = DiffBy1_SCCS[i] - DiffBy2_SCCS[i];
    Weight_Data[Count] = DataWeight_DiffBy[i];
    Sigma[Count] = Sigma_DiffBy[i];
    ++Count;
  }
  for (int i = 0; i < NumDiffVector; ++i) {
    d[Count] = DiffBz1_SCCS[i] - DiffBz2_SCCS[i];
    Weight_Data[Count] = DataWeight_DiffBz[i];
    Sigma[Count] = Sigma_DiffBz[i];
    ++Count;
  }

  std::vector<double> Coefficients(Num_Coef);
  std::vector<double> Misfit(Max_Iteration);
  std::vector<double> d_Predict(NumData);
  std::vector<double> LocalData(NumLocal), LocalWeight(NumLocal),
      LocalSigma(NumLocal);
  for (int local = 0; local < NumLocal; ++local) {
    const int global = GlobalRows[local];
    LocalData[local] = d[global];
    LocalWeight[local] = Weight_Data[global];
    LocalSigma[local] = Sigma[global];
  }

  RSCHAConsole::section(5, "Inversion iterations");
  PreCond_BICGSTAB(NumData, Num_Coef, Matrix, GlobalRows, m0, LocalData,
                   epsilon, damping, Max_Iteration, LocalWeight, Weight_Model,
                   Prior_Model, LocalSigma, RobustSwitch, RobustConstant,
                   Coefficients, Misfit, d_Predict);

  std::vector<double> ReconsDiffTMI1All(NumDiffTMI, 0.0),
      ReconsDiffTMI2All(NumDiffTMI, 0.0);
  std::vector<double> ReconsDiffVector1All(3 * NumDiffVector, 0.0),
      ReconsDiffVector2All(3 * NumDiffVector, 0.0);
  for (int local = 0; local < static_cast<int>(OwnedDiffTMI.size()); ++local) {
    const int global = OwnedDiffTMI[local];
    for (int j = 0; j < Num_Coef; ++j) {
      ReconsDiffTMI1All[global] += MatrixDiffTMI1[local][j] * Coefficients[j];
      ReconsDiffTMI2All[global] += MatrixDiffTMI2[local][j] * Coefficients[j];
    }
  }
  for (int local = 0; local < static_cast<int>(OwnedDiffVector.size()); ++local) {
    const int global = OwnedDiffVector[local];
    for (int component = 0; component < 3; ++component) {
      const int local_row = component * OwnedDiffVector.size() + local;
      const int global_row = component * NumDiffVector + global;
      for (int j = 0; j < Num_Coef; ++j) {
        ReconsDiffVector1All[global_row] +=
            MatrixDiffVector1[local_row][j] * Coefficients[j];
        ReconsDiffVector2All[global_row] +=
            MatrixDiffVector2[local_row][j] * Coefficients[j];
      }
    }
  }
  MPIParallel::allreduce_sum(ReconsDiffTMI1All);
  MPIParallel::allreduce_sum(ReconsDiffTMI2All);
  MPIParallel::allreduce_sum(ReconsDiffVector1All);
  MPIParallel::allreduce_sum(ReconsDiffVector2All);

  std::vector<double> n, Spectrum;
  RSCHA_Spectrum(Coefficients, K_Max_Internal, K_Max_External, P_Max, M_Max, nk,
                 Theta0, Min_Radius, Min_Radius, n, Spectrum);

  if (!mpi.is_root())
    return 0;

  RSCHAConsole::section(6, "Reconstruction and spectrum");
  std::vector<double> ReconsTMI(NumTMI, 0.0);
  std::vector<double> ReconsBx(NumVector, 0.0), ReconsBy(NumVector, 0.0),
      ReconsBz(NumVector, 0.0);
  std::vector<double> ReconsBx_SCCS(NumVector, 0.0),
      ReconsBy_SCCS(NumVector, 0.0), ReconsBz_SCCS(NumVector, 0.0);
  std::vector<double> ReconsDiffTMI1(NumDiffTMI, 0.0),
      ReconsDiffTMI2(NumDiffTMI, 0.0);
  std::vector<double> ReconsDiffBx1(NumDiffVector, 0.0),
      ReconsDiffBy1(NumDiffVector, 0.0), ReconsDiffBz1(NumDiffVector, 0.0);
  std::vector<double> ReconsDiffBx2(NumDiffVector, 0.0),
      ReconsDiffBy2(NumDiffVector, 0.0), ReconsDiffBz2(NumDiffVector, 0.0);
  std::vector<double> ReconsDiffBx1_SCCS(NumDiffVector, 0.0),
      ReconsDiffBy1_SCCS(NumDiffVector, 0.0),
      ReconsDiffBz1_SCCS(NumDiffVector, 0.0);
  std::vector<double> ReconsDiffBx2_SCCS(NumDiffVector, 0.0),
      ReconsDiffBy2_SCCS(NumDiffVector, 0.0),
      ReconsDiffBz2_SCCS(NumDiffVector, 0.0);

  Count = 0;
  for (int i = 0; i < NumTMI; ++i) {
    ReconsTMI[i] = d_Predict[Count];
    ++Count;
  }
  for (int i = 0; i < NumVector; ++i) {
    ReconsBx_SCCS[i] = d_Predict[Count];
    ++Count;
  }
  for (int i = 0; i < NumVector; ++i) {
    ReconsBy_SCCS[i] = d_Predict[Count];
    ++Count;
  }
  for (int i = 0; i < NumVector; ++i) {
    ReconsBz_SCCS[i] = d_Predict[Count];
    ++Count;
  }

  ReconsDiffTMI1 = ReconsDiffTMI1All;
  ReconsDiffTMI2 = ReconsDiffTMI2All;

  for (int i = 0; i < NumDiffVector; ++i) {
    ReconsDiffBx1_SCCS[i] = ReconsDiffVector1All[i];
    ReconsDiffBy1_SCCS[i] = ReconsDiffVector1All[i + NumDiffVector];
    ReconsDiffBz1_SCCS[i] = ReconsDiffVector1All[i + 2 * NumDiffVector];
    ReconsDiffBx2_SCCS[i] = ReconsDiffVector2All[i];
    ReconsDiffBy2_SCCS[i] = ReconsDiffVector2All[i + NumDiffVector];
    ReconsDiffBz2_SCCS[i] = ReconsDiffVector2All[i + 2 * NumDiffVector];
  }

  SCHA_3CompsTrans(LongitudeVector, LatitudeVector, Lon0, Lat0, ReconsBx_SCCS,
                   ReconsBy_SCCS, ReconsBz_SCCS, 2, ReconsBx, ReconsBy,
                   ReconsBz);
  SCHA_3CompsTrans(LongitudeDiffVector1, LatitudeDiffVector1, Lon0, Lat0,
                   ReconsDiffBx1_SCCS, ReconsDiffBy1_SCCS, ReconsDiffBz1_SCCS,
                   2, ReconsDiffBx1, ReconsDiffBy1, ReconsDiffBz1);
  SCHA_3CompsTrans(LongitudeDiffVector2, LatitudeDiffVector2, Lon0, Lat0,
                   ReconsDiffBx2_SCCS, ReconsDiffBy2_SCCS, ReconsDiffBz2_SCCS,
                   2, ReconsDiffBx2, ReconsDiffBy2, ReconsDiffBz2);

  RSCHAConsole::section(7, "Save model and reconstructed data");
  if (!SaveModellingResults(

      ReconsTMIFile, ReconsVectorFile, ReconsDiffTMIFile, ReconsDiffVectorFile,

      ModellingResultsFile, ModelSpectrumFile,

      IncludeTMI, IncludeVector, IncludeDiffTMI, IncludeDiffVector,

      Lon0, Lat0, Theta0, K_Max_Internal, K_Max_External, P_Max, M_Max,
      Max_Radius, Min_Radius,

      LongitudeTMI, LatitudeTMI, RadiusTMI, TMI, ReconsTMI,

      LongitudeVector, LatitudeVector, RadiusVector, Bx, By, Bz, ReconsBx,
      ReconsBy, ReconsBz,

      LongitudeDiffTMI1, LatitudeDiffTMI1, RadiusDiffTMI1, DiffTMI1,
      ReconsDiffTMI1, LongitudeDiffTMI2, LatitudeDiffTMI2, RadiusDiffTMI2,
      DiffTMI2, ReconsDiffTMI2,

      LongitudeDiffVector1, LatitudeDiffVector1, RadiusDiffVector1, DiffBx1,
      DiffBy1, DiffBz1, ReconsDiffBx1, ReconsDiffBy1, ReconsDiffBz1,
      LongitudeDiffVector2, LatitudeDiffVector2, RadiusDiffVector2, DiffBx2,
      DiffBy2, DiffBz2, ReconsDiffBx2, ReconsDiffBy2, ReconsDiffBz2,
      Coefficients, n, Spectrum))
    return 1;

  RSCHAConsole::complete("R-SCHA inversion complete");

  return 0;
}
catch (const std::exception &error) {
  std::cerr << "Calculation failed on rank " << mpi.rank() << ": " << error.what() << std::endl;
  MPI_Abort(MPI_COMM_WORLD, 1);
  return 1;
}
}
