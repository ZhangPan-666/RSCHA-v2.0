#include "Basis Function.h"
#include "RSCHAConsole.h"
#include "CoordinateTransformation.h"
#include "MPIParallel.h"
#include "R_SCHA.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <vector>

using namespace CoordinateTransformation;
using namespace R_SCHA;

int main(int argc, char **argv) {

  MPIParallel::Environment mpi(argc, argv);
  try {
  if (!mpi.is_root())
    std::cout.setstate(std::ios_base::failbit);

  RSCHAConsole::banner("R-SCHA V/B/T forward calculation", "MPI");
  RSCHAConsole::section(1, "Model and coordinate input");

  const char ModellingParametersFile[] = "Results_RSCHAModel.dat";
  const char CoordinatesFile[] = "Parameters_Coordinates.dat";
  const char ForwardResultsFile[] = "Results_ForwardResults.dat";

  vector<vector<double>> MatrixTemp;

  if (readmatrix(ModellingParametersFile, MatrixTemp) != 1) {
    cout << "An error occurred while reading the modelling parameters file."
         << endl;
    throw std::runtime_error("Input validation failed.");
  }

  if (MatrixTemp.size() <= 9 || MatrixTemp.front().size() != 1) {
    cout << "An error occurred while reading the modelling parameters file."
         << endl;
    cout << "Please check the matrix dimensions in the input text file. "
            "(Expected: 9 parameters followed by coefficients, one column)"
         << endl;
    cout << ModellingParametersFile << endl;
    throw std::runtime_error("Input validation failed.");
  }

  const double Lon0 = MatrixTemp[0][0];
  const double Lat0 = MatrixTemp[1][0];
  const double Theta0 = MatrixTemp[2][0];
  const int K_Max_Internal = checked_index(MatrixTemp[3][0]);
  const int K_Max_External = checked_index(MatrixTemp[4][0]);
  const int P_Max = checked_index(MatrixTemp[5][0]);
  const int M_Max = checked_index(MatrixTemp[6][0]);
  const double Max_Radius = MatrixTemp[7][0];
  const double Min_Radius = MatrixTemp[8][0];
  check_geometry(Lon0, Lat0, Theta0, Min_Radius, Max_Radius);
  const int Num_Coef = coefficient_count(K_Max_Internal, K_Max_External, P_Max, M_Max);

  const int Num_CoefIn = (K_Max_Internal + 1) * (K_Max_Internal + 2) / 2;
  const int Num_CoefEx = (K_Max_External + 1) * (K_Max_External + 2) / 2;
  const int Num_CoefMe = P_Max * (P_Max + 1);
  const int Num_CoefCo = M_Max + 1;

  if (MatrixTemp.size() != static_cast<std::size_t>(Num_Coef) + 9) {
    cout << "An error occurred while reading the modelling parameters file."
         << endl;
    cout << "The modelling parameters file has been damaged" << endl;
    cout << ModellingParametersFile << endl;
    throw std::runtime_error("Input validation failed.");
  }

  vector<double> Coefficients(Num_Coef);
  for (int i = 0; i < Num_Coef; ++i) {
    Coefficients[i] = MatrixTemp[i + 9][0];
  }

  MatrixTemp.clear();
  MatrixTemp.shrink_to_fit();

  cout << "Modelling parameters loaded from:" << "\n     "
       << ModellingParametersFile << endl
       << endl;

  if (readmatrix(CoordinatesFile, MatrixTemp) != 1) {
    cout << "An error occurred while reading the forward coordinate file."
         << endl;
    throw std::runtime_error("Input validation failed.");
  }

  if (MatrixTemp.front().size() != 3) {
    cout << "An error occurred while reading the forward coordinate file."
         << endl;
    cout << "Please check the matrix dimensions in the input text file. "
            "(Expected: 3 columns)"
         << endl;
    cout << CoordinatesFile << endl;
    throw std::runtime_error("Input validation failed.");
  }

  int Num_Coordinate = static_cast<int>(MatrixTemp.size());
  vector<double> Longitude(Num_Coordinate), Latitude(Num_Coordinate),
      Radius(Num_Coordinate);

  for (int i = 0; i < Num_Coordinate; ++i) {
    Longitude[i] = MatrixTemp[i][0];
    Latitude[i] = MatrixTemp[i][1];
    Radius[i] = MatrixTemp[i][2];
    if (Longitude[i] < -180.0 || Longitude[i] > 180.0 ||
        Latitude[i] < -90.0 || Latitude[i] > 90.0 || Radius[i] <= 0.0)
      throw std::invalid_argument("Invalid forward coordinates.");
  }

  MatrixTemp.clear();
  MatrixTemp.shrink_to_fit();

  cout << "Forward coordinates loaded from:" << "\n     " << CoordinatesFile
       << endl;

  RSCHAConsole::section(2, "Model and coordinate settings");
  RSCHAConsole::model(Lon0, Lat0, Theta0, K_Max_Internal, K_Max_External,
                      P_Max, M_Max, Min_Radius, Max_Radius);
  RSCHAConsole::row("Input coordinates", Num_Coordinate);
  RSCHAConsole::row("Model coefficients", Num_Coef);
  RSCHAConsole::row("OpenMP threads per process", omp_get_max_threads());
  RSCHAConsole::section(3, "Regional basis and coordinates");

  const int K_Max = std::max(K_Max_Internal, K_Max_External);
  std::vector<double> m_Legendre(std::max(Num_CoefIn, Num_CoefEx));
  std::vector<double> nk(m_Legendre.size());
  int Count = 0;
  for (int i = 0; i <= K_Max; ++i)
    for (int j = 0; j <= i; ++j)
      m_Legendre[Count++] = j;
  R_SCHA_nk(K_Max, Theta0, nk);

  constexpr int bulk_preparation_threshold = 32;
  const int local_coordinate_count =
      Num_Coordinate <= mpi.rank()
          ? 0
          : 1 + (Num_Coordinate - 1 - mpi.rank()) / mpi.size();
  if (local_coordinate_count >= bulk_preparation_threshold)
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

  vector<double> Lamda(Num_Coordinate), Theta(Num_Coordinate);
  SCHA_CT(Longitude, Latitude, Lon0, Lat0, 1, Lamda, Theta);

  int retained = 0;
  for (int i = 0; i < Num_Coordinate; ++i) {
    if (Theta[i] >= Theta0 || Theta[i] < 1.0e-3 ||
        Radius[i] < Min_Radius || Radius[i] > Max_Radius)
      continue;
    if (retained != i) {
      Longitude[retained] = Longitude[i];
      Latitude[retained] = Latitude[i];
      Lamda[retained] = Lamda[i];
      Theta[retained] = Theta[i];
      Radius[retained] = Radius[i];
    }
    ++retained;
  }
  const int EraseCount = Num_Coordinate - retained;
  Num_Coordinate = retained;
  Longitude.resize(retained);
  Latitude.resize(retained);
  Lamda.resize(retained);
  Theta.resize(retained);
  Radius.resize(retained);

  if (EraseCount > 0) {
    cout << endl;
    cout << EraseCount
         << " Forward data points were outside the modelling area or " << endl;
    cout << "    too close to the centre of the modeling area and have been "
            "deleted."
         << endl;
  }

  if (Num_Coordinate == 0)
    throw std::invalid_argument("No forward coordinates remain in the modelling region.");

  vector<double> KernelV(Num_Coef);
  vector<double> KernelBx(Num_Coef), KernelBy(Num_Coef), KernelBz(Num_Coef);
  vector<double> KernelTxx(Num_Coef), KernelTxy(Num_Coef), KernelTyy(Num_Coef);
  vector<double> KernelTzx(Num_Coef), KernelTzy(Num_Coef), KernelTzz(Num_Coef);

  vector<double> V_RSCHA(Num_Coordinate);
  vector<double> Bx_RSCHA(Num_Coordinate), By_RSCHA(Num_Coordinate),
      Bz_RSCHA(Num_Coordinate);
  vector<double> Txx_RSCHA(Num_Coordinate), Txy_RSCHA(Num_Coordinate),
      Tyy_RSCHA(Num_Coordinate);
  vector<double> Tzx_RSCHA(Num_Coordinate), Tzy_RSCHA(Num_Coordinate),
      Tzz_RSCHA(Num_Coordinate);

  RSCHAConsole::section(4, "V/B/T forward calculation");
  RSCHAConsole::row("Evaluation coordinates", Num_Coordinate);
  RSCHAConsole::row("MPI processes", mpi.size());

  const std::vector<int> OwnedCoordinates = MPIParallel::owned_indices(Num_Coordinate, mpi.rank(), mpi.size());
  MPIParallel::phase_begin("forward local VBT calculation");
  double forward_started = MPIParallel::time();

  // Every rank enters the same progress collectives, including ranks with no rows.
  const int max_local_coordinate_count =
      Num_Coordinate / mpi.size() + (Num_Coordinate % mpi.size() != 0);
  const int progress_interval = (max_local_coordinate_count - 1) / 25 + 1;
  if (mpi.is_root()) {
    std::printf("[Forward progress] 0/%d (0.0%%); starting\n", Num_Coordinate);
    std::fflush(stdout);
  }

  for (int local_position = 0; local_position < max_local_coordinate_count; ++local_position) {
    if (static_cast<std::size_t>(local_position) < OwnedCoordinates.size()) {
      const int dataidx = OwnedCoordinates[local_position];

      RSCHA_TensorKernel(Radius[dataidx], Theta[dataidx], Theta0, Lamda[dataidx],
                         Min_Radius, Max_Radius, Num_CoefIn, Num_CoefEx, Num_CoefMe,
                         Num_CoefCo, nk, m_Legendre, p, m_Mehler, m_Constant,
                         KernelV, KernelBx, KernelBy, KernelBz, KernelTxx,
                         KernelTxy, KernelTyy, KernelTzx, KernelTzy, KernelTzz);

      double v = 0.0, bx = 0.0, by = 0.0, bz = 0.0;
      double txx = 0.0, txy = 0.0, tyy = 0.0;
      double tzx = 0.0, tzy = 0.0, tzz = 0.0;
      for (int coefidx = 0; coefidx < Num_Coef; ++coefidx) {
        const double coefficient = Coefficients[coefidx];
        v += KernelV[coefidx] * coefficient;
        bx += KernelBx[coefidx] * coefficient;
        by += KernelBy[coefidx] * coefficient;
        bz += KernelBz[coefidx] * coefficient;
        txx += KernelTxx[coefidx] * coefficient;
        txy += KernelTxy[coefidx] * coefficient;
        tyy += KernelTyy[coefidx] * coefficient;
        tzx += KernelTzx[coefidx] * coefficient;
        tzy += KernelTzy[coefidx] * coefficient;
        tzz += KernelTzz[coefidx] * coefficient;
      }
      V_RSCHA[dataidx] = v;
      Bx_RSCHA[dataidx] = bx;
      By_RSCHA[dataidx] = by;
      Bz_RSCHA[dataidx] = bz;
      Txx_RSCHA[dataidx] = txx;
      Txy_RSCHA[dataidx] = txy;
      Tyy_RSCHA[dataidx] = tyy;
      Tzx_RSCHA[dataidx] = tzx;
      Tzy_RSCHA[dataidx] = tzy;
      Tzz_RSCHA[dataidx] = tzz;
    }

    if ((local_position + 1) % progress_interval == 0 ||
        local_position + 1 == max_local_coordinate_count) {
      const int local_completed = std::min(
          local_position + 1, static_cast<int>(OwnedCoordinates.size()));
      int global_completed = 0;
      MPIParallel::check(MPI_Allreduce(&local_completed, &global_completed, 1,
          MPI_INT, MPI_SUM, MPI_COMM_WORLD), "MPI_Allreduce forward progress");
      if (mpi.is_root()) {
        const double elapsed = MPIParallel::time() - forward_started;
        const double fraction = static_cast<double>(global_completed) /
            static_cast<double>(Num_Coordinate);
        const int bar_width = 40;
        const int filled = std::min(bar_width,
            static_cast<int>(fraction * static_cast<double>(bar_width) + 0.5));
        const double rate = elapsed > 0.0 ? global_completed / elapsed : 0.0;
        const double eta = rate > 0.0 ? (Num_Coordinate - global_completed) / rate : 0.0;
        std::printf("[Forward progress] [");
        for (int i = 0; i < bar_width; ++i) std::printf(i < filled ? "#" : "-");
        std::printf("] %d/%d (%6.2f%%); %.1f points/s; ETA %.1f min\n",
            global_completed, Num_Coordinate, 100.0 * fraction, rate, eta / 60.0);
        std::fflush(stdout);
      }
    }
  }
  MPIParallel::phase_complete("forward local VBT calculation", forward_started);
  MPIParallel::phase_begin("forward owned-value gather and broadcast");
  forward_started = MPIParallel::time();
  MPIParallel::IndexedGather gather(Num_Coordinate, OwnedCoordinates);
  std::vector<double> local_values(OwnedCoordinates.size());
  std::vector<double> *fields[] = {&V_RSCHA, &Bx_RSCHA, &By_RSCHA, &Bz_RSCHA,
      &Txx_RSCHA, &Txy_RSCHA, &Tyy_RSCHA, &Tzx_RSCHA, &Tzy_RSCHA, &Tzz_RSCHA};
  for (std::vector<double> *field : fields) {
    for (std::size_t i = 0; i < OwnedCoordinates.size(); ++i)
      local_values[i] = (*field)[OwnedCoordinates[i]];
    gather.collect(local_values, *field);
  }
  MPIParallel::phase_complete("forward owned-value gather and broadcast", forward_started);

  if (!mpi.is_root())
    return 0;
  RSCHAConsole::done("Forward kernels complete");
  RSCHAConsole::section(5, "Rotate components and save results");
  vector<double> V_LocalCart(Num_Coordinate);
  vector<double> Bx_LocalCart(Num_Coordinate), By_LocalCart(Num_Coordinate),
      Bz_LocalCart(Num_Coordinate);
  vector<double> Txx_LocalCart(Num_Coordinate), Txy_LocalCart(Num_Coordinate),
      Tyy_LocalCart(Num_Coordinate);
  vector<double> Tzx_LocalCart(Num_Coordinate), Tzy_LocalCart(Num_Coordinate),
      Tzz_LocalCart(Num_Coordinate);

  Rotate2RankTensor(Lon0, Lat0, Longitude, Latitude, Lamda, Theta, V_RSCHA,
                    Bx_RSCHA, By_RSCHA, Bz_RSCHA, Txx_RSCHA, Txy_RSCHA,
                    Tyy_RSCHA, Tzx_RSCHA, Tzy_RSCHA, Tzz_RSCHA, V_LocalCart,
                    Bx_LocalCart, By_LocalCart, Bz_LocalCart, Txx_LocalCart,
                    Txy_LocalCart, Tyy_LocalCart, Tzx_LocalCart, Tzy_LocalCart,
                    Tzz_LocalCart);

  ofstream fp_ForwardResults(ForwardResultsFile, ios::out);
  if (!fp_ForwardResults) {
    cout << "Fail to open file: " << ForwardResultsFile << "\n";
    cout << "An error occurred while saving the file." << endl;
    return 1;
  } else {
    fp_ForwardResults.precision(17);
    for (int dataidx = 0; dataidx < Num_Coordinate; ++dataidx) {
      fp_ForwardResults << Longitude[dataidx] << ",";
      fp_ForwardResults << Latitude[dataidx] << ",";
      fp_ForwardResults << Radius[dataidx] << ",";
      fp_ForwardResults << V_LocalCart[dataidx] * 1.0e-6 << ",";
      fp_ForwardResults << Bx_LocalCart[dataidx] << ",";
      fp_ForwardResults << By_LocalCart[dataidx] << ",";
      fp_ForwardResults << Bz_LocalCart[dataidx] << ",";
      fp_ForwardResults << Txx_LocalCart[dataidx] * 1.0e3 << ",";
      fp_ForwardResults << Txy_LocalCart[dataidx] * 1.0e3 << ",";
      fp_ForwardResults << Tyy_LocalCart[dataidx] * 1.0e3 << ",";
      fp_ForwardResults << Tzx_LocalCart[dataidx] * 1.0e3 << ",";
      fp_ForwardResults << Tzy_LocalCart[dataidx] * 1.0e3 << ",";
      fp_ForwardResults << Tzz_LocalCart[dataidx] * 1.0e3 << '\n';
    }
    fp_ForwardResults.close();
    if (!fp_ForwardResults)
      return 1;

    cout << "\nForwarded 0-, 1- and 2-rank gradient tensor data stored in: "
         << ForwardResultsFile << endl;
  }
  RSCHAConsole::complete("R-SCHA forward complete");

  return 0;
}
catch (const std::exception &error) {
  std::cerr << "Calculation failed: " << error.what() << std::endl;
  MPI_Abort(MPI_COMM_WORLD, 1);
  return 1;
}
}
