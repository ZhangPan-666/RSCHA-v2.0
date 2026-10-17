#include "Basis Function.h"
#include "RSCHAConsole.h"
#include "CoordinateTransformation.h"
#include "R_SCHA.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <vector>

using namespace CoordinateTransformation;
using namespace R_SCHA;

int main() try {

  RSCHAConsole::banner("R-SCHA V/B/T forward calculation", "OpenMP");
  RSCHAConsole::section(1, "Model and coordinate input");

  const char ModellingParametersFile[] = "Results_RSCHAModel.dat";
  const char CoordinatesFile[] = "Parameters_Coordinates.dat";
  const char ForwardResultsFile[] = "Results_ForwardResults.dat";

  vector<vector<double>> MatrixTemp;

  if (readmatrix(ModellingParametersFile, MatrixTemp) != 1) {
    cout << "An error occurred while reading the modelling parameters file."
         << endl;
    return -13;
  }

  if (MatrixTemp.size() <= 9 || MatrixTemp.front().size() != 1) {
    cout << "An error occurred while reading the modelling parameters file."
         << endl;
    cout << "Please check the matrix dimensions in the input text file. "
            "(Expected: 9 parameters followed by coefficients, one column)"
         << endl;
    cout << ModellingParametersFile << endl;
    return -13;
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
    return -13;
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
    return -14;
  }

  if (MatrixTemp.front().size() != 3) {
    cout << "An error occurred while reading the forward coordinate file."
         << endl;
    cout << "Please check the matrix dimensions in the input text file. "
            "(Expected: 3 columns)"
         << endl;
    cout << CoordinatesFile << endl;
    return -14;
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
  if (Num_Coordinate >= bulk_preparation_threshold) {
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < static_cast<int>(nk.size()); ++i)
      Legendre::PrepareLegacySchmidtMode(m_Legendre[i], nk[i]);
  }

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

  vector<double> V_RSCHA(Num_Coordinate);
  vector<double> Bx_RSCHA(Num_Coordinate), By_RSCHA(Num_Coordinate),
      Bz_RSCHA(Num_Coordinate);
  vector<double> Txx_RSCHA(Num_Coordinate), Txy_RSCHA(Num_Coordinate),
      Tyy_RSCHA(Num_Coordinate);
  vector<double> Tzx_RSCHA(Num_Coordinate), Tzy_RSCHA(Num_Coordinate),
      Tzz_RSCHA(Num_Coordinate);

  vector<double> V_LocalCart(Num_Coordinate);
  vector<double> Bx_LocalCart(Num_Coordinate), By_LocalCart(Num_Coordinate),
      Bz_LocalCart(Num_Coordinate);
  vector<double> Txx_LocalCart(Num_Coordinate), Txy_LocalCart(Num_Coordinate),
      Tyy_LocalCart(Num_Coordinate);
  vector<double> Tzx_LocalCart(Num_Coordinate), Tzy_LocalCart(Num_Coordinate),
      Tzz_LocalCart(Num_Coordinate);

  RSCHAConsole::section(4, "V/B/T forward calculation");
  RSCHAConsole::row("Evaluation coordinates", Num_Coordinate);

  const int progress_interval = (Num_Coordinate - 1) / 25 + 1;
  const auto forward_started = std::chrono::steady_clock::now();
  std::atomic<int> progress_completed(0);
  int progress_last_reported = 0;
  std::printf("[Forward progress] 0/%d (0.0%%); starting\n", Num_Coordinate);
  std::fflush(stdout);
  auto mark_progress = [&]() {
    const int completed = progress_completed.fetch_add(1, std::memory_order_relaxed) + 1;
    if (completed % progress_interval != 0 && completed != Num_Coordinate)
      return;
#pragma omp critical(rscha_forward_progress)
    {
      if (completed > progress_last_reported) {
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - forward_started).count();
        const double fraction = static_cast<double>(completed) / Num_Coordinate;
        const int filled = std::min(40, static_cast<int>(40 * fraction + 0.5));
        const double rate = elapsed > 0.0 ? completed / elapsed : 0.0;
        const double eta = rate > 0.0 ? (Num_Coordinate - completed) / rate : 0.0;
        std::printf("[Forward progress] [");
        for (int i = 0; i < 40; ++i)
          std::printf(i < filled ? "#" : "-");
        std::printf("] %d/%d (%6.2f%%); %.1f points/s; ETA %.1f min\n",
            completed, Num_Coordinate, 100 * fraction, rate, eta / 60);
        std::fflush(stdout);
        progress_last_reported = completed;
      }
    }
  };

  auto forward_coordinate =
      [&](int dataidx, vector<double> &KernelV, vector<double> &KernelBx,
          vector<double> &KernelBy, vector<double> &KernelBz,
          vector<double> &KernelTxx, vector<double> &KernelTxy,
          vector<double> &KernelTyy, vector<double> &KernelTzx,
          vector<double> &KernelTzy, vector<double> &KernelTzz) {
        RSCHA_TensorKernel(
            Radius[dataidx], Theta[dataidx], Theta0, Lamda[dataidx],
            Min_Radius, Max_Radius, Num_CoefIn, Num_CoefEx, Num_CoefMe,
            Num_CoefCo, nk, m_Legendre, p, m_Mehler, m_Constant, KernelV,
            KernelBx, KernelBy, KernelBz, KernelTxx, KernelTxy, KernelTyy,
            KernelTzx, KernelTzy, KernelTzz);

        double sumV = 0.0, sumBx = 0.0, sumBy = 0.0, sumBz = 0.0;
        double sumTxx = 0.0, sumTxy = 0.0, sumTyy = 0.0;
        double sumTzx = 0.0, sumTzy = 0.0, sumTzz = 0.0;
#pragma omp simd reduction(+ : sumV, sumBx, sumBy, sumBz, sumTxx, sumTxy,     \
                               sumTyy, sumTzx, sumTzy, sumTzz)
        for (int coefidx = 0; coefidx < Num_Coef; ++coefidx) {
          const double coefficient = Coefficients[coefidx];
          sumV += KernelV[coefidx] * coefficient;
          sumBx += KernelBx[coefidx] * coefficient;
          sumBy += KernelBy[coefidx] * coefficient;
          sumBz += KernelBz[coefidx] * coefficient;
          sumTxx += KernelTxx[coefidx] * coefficient;
          sumTxy += KernelTxy[coefidx] * coefficient;
          sumTyy += KernelTyy[coefidx] * coefficient;
          sumTzx += KernelTzx[coefidx] * coefficient;
          sumTzy += KernelTzy[coefidx] * coefficient;
          sumTzz += KernelTzz[coefidx] * coefficient;
        }
        V_RSCHA[dataidx] = sumV;
        Bx_RSCHA[dataidx] = sumBx;
        By_RSCHA[dataidx] = sumBy;
        Bz_RSCHA[dataidx] = sumBz;
        Txx_RSCHA[dataidx] = sumTxx;
        Txy_RSCHA[dataidx] = sumTxy;
        Tyy_RSCHA[dataidx] = sumTyy;
        Tzx_RSCHA[dataidx] = sumTzx;
        Tzy_RSCHA[dataidx] = sumTzy;
        Tzz_RSCHA[dataidx] = sumTzz;
      };

  omp_set_max_active_levels(1);
  if (Num_Coordinate < omp_get_max_threads()) {
    vector<double> KernelV(Num_Coef), KernelBx(Num_Coef),
        KernelBy(Num_Coef), KernelBz(Num_Coef), KernelTxx(Num_Coef),
        KernelTxy(Num_Coef), KernelTyy(Num_Coef), KernelTzx(Num_Coef),
        KernelTzy(Num_Coef), KernelTzz(Num_Coef);
    for (int dataidx = 0; dataidx < Num_Coordinate; ++dataidx) {
      forward_coordinate(dataidx, KernelV, KernelBx, KernelBy, KernelBz,
                         KernelTxx, KernelTxy, KernelTyy, KernelTzx,
                         KernelTzy, KernelTzz);
      mark_progress();
    }
  } else {
#pragma omp parallel
    {
      vector<double> KernelV(Num_Coef), KernelBx(Num_Coef),
          KernelBy(Num_Coef), KernelBz(Num_Coef), KernelTxx(Num_Coef),
          KernelTxy(Num_Coef), KernelTyy(Num_Coef), KernelTzx(Num_Coef),
          KernelTzy(Num_Coef), KernelTzz(Num_Coef);

#pragma omp for schedule(dynamic, 1)
      for (int dataidx = 0; dataidx < Num_Coordinate; ++dataidx) {
        forward_coordinate(dataidx, KernelV, KernelBx, KernelBy, KernelBz,
                           KernelTxx, KernelTxy, KernelTyy, KernelTzx,
                           KernelTzy, KernelTzz);
        mark_progress();
      }
    }
  }
  RSCHAConsole::done("Forward kernels complete");
  RSCHAConsole::section(5, "Rotate components and save results");
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
  return 1;
}
