#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>

namespace CoordinateTransformation {

using namespace DegreeTrig;

const double rad = 4.0 * atan(1.0) / 180.0;

inline void SCHA_CT(const std::vector<double> &Lon_Original,
                    const std::vector<double> &Lat_Original, double Lon0,
                    double Lat0, int CTType, std::vector<double> &Lon_Converted,
                    std::vector<double> &Lat_Converted) {
  const std::size_t Num_Coordinate = Lon_Original.size();

  if (CTType == 1) {

    for (std::size_t i = 0; i < Num_Coordinate; ++i) {
      const double x = cosd(Lat_Original[i]) * cosd(Lon_Original[i]);
      const double y = cosd(Lat_Original[i]) * sind(Lon_Original[i]);
      const double z = sind(Lat_Original[i]);

      const double x1 = x * cosd(Lon0) + y * sind(Lon0);
      const double y1 = -x * sind(Lon0) + y * cosd(Lon0);
      const double z1 = z;

      const double x2 = x1 * cosd(90.0 - Lat0) - z1 * sind(90.0 - Lat0);
      const double y2 = y1;
      const double z2 = x1 * sind(90.0 - Lat0) + z1 * cosd(90.0 - Lat0);

      Lon_Converted[i] = std::atan2(y2, x2) / rad;
      Lat_Converted[i] = (2.0 * std::atan(1.0) -
                          std::atan2(z2, std::sqrt(x2 * x2 + y2 * y2))) /
                         rad;
    }
  } else if (CTType == 2) {

    for (std::size_t i = 0; i < Num_Coordinate; ++i) {
      const double x = cosd(90.0 - Lat_Original[i]) * cosd(Lon_Original[i]);
      const double y = cosd(90.0 - Lat_Original[i]) * sind(Lon_Original[i]);
      const double z = sind(90.0 - Lat_Original[i]);

      const double x1 = x * cosd(90.0 - Lat0) + z * sind(90.0 - Lat0);
      const double y1 = y;
      const double z1 = -x * sind(90.0 - Lat0) + z * cosd(90.0 - Lat0);

      const double x2 = x1 * cosd(Lon0) - y1 * sind(Lon0);
      const double y2 = x1 * sind(Lon0) + y1 * cosd(Lon0);
      const double z2 = z1;

      Lon_Converted[i] = std::atan2(y2, x2) / rad;
      Lat_Converted[i] = std::atan2(z2, std::sqrt(x2 * x2 + y2 * y2)) / rad;
    }
  }
}

inline void SCHA_3CompsTrans(const std::vector<double> &Longitude,
                             const std::vector<double> &Latitude, double Lon0,
                             double Lat0,
                             const std::vector<double> &Bx_Original,
                             const std::vector<double> &By_Original,
                             const std::vector<double> &Bz_Original, int CTType,
                             std::vector<double> &Bx_Converted,
                             std::vector<double> &By_Converted,
                             std::vector<double> &Bz_Converted) {
  const std::size_t Num_Coordinate = Longitude.size();

  std::vector<double> Lamda_SCCS(Num_Coordinate), Theta_SCCS(Num_Coordinate);
  SCHA_CT(Longitude, Latitude, Lon0, Lat0, 1, Lamda_SCCS, Theta_SCCS);
  for (std::size_t i = 0; i < Num_Coordinate; ++i)
    Theta_SCCS[i] = 90.0 - Theta_SCCS[i];

  if (CTType == 1) {

    for (std::size_t i = 0; i < Num_Coordinate; ++i) {

      const double TempX1 =
          -sind(Latitude[i]) * cosd(Longitude[i]) * Bx_Original[i] -
          sind(Longitude[i]) * By_Original[i] -
          cosd(Latitude[i]) * cosd(Longitude[i]) * Bz_Original[i];

      const double TempY1 =
          -sind(Latitude[i]) * sind(Longitude[i]) * Bx_Original[i] +
          cosd(Longitude[i]) * By_Original[i] -
          cosd(Latitude[i]) * sind(Longitude[i]) * Bz_Original[i];

      const double TempZ1 = cosd(Latitude[i]) * Bx_Original[i] -
                            sind(Latitude[i]) * Bz_Original[i];

      const double TempX2 = cosd(Lon0) * TempX1 + sind(Lon0) * TempY1;
      const double TempY2 = -sind(Lon0) * TempX1 + cosd(Lon0) * TempY1;
      const double TempZ2 = TempZ1;

      const double TempX3 =
          cosd(90.0 - Lat0) * TempX2 - sind(90.0 - Lat0) * TempZ2;
      const double TempY3 = TempY2;
      const double TempZ3 =
          sind(90.0 - Lat0) * TempX2 + cosd(90.0 - Lat0) * TempZ2;

      Bx_Converted[i] = -sind(Theta_SCCS[i]) * cosd(Lamda_SCCS[i]) * TempX3 -
                        sind(Theta_SCCS[i]) * sind(Lamda_SCCS[i]) * TempY3 +
                        cosd(Theta_SCCS[i]) * TempZ3;

      By_Converted[i] =
          -sind(Lamda_SCCS[i]) * TempX3 + cosd(Lamda_SCCS[i]) * TempY3;

      Bz_Converted[i] = -cosd(Theta_SCCS[i]) * cosd(Lamda_SCCS[i]) * TempX3 -
                        cosd(Theta_SCCS[i]) * sind(Lamda_SCCS[i]) * TempY3 -
                        sind(Theta_SCCS[i]) * TempZ3;
    }
  } else if (CTType == 2) {

    for (std::size_t i = 0; i < Num_Coordinate; ++i) {
      const double TempX1 =
          -sind(Theta_SCCS[i]) * cosd(Lamda_SCCS[i]) * Bx_Original[i] -
          sind(Lamda_SCCS[i]) * By_Original[i] -
          cosd(Theta_SCCS[i]) * cosd(Lamda_SCCS[i]) * Bz_Original[i];

      const double TempY1 =
          -sind(Theta_SCCS[i]) * sind(Lamda_SCCS[i]) * Bx_Original[i] +
          cosd(Lamda_SCCS[i]) * By_Original[i] -
          cosd(Theta_SCCS[i]) * sind(Lamda_SCCS[i]) * Bz_Original[i];

      const double TempZ1 = cosd(Theta_SCCS[i]) * Bx_Original[i] -
                            sind(Theta_SCCS[i]) * Bz_Original[i];

      const double TempX2 =
          cosd(90.0 - Lat0) * TempX1 + sind(90.0 - Lat0) * TempZ1;
      const double TempY2 = TempY1;
      const double TempZ2 =
          -sind(90.0 - Lat0) * TempX1 + cosd(90.0 - Lat0) * TempZ1;

      const double TempX3 = cosd(Lon0) * TempX2 - sind(Lon0) * TempY2;
      const double TempY3 = sind(Lon0) * TempX2 + cosd(Lon0) * TempY2;
      const double TempZ3 = TempZ2;

      Bx_Converted[i] = -sind(Latitude[i]) * cosd(Longitude[i]) * TempX3 -
                        sind(Latitude[i]) * sind(Longitude[i]) * TempY3 +
                        cosd(Latitude[i]) * TempZ3;

      By_Converted[i] =
          -sind(Longitude[i]) * TempX3 + cosd(Longitude[i]) * TempY3;

      Bz_Converted[i] = -cosd(Latitude[i]) * cosd(Longitude[i]) * TempX3 -
                        cosd(Latitude[i]) * sind(Longitude[i]) * TempY3 -
                        sind(Latitude[i]) * TempZ3;
    }
  }
}

inline void Rotate2RankTensor(
    double Lon0, double Lat0, const std::vector<double> &Longitude,
    const std::vector<double> &Latitude, const std::vector<double> &Lamda,
    const std::vector<double> &Theta,

    const std::vector<double> &V_RSCHA, const std::vector<double> &Bx_RSCHA,
    const std::vector<double> &By_RSCHA, const std::vector<double> &Bz_RSCHA,
    const std::vector<double> &Txx_RSCHA, const std::vector<double> &Txy_RSCHA,
    const std::vector<double> &Tyy_RSCHA, const std::vector<double> &Tzx_RSCHA,
    const std::vector<double> &Tzy_RSCHA, const std::vector<double> &Tzz_RSCHA,

    std::vector<double> &V_LocalCart, std::vector<double> &Bx_LocalCart,
    std::vector<double> &By_LocalCart, std::vector<double> &Bz_LocalCart,
    std::vector<double> &Txx_LocalCart, std::vector<double> &Txy_LocalCart,
    std::vector<double> &Tyy_LocalCart, std::vector<double> &Tzx_LocalCart,
    std::vector<double> &Tzy_LocalCart, std::vector<double> &Tzz_LocalCart) {
  const std::size_t Num_Coordinate = Longitude.size();

  std::vector<double> Vx_RSCHA(Num_Coordinate), Vy_RSCHA(Num_Coordinate),
      Vz_RSCHA(Num_Coordinate);
  std::vector<double> Vxx_RSCHA(Num_Coordinate), Vxy_RSCHA(Num_Coordinate),
      Vyy_RSCHA(Num_Coordinate);
  std::vector<double> Vzx_RSCHA(Num_Coordinate), Vzy_RSCHA(Num_Coordinate),
      Vzz_RSCHA(Num_Coordinate);

  std::vector<double> VR_RSCHA(Num_Coordinate), VTheta_RSCHA(Num_Coordinate),
      VLamda_RSCHA(Num_Coordinate);
  std::vector<double> VRR_RSCHA(Num_Coordinate), VRTheta_RSCHA(Num_Coordinate),
      VRLamda_RSCHA(Num_Coordinate);
  std::vector<double> VThetaTheta_RSCHA(Num_Coordinate),
      VThetaLamda_RSCHA(Num_Coordinate), VLamdaLamda_RSCHA(Num_Coordinate);

  std::vector<double> Vx2(Num_Coordinate), Vy2(Num_Coordinate),
      Vz2(Num_Coordinate);
  std::vector<double> Vxx2(Num_Coordinate), Vxy2(Num_Coordinate),
      Vyy2(Num_Coordinate);
  std::vector<double> Vzx2(Num_Coordinate), Vzy2(Num_Coordinate),
      Vzz2(Num_Coordinate);

  std::vector<double> Vx1(Num_Coordinate), Vy1(Num_Coordinate),
      Vz1(Num_Coordinate);
  std::vector<double> Vxx1(Num_Coordinate), Vxy1(Num_Coordinate),
      Vyy1(Num_Coordinate);
  std::vector<double> Vzx1(Num_Coordinate), Vzy1(Num_Coordinate),
      Vzz1(Num_Coordinate);

  std::vector<double> Vx(Num_Coordinate), Vy(Num_Coordinate),
      Vz(Num_Coordinate);
  std::vector<double> Vxx(Num_Coordinate), Vxy(Num_Coordinate),
      Vyy(Num_Coordinate);
  std::vector<double> Vzx(Num_Coordinate), Vzy(Num_Coordinate),
      Vzz(Num_Coordinate);

  std::vector<double> VR_SCS(Num_Coordinate), VTheta_SCS(Num_Coordinate),
      VLamda_SCS(Num_Coordinate);
  std::vector<double> VRR_SCS(Num_Coordinate), VRTheta_SCS(Num_Coordinate),
      VRLamda_SCS(Num_Coordinate);
  std::vector<double> VThetaTheta_SCS(Num_Coordinate),
      VThetaLamda_SCS(Num_Coordinate), VLamdaLamda_SCS(Num_Coordinate);

  V_LocalCart = V_RSCHA;

  const double CoLat0 = 90.0 - Lat0;

  for (std::size_t i = 0; i < Num_Coordinate; ++i) {
    Vx_RSCHA[i] = -Bx_RSCHA[i];
    Vy_RSCHA[i] = -By_RSCHA[i];
    Vz_RSCHA[i] = -Bz_RSCHA[i];

    Vxx_RSCHA[i] = -Txx_RSCHA[i];
    Vxy_RSCHA[i] = -Txy_RSCHA[i];
    Vyy_RSCHA[i] = -Tyy_RSCHA[i];
    Vzx_RSCHA[i] = -Tzx_RSCHA[i];
    Vzy_RSCHA[i] = -Tzy_RSCHA[i];
    Vzz_RSCHA[i] = -Tzz_RSCHA[i];

    VR_RSCHA[i] = -Vz_RSCHA[i];
    VTheta_RSCHA[i] = -Vx_RSCHA[i];
    VLamda_RSCHA[i] = Vy_RSCHA[i];

    VRR_RSCHA[i] = Vzz_RSCHA[i];
    VRTheta_RSCHA[i] = Vzx_RSCHA[i];
    VRLamda_RSCHA[i] = -Vzy_RSCHA[i];
    VThetaTheta_RSCHA[i] = Vxx_RSCHA[i];
    VThetaLamda_RSCHA[i] = -Vxy_RSCHA[i];
    VLamdaLamda_RSCHA[i] = Vyy_RSCHA[i];
  }

  for (std::size_t i = 0; i < Num_Coordinate; ++i) {

    Vx2[i] = VR_RSCHA[i] * sind(Theta[i]) * cosd(Lamda[i]) +
             VTheta_RSCHA[i] * cosd(Theta[i]) * cosd(Lamda[i]) -
             VLamda_RSCHA[i] * sind(Lamda[i]);

    Vy2[i] = VR_RSCHA[i] * sind(Theta[i]) * sind(Lamda[i]) +
             VTheta_RSCHA[i] * cosd(Theta[i]) * sind(Lamda[i]) +
             VLamda_RSCHA[i] * cosd(Lamda[i]);

    Vz2[i] = VR_RSCHA[i] * cosd(Theta[i]) - VTheta_RSCHA[i] * sind(Theta[i]);

    Vxx2[i] =
        VRR_RSCHA[i] * std::pow(sind(Theta[i]), 2) *
            std::pow(cosd(Lamda[i]), 2) +
        VRTheta_RSCHA[i] * sind(2 * Theta[i]) * std::pow(cosd(Lamda[i]), 2) -
        VRLamda_RSCHA[i] * sind(Theta[i]) * sind(2 * Lamda[i]) +
        VThetaTheta_RSCHA[i] * std::pow(cosd(Theta[i]), 2) *
            std::pow(cosd(Lamda[i]), 2) -
        VThetaLamda_RSCHA[i] * cosd(Theta[i]) * sind(2 * Lamda[i]) +
        VLamdaLamda_RSCHA[i] * std::pow(sind(Lamda[i]), 2);

    Vxy2[i] =
        0.5 * VRR_RSCHA[i] * std::pow(sind(Theta[i]), 2) * sind(2 * Lamda[i]) +
        0.5 * VRTheta_RSCHA[i] * sind(2 * Theta[i]) * sind(2 * Lamda[i]) +
        VRLamda_RSCHA[i] * sind(Theta[i]) * cosd(2 * Lamda[i]) +
        0.5 * VThetaTheta_RSCHA[i] * std::pow(cosd(Theta[i]), 2) *
            sind(2 * Lamda[i]) +
        VThetaLamda_RSCHA[i] * cosd(Theta[i]) * cosd(2 * Lamda[i]) -
        0.5 * VLamdaLamda_RSCHA[i] * sind(2 * Lamda[i]);

    Vyy2[i] =
        VRR_RSCHA[i] * std::pow(sind(Theta[i]), 2) *
            std::pow(sind(Lamda[i]), 2) +
        VRTheta_RSCHA[i] * sind(2 * Theta[i]) * std::pow(sind(Lamda[i]), 2) +
        VRLamda_RSCHA[i] * sind(Theta[i]) * sind(2 * Lamda[i]) +
        VThetaTheta_RSCHA[i] * std::pow(cosd(Theta[i]), 2) *
            std::pow(sind(Lamda[i]), 2) +
        VThetaLamda_RSCHA[i] * cosd(Theta[i]) * sind(2 * Lamda[i]) +
        VLamdaLamda_RSCHA[i] * std::pow(cosd(Lamda[i]), 2);

    Vzx2[i] = 0.5 * VRR_RSCHA[i] * sind(2 * Theta[i]) * cosd(Lamda[i]) +
              VRTheta_RSCHA[i] * cosd(2 * Theta[i]) * cosd(Lamda[i]) -
              VRLamda_RSCHA[i] * cosd(Theta[i]) * sind(Lamda[i]) -
              0.5 * VThetaTheta_RSCHA[i] * sind(2 * Theta[i]) * cosd(Lamda[i]) +
              VThetaLamda_RSCHA[i] * sind(Theta[i]) * sind(Lamda[i]);

    Vzy2[i] = 0.5 * VRR_RSCHA[i] * sind(2 * Theta[i]) * sind(Lamda[i]) +
              VRTheta_RSCHA[i] * cosd(2 * Theta[i]) * sind(Lamda[i]) +
              VRLamda_RSCHA[i] * cosd(Theta[i]) * cosd(Lamda[i]) -
              0.5 * VThetaTheta_RSCHA[i] * sind(2 * Theta[i]) * sind(Lamda[i]) -
              VThetaLamda_RSCHA[i] * sind(Theta[i]) * cosd(Lamda[i]);

    Vzz2[i] = VRR_RSCHA[i] * std::pow(cosd(Theta[i]), 2) -
              VRTheta_RSCHA[i] * sind(2 * Theta[i]) +
              VThetaTheta_RSCHA[i] * std::pow(sind(Theta[i]), 2);
  }

  for (std::size_t i = 0; i < Num_Coordinate; ++i) {
    Vx1[i] = Vx2[i] * cosd(CoLat0) + Vz2[i] * sind(CoLat0);
    Vy1[i] = Vy2[i];
    Vz1[i] = -Vx2[i] * sind(CoLat0) + Vz2[i] * cosd(CoLat0);

    Vxx1[i] = Vxx2[i] * std::pow(cosd(CoLat0), 2) + Vzx2[i] * sind(2 * CoLat0) +
              Vzz2[i] * std::pow(sind(CoLat0), 2);
    Vxy1[i] = Vxy2[i] * cosd(CoLat0) + Vzy2[i] * sind(CoLat0);
    Vyy1[i] = Vyy2[i];
    Vzx1[i] = -0.5 * Vxx2[i] * sind(2 * CoLat0) + Vzx2[i] * cosd(2 * CoLat0) +
              0.5 * Vzz2[i] * sind(2 * CoLat0);
    Vzy1[i] = -Vxy2[i] * sind(CoLat0) + Vzy2[i] * cosd(CoLat0);
    Vzz1[i] = Vxx2[i] * std::pow(sind(CoLat0), 2) - Vzx2[i] * sind(2 * CoLat0) +
              Vzz2[i] * std::pow(cosd(CoLat0), 2);
  }

  for (std::size_t i = 0; i < Num_Coordinate; ++i) {
    Vx[i] = Vx1[i] * cosd(Lon0) - Vy1[i] * sind(Lon0);
    Vy[i] = Vx1[i] * sind(Lon0) + Vy1[i] * cosd(Lon0);
    Vz[i] = Vz1[i];

    Vxx[i] = Vxx1[i] * std::pow(cosd(Lon0), 2) - Vxy1[i] * sind(2 * Lon0) +
             Vyy1[i] * std::pow(sind(Lon0), 2);
    Vxy[i] = 0.5 * Vxx1[i] * sind(2 * Lon0) + Vxy1[i] * cosd(2 * Lon0) -
             0.5 * Vyy1[i] * sind(2 * Lon0);
    Vyy[i] = Vxx1[i] * std::pow(sind(Lon0), 2) + Vxy1[i] * sind(2 * Lon0) +
             Vyy1[i] * std::pow(cosd(Lon0), 2);
    Vzx[i] = Vzx1[i] * cosd(Lon0) - Vzy1[i] * sind(Lon0);
    Vzy[i] = Vzx1[i] * sind(Lon0) + Vzy1[i] * cosd(Lon0);
    Vzz[i] = Vzz1[i];
  }

  for (std::size_t i = 0; i < Num_Coordinate; ++i) {
    const double colat = 90.0 - Latitude[i];

    VR_SCS[i] = Vx[i] * sind(colat) * cosd(Longitude[i]) +
                Vy[i] * sind(colat) * sind(Longitude[i]) + Vz[i] * cosd(colat);

    VTheta_SCS[i] = Vx[i] * cosd(colat) * cosd(Longitude[i]) +
                    Vy[i] * cosd(colat) * sind(Longitude[i]) -
                    Vz[i] * sind(colat);

    VLamda_SCS[i] = -Vx[i] * sind(Longitude[i]) + Vy[i] * cosd(Longitude[i]);

    VRR_SCS[i] =
        Vxx[i] * std::pow(sind(colat), 2) * std::pow(cosd(Longitude[i]), 2) +
        Vxy[i] * std::pow(sind(colat), 2) * sind(2 * Longitude[i]) +
        Vzx[i] * sind(2 * colat) * cosd(Longitude[i]) +
        Vyy[i] * std::pow(sind(colat), 2) * std::pow(sind(Longitude[i]), 2) +
        Vzy[i] * sind(2 * colat) * sind(Longitude[i]) +
        Vzz[i] * std::pow(cosd(colat), 2);

    VRTheta_SCS[i] =
        0.5 * Vxx[i] * sind(2 * colat) * std::pow(cosd(Longitude[i]), 2) +
        0.5 * Vxy[i] * sind(2 * colat) * sind(2 * Longitude[i]) +
        Vzx[i] * cosd(2 * colat) * cosd(Longitude[i]) +
        0.5 * Vyy[i] * sind(2 * colat) * std::pow(sind(Longitude[i]), 2) +
        Vzy[i] * cosd(2 * colat) * sind(Longitude[i]) -
        0.5 * Vzz[i] * sind(2 * colat);

    VRLamda_SCS[i] = -0.5 * Vxx[i] * sind(colat) * sind(2 * Longitude[i]) +
                     Vxy[i] * sind(colat) * cosd(2 * Longitude[i]) -
                     Vzx[i] * cosd(colat) * sind(Longitude[i]) +
                     0.5 * Vyy[i] * sind(colat) * sind(2 * Longitude[i]) +
                     Vzy[i] * cosd(colat) * cosd(Longitude[i]);

    VThetaTheta_SCS[i] =
        Vxx[i] * std::pow(cosd(colat), 2) * std::pow(cosd(Longitude[i]), 2) +
        Vxy[i] * std::pow(cosd(colat), 2) * sind(2 * Longitude[i]) -
        Vzx[i] * sind(2 * colat) * cosd(Longitude[i]) +
        Vyy[i] * std::pow(cosd(colat), 2) * std::pow(sind(Longitude[i]), 2) -
        Vzy[i] * sind(2 * colat) * sind(Longitude[i]) +
        Vzz[i] * std::pow(sind(colat), 2);

    VThetaLamda_SCS[i] = -0.5 * Vxx[i] * cosd(colat) * sind(2 * Longitude[i]) +
                         Vxy[i] * cosd(colat) * cosd(2 * Longitude[i]) +
                         Vzx[i] * sind(colat) * sind(Longitude[i]) +
                         0.5 * Vyy[i] * cosd(colat) * sind(2 * Longitude[i]) -
                         Vzy[i] * sind(colat) * cosd(Longitude[i]);

    VLamdaLamda_SCS[i] = Vxx[i] * std::pow(sind(Longitude[i]), 2) -
                         Vxy[i] * sind(2 * Longitude[i]) +
                         Vyy[i] * std::pow(cosd(Longitude[i]), 2);
  }

  std::vector<double> Vx_LocalCart(Num_Coordinate),
      Vy_LocalCart(Num_Coordinate), Vz_LocalCart(Num_Coordinate);
  std::vector<double> Vxx_LocalCart(Num_Coordinate),
      Vxy_LocalCart(Num_Coordinate), Vyy_LocalCart(Num_Coordinate);
  std::vector<double> Vzx_LocalCart(Num_Coordinate),
      Vzy_LocalCart(Num_Coordinate), Vzz_LocalCart(Num_Coordinate);

  for (std::size_t i = 0; i < Num_Coordinate; ++i) {
    Vx_LocalCart[i] = -VTheta_SCS[i];
    Vy_LocalCart[i] = VLamda_SCS[i];
    Vz_LocalCart[i] = -VR_SCS[i];

    Vxx_LocalCart[i] = VThetaTheta_SCS[i];
    Vxy_LocalCart[i] = -VThetaLamda_SCS[i];
    Vyy_LocalCart[i] = VLamdaLamda_SCS[i];
    Vzx_LocalCart[i] = VRTheta_SCS[i];
    Vzy_LocalCart[i] = -VRLamda_SCS[i];
    Vzz_LocalCart[i] = VRR_SCS[i];
  }

  for (std::size_t i = 0; i < Num_Coordinate; ++i) {
    Bx_LocalCart[i] = -Vx_LocalCart[i];
    By_LocalCart[i] = -Vy_LocalCart[i];
    Bz_LocalCart[i] = -Vz_LocalCart[i];

    Txx_LocalCart[i] = -Vxx_LocalCart[i];
    Txy_LocalCart[i] = -Vxy_LocalCart[i];
    Tyy_LocalCart[i] = -Vyy_LocalCart[i];
    Tzx_LocalCart[i] = -Vzx_LocalCart[i];
    Tzy_LocalCart[i] = -Vzy_LocalCart[i];
    Tzz_LocalCart[i] = -Vzz_LocalCart[i];
  }
}

}
