#pragma once
#include "SHRSCHA/AnalyticForward.h"

namespace namespace_ConvertSH2RSCHA {
inline void RotateSH(double Lon0, double Lat0, double RE, std::vector<int>& n,
                     std::vector<int>& m, std::vector<double>& gnm,
                     std::vector<double>& hnm, int n_Start, int n_End) {
  if(n_Start!=1 || !(RE>0) || !std::isfinite(RE))
    throw std::invalid_argument("Magnetic SH rotation requires n_Start=1 and positive reference radius");
  const auto count=SHRotation::triangular_size(n_End);
  if(n.size()!=count || m.size()!=count || gnm.size()!=count || hnm.size()!=count)
    throw std::invalid_argument("SH input must be a complete triangular table");
  for(int degree=0;degree<=n_End;++degree)for(int order=0;order<=degree;++order) {
    const auto ix=SHRotation::index(degree,order);
    if(n[ix]!=degree || m[ix]!=order)throw std::invalid_argument("Invalid SH degree/order table");
  }
  gnm[0]=hnm[0]=0;
  SHRotation::geographic_to_cap(Lon0,Lat0,n_End,gnm,hnm);
}

inline void ConvertSH2RSCHA(double Lon0, double Lat0, double Theta0,
                            int K_Max_Internal, int K_Max_External, int P_Max,
                            int M_Max, double Min_Radius, double Max_Radius,
                            double RE, std::vector<double>& nk,
                            std::vector<int>& n, std::vector<int>& m,
                            std::vector<double>& gnm, std::vector<double>& hnm,
                            int n_Start, int n_End,
                            std::vector<double>& Coefficients) {
  RotateSH(Lon0,Lat0,RE,n,m,gnm,hnm,n_Start,n_End);
  SHRSCHA::RegionalModel model;
  model.longitude=Lon0;model.latitude=Lat0;
  model.geometry.reference=RE;model.geometry.lower=Min_Radius;
  model.geometry.upper=Max_Radius;model.geometry.aperture_deg=Theta0;
  model.internal=K_Max_Internal;model.external=K_Max_External;
  model.mehler=P_Max;model.constant=M_Max;model.roots=nk;
  SHRSCHA::SHCoefficients internal(n_End),external(n_End);
  internal.g=gnm;internal.h=hnm;
  model=SHRSCHA::project_model(model,internal,external,false);
  nk=std::move(model.roots);Coefficients=std::move(model.coefficients);
}
} // namespace namespace_ConvertSH2RSCHA
