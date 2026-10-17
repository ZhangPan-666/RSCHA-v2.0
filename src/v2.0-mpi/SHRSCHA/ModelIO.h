#pragma once
#include "AnalyticForward.h"
#include "../RSCHACompatibility.h"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace SHRSCHA {
inline std::vector<std::vector<double>> read_table(const RSCHACompat::path& path) {
  std::ifstream input(path.c_str());if(!input)throw std::runtime_error("Cannot read model input");
  std::vector<std::vector<double>> table;std::string line;
  while(std::getline(input,line)) {
    std::replace(line.begin(),line.end(),',',' ');std::istringstream stream(line);
    std::vector<double> row;double value=0;
    while(stream>>value) {if(!std::isfinite(value))throw std::runtime_error("Non-finite model input");row.push_back(value);}
    if(!stream.eof())throw std::runtime_error("Malformed model input");
    if(!row.empty())table.push_back(row);
  }
  return table;
}
inline SHCoefficients read_sh(const RSCHACompat::path& path) {
  const auto table=read_table(path);
  if(table.empty() || table.back().size()!=4)throw std::runtime_error("SH input must have four columns");
  const int N=R_SCHA::checked_index(table.back()[0]);SHCoefficients result(N);
  if(table.size()!=result.g.size())throw std::runtime_error("SH table must include the complete n=0..N triangle");
  for(int n=0;n<=N;++n)for(int m=0;m<=n;++m) {
    const auto ix=SHRotation::index(n,m);
    if(table[ix].size()!=4 || table[ix][0]!=n || table[ix][1]!=m)
      throw std::runtime_error("SH rows must follow degree/order sequence");
    result.g[ix]=table[ix][2];result.h[ix]=m?table[ix][3]:0;
  }
  return result;
}
inline RegionalModel read_regional(const RSCHACompat::path& path,bool coefficients=true) {
  const auto table=read_table(path);
  if(table.size()<9)throw std::runtime_error("Regional input needs nine geometry/truncation rows");
  for(const auto& row:table)if(row.size()!=1)throw std::runtime_error("Regional input must contain one column");
  RegionalModel model;model.longitude=table[0][0];model.latitude=table[1][0];
  model.geometry.aperture_deg=table[2][0];model.internal=R_SCHA::checked_index(table[3][0]);
  model.external=R_SCHA::checked_index(table[4][0]);model.mehler=R_SCHA::checked_index(table[5][0]);
  model.constant=R_SCHA::checked_index(table[6][0]);model.geometry.upper=table[7][0];model.geometry.lower=table[8][0];
  model.geometry.validate();R_SCHA::check_geometry(model.longitude,model.latitude,model.geometry.aperture_deg,
                                               model.geometry.lower,model.geometry.upper);
  if(coefficients) {
    if(table.size()!=static_cast<std::size_t>(9+model.count()))throw std::runtime_error("Invalid regional coefficient count");
    for(std::size_t i=9;i<table.size();++i)model.coefficients.push_back(table[i][0]);
  }else if(table.size()!=9)throw std::runtime_error("Modelling parameters must contain exactly nine rows");
  return model;
}
inline void write_regional(const RSCHACompat::path& path,const RegionalModel& model) {
  std::ofstream output(path.c_str());if(!output)throw std::runtime_error("Cannot write regional model");
  output<<std::setprecision(17)<<model.longitude<<'\n'<<model.latitude<<'\n'<<model.geometry.aperture_deg<<'\n'
        <<model.internal<<'\n'<<model.external<<'\n'<<model.mehler<<'\n'<<model.constant<<'\n'
        <<model.geometry.upper<<'\n'<<model.geometry.lower<<'\n';
  for(double value:model.coefficients)output<<value<<'\n';
  output.close();if(!output)throw std::runtime_error("Regional model write failed");
}
inline void write_sh(const RSCHACompat::path& path,const SHCoefficients& sh) {
  std::ofstream output(path.c_str());if(!output)throw std::runtime_error("Cannot write SH model");
  output<<std::setprecision(17);
  for(int n=0;n<=sh.degree;++n)for(int m=0;m<=n;++m) {
    const auto ix=SHRotation::index(n,m);output<<n<<','<<m<<','<<sh.g[ix]<<','<<sh.h[ix]<<'\n';
  }
  output.close();if(!output)throw std::runtime_error("SH model write failed");
}
} // namespace SHRSCHA
