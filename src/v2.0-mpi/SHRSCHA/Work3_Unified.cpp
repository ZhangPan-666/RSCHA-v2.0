#include "ModelIO.h"
#include <iostream>

#ifdef _WIN32
using ArgumentChar=wchar_t;
#else
using ArgumentChar=char;
#endif
int run(int argc,ArgumentChar** argv) {
  if(argc>3)throw std::invalid_argument("Usage: Work3_Unified [input_directory] [output_directory]");
  const auto input=argc>1?RSCHACompat::path(argv[1]):RSCHACompat::current_path();
  const auto output=argc>2?RSCHACompat::path(argv[2]):input;
  auto model=SHRSCHA::read_regional(input/"Parameters_Modelling.dat",false);
  const auto internal=SHRSCHA::read_sh(input/"Parameters_SHModel.dat");
  const auto external_path=input/"Parameters_SHModel_External.dat";
  const bool has_external=RSCHACompat::exists(external_path);
  const auto external=has_external?SHRSCHA::read_sh(external_path):SHRSCHA::SHCoefficients(internal.degree);
  const auto side_path=input/"Parameters_SideProjection.dat";
  const bool explicit_side=RSCHACompat::exists(side_path);
  if(explicit_side) {
    const auto side=SHRSCHA::read_table(side_path);
    if(side.size()!=2 || side[0].size()!=1 || side[1].size()!=1)
      throw std::invalid_argument("Side projection requires two rows: radial_p_count, azimuth_m_max");
    const int radial=R_SCHA::checked_index(side[0][0]),order=R_SCHA::checked_index(side[1][0]);
    if(radial>1024 || order>1024)throw std::invalid_argument("Side projection exceeds storage guard");
    model.side_radial_limit=radial;model.side_order_limit=order;
    model.mehler=std::max(radial,order);model.constant=std::max(model.constant,order);
  }
  model=SHRSCHA::project_model(model,internal,external);
  std::vector<double> degree,spectrum;
  R_SCHA::RSCHA_Spectrum(model.coefficients,model.internal,model.external,model.mehler,model.constant,
      model.roots,model.geometry.aperture_deg,model.geometry.lower,model.geometry.lower,degree,spectrum);
  if(SHRSCHA::execution_rank()!=0)return 0;
  RSCHACompat::create_directories(output);
  SHRSCHA::write_regional(output/"Results_RSCHAModel.dat",model);
  std::ofstream file((output/"Results_RSCHASpectrum.dat").c_str());if(!file)throw std::runtime_error("Cannot write spectrum");
  file<<std::setprecision(17);for(std::size_t i=0;i<degree.size();++i)file<<degree[i]<<','<<spectrum[i]<<'\n';
  file.close();if(!file)throw std::runtime_error("Spectrum write failed");
  std::ofstream meta((output/"Projection_Metadata.txt").c_str());if(!meta)throw std::runtime_error("Cannot write projection metadata");
  meta<<"method=unified analytical boundary projection\nrotation=Wigner passive geographic-to-cap\n"
      <<"SH_n_min=1\nSH_degree="<<std::max(internal.degree,external.degree)
      <<"\nexternal_SH_input="<<(has_external?"present":"zero")
      <<"\nside_radius_factor=R_E/a\nR_E="<<std::setprecision(17)<<model.geometry.reference
      <<"\nspectrum_scope=internal and external Legendre families only\n"
      <<"side_configuration="<<(explicit_side?"independent active limits, legacy padded storage":"legacy coupled limits")
      <<"\nside_radial_p="<<(explicit_side?model.side_radial_limit:model.mehler)
      <<"\nside_azimuth_m="<<(explicit_side?model.side_order_limit:model.mehler)
      <<"\nside_storage_P="<<model.mehler<<"\nconstant_m="<<model.constant
      <<"\nside_high_orders_truncated="<<(std::max(model.constant,explicit_side?model.side_order_limit:model.mehler)<std::max(internal.degree,external.degree)?"yes":"no")<<'\n';
  meta.close();if(!meta)throw std::runtime_error("Metadata write failed");
  std::cout<<"Unified SH -> R-SCHA projection complete: "<<model.coefficients.size()<<" coefficients.\n";
  return 0;
}
#ifdef _WIN32
int wmain(int argc,wchar_t** argv){return SHRSCHA::guarded_application([&](){return run(argc,argv);});}
#else
int main(int argc,char** argv){return SHRSCHA::guarded_application([&](){return run(argc,argv);});}
#endif
