#pragma once
#include "BoundaryProjection.h"
#include "SHRotation.h"
#include "Execution.h"
#include "../R_SCHA.h"

namespace SHRSCHA {
// Work3 needs only the coefficient map. Avoid assembling the cap SH Gram
// and reverse moments, which can be expensive at high global degree.
inline Matrix analytic_forward_block(const Geometry& geometry,int m,int N,
    const std::vector<double>& roots,int internal_count,int external_count,int P,int n_min=1) {
  geometry.validate();
  const int first=std::max(m,n_min),L=std::max(0,N-first+1),K=std::max(internal_count,external_count);
  if(m<0 || internal_count<0 || external_count<0 || P<0 || n_min<0 || n_min>1 ||
      roots.size()<static_cast<std::size_t>(K))throw std::invalid_argument("Invalid analytical forward block");
  Matrix M(internal_count+external_count+P+1,2*L);
  if(!L)return M;
  const auto side=sh_angular(m,N,geometry.theta());
  for(int k=0;k<K;++k) {
    const double nu=roots[k];
    const double derivative=Legendre::BoundaryDegreeDerivativeLegacySchmidt(m,nu,geometry.theta());
    for(int j=0;j<L;++j) {
      const int n=first+j;const double gap=n*double(n+1)-nu*(nu+1);
      double overlap=0;
      if(std::fabs(gap)>1e-7*std::max({1.0,n*double(n+1),nu*(nu+1)})) {
        if(!std::isfinite(derivative) || derivative==0)
          throw std::runtime_error("Boundary degree derivative is unresolved");
        overlap=(2*nu+1)/gap*side.value[n]/derivative;
      }else {
        const auto q=gauss_interval(std::cos(geometry.theta()),1,std::max(128,2*(N+K)+32));
        long double cross=0,norm2=0;
        for(std::size_t point=0;point<q.x.size();++point) {
          const double theta=std::acos(q.x[point]);
          const auto rs=Legendre::EvaluateLegacySchmidt(m,nu,theta);
          if(rs.status!=0)throw std::runtime_error("Coincident angular overlap failed");
          const double value=static_cast<double>(rs.value);
          cross+=q.w[point]*sh_angular(m,N,theta).value[n]*value;
          norm2+=q.w[point]*value*value;
        }
        if(!(norm2>0))throw std::runtime_error("Coincident angular norm is unresolved");
        overlap=static_cast<double>(cross/norm2);
      }
      const auto radial=legendre_radial(geometry,n,nu);
      if(k<internal_count) {M(k,j)=radial[0]*overlap;M(k,L+j)=radial[1]*overlap;}
      if(k<external_count) {M(internal_count+k,j)=radial[2]*overlap;M(internal_count+k,L+j)=radial[3]*overlap;}
    }
  }
  for(int p=1;p<=P+1;++p) {
    const int radial_p=p<=P?p:0,row=internal_count+external_count+p-1;
    double value=0,derivative=0;
    if(radial_p)Mehler::Mehler_1rank(m,radial_p*std::acos(-1.0)/geometry.shell_log(),
      geometry.aperture_deg,geometry.aperture_deg,value,derivative);
    else Constant::Constant_1rank(m,geometry.aperture_deg,geometry.aperture_deg,value,derivative);
    if(!std::isfinite(value) || value==0)throw std::runtime_error("Side basis normalization is unresolved");
    for(int j=0;j<L;++j) {
      const int n=first+j;const double factor=geometry.reference/geometry.lower*side.value[n]/value;
      M(row,j)=factor*radial_moment(geometry,n,radial_p,false);
      M(row,L+j)=factor*radial_moment(geometry,n,radial_p,true);
    }
  }
  finite(M);return M;
}
struct SHCoefficients {
  int degree=0;
  std::vector<double> g,h;
  explicit SHCoefficients(int N=0):degree(N),g(SHRotation::triangular_size(N)),h(g.size()){}
};
struct RegionalModel {
  double longitude=0,latitude=90;
  Geometry geometry;
  int internal=0,external=0,mehler=0,constant=0;
  // Optional active side truncation; legacy storage still uses a square Mehler block.
  int side_radial_limit=-1,side_order_limit=-1;
  std::vector<double> coefficients,roots;
  int count()const{return R_SCHA::coefficient_count(internal,external,mehler,constant);}
};
inline RegionalModel project_model(RegionalModel model,SHCoefficients internal,
                                   SHCoefficients external,bool geographic=true) {
  model.geometry.validate();
  for(const auto* source:{&internal,&external}) {
    if(source->g.size()!=SHRotation::triangular_size(source->degree) || source->h.size()!=source->g.size())
      throw std::invalid_argument("SH input must contain a complete triangular coefficient table");
    for(std::size_t i=0;i<source->g.size();++i)
      if(!std::isfinite(source->g[i]) || !std::isfinite(source->h[i]))
        throw std::invalid_argument("SH input must be finite");
  }
  if(model.internal>2048 || model.external>2048 || model.mehler>1024 || model.constant>2048)
    throw std::invalid_argument("Regional truncation exceeds this implementation's guard");
  if(model.side_radial_limit>model.mehler || model.side_order_limit>model.mehler ||
      model.side_radial_limit < -1 || model.side_order_limit < -1)
    throw std::invalid_argument("Active side limits exceed legacy Mehler storage");
  const int N=std::max(internal.degree,external.degree);
  internal.g.resize(SHRotation::triangular_size(N));internal.h.resize(internal.g.size());
  external.g.resize(internal.g.size());external.h.resize(internal.g.size());
  internal.g[0]=internal.h[0]=external.g[0]=external.h[0]=0; // magnetic n_min=1
  if(geographic) {
    SHRotation::geographic_to_cap(model.longitude,model.latitude,N,internal.g,internal.h);
    SHRotation::geographic_to_cap(model.longitude,model.latitude,N,external.g,external.h);
  }
  const int K=std::max(model.internal,model.external),ti=(model.internal+1)*(model.internal+2)/2;
  const int te=(model.external+1)*(model.external+2)/2,tm=model.mehler*(model.mehler+1);
  model.coefficients.assign(model.count(),0);
  const auto root_count=static_cast<std::size_t>((K+1)*(K+2)/2);
  // A supplied root cache follows the legacy API contract: it must correspond
  // to this aperture and truncation. File-based callers leave it empty.
  if(model.roots.size()!=root_count) {
    model.roots.resize(root_count);
    R_SCHA::R_SCHA_nk(K,model.geometry.aperture_deg,model.roots);
  }
  for(double root:model.roots)if(!std::isfinite(root) || root<=0)
    throw std::invalid_argument("Invalid regional root cache");
  const int maximum=std::max({model.internal,model.external,model.mehler,model.constant});
  for(int m=execution_rank();m<=maximum;m+=execution_size()) {
    const int ki=std::max(0,model.internal-m+1),ke=std::max(0,model.external-m+1);
    const int radial=model.side_radial_limit>=0?model.side_radial_limit:model.mehler;
    const int side_order=model.side_order_limit>=0?model.side_order_limit:model.mehler;
    const int P=m<=side_order?radial:0;
    std::vector<double> roots;
    for(int k=m;k<=K;++k)roots.push_back(model.roots[k*(k+1)/2+m]);
    const auto M=analytic_forward_block(model.geometry,m,N,roots,ki,ke,P);
    const int first=std::max(1,m),L=std::max(0,N-first+1);
    std::vector<double> xc(2*L),xs(2*L);
    for(int j=0;j<L;++j) {
      const auto ix=SHRotation::index(first+j,m);
      xc[j]=internal.g[ix];xs[j]=internal.h[ix];xc[L+j]=external.g[ix];xs[L+j]=external.h[ix];
    }
    const auto rc=multiply(M,xc),rs=multiply(M,xs);
    for(int k=m;k<=model.internal;++k) {
      const int ix=k*(k+1)/2+m;model.coefficients[ix]=rc[k-m];model.coefficients[ti+ix]=rs[k-m];
    }
    for(int k=m;k<=model.external;++k) {
      const int ix=k*(k+1)/2+m;model.coefficients[2*ti+ix]=rc[ki+k-m];model.coefficients[2*ti+te+ix]=rs[ki+k-m];
    }
    for(int p=1;p<=P;++p) {
      const int ix=(p-1)*(model.mehler+1)+m;
      model.coefficients[2*(ti+te)+ix]=rc[ki+ke+p-1];
      model.coefficients[2*(ti+te)+tm+ix]=rs[ki+ke+p-1];
    }
    if(m<=model.constant) {
      model.coefficients[2*(ti+te+tm)+m]=rc.back();
      model.coefficients[2*(ti+te+tm)+model.constant+1+m]=rs.back();
    }
  }
  combine_coefficients(model.coefficients);
  return model;
}
} // namespace SHRSCHA
