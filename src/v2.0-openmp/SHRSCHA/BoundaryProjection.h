#pragma once
#include "Numerics.h"
#include "../Basis Function.h"
#include <array>
#include <string>

namespace SHRSCHA {
struct Geometry {
  double reference=6371200,lower=6356200,upper=6871200,aperture_deg=27;
  double theta() const {return aperture_deg*std::acos(-1.0)/180;}
  double shell_log() const {return std::log1p((upper-lower)/lower);}
  void validate() const {
    for(double v:{reference,lower,upper,aperture_deg})
      if(!std::isfinite(v))throw std::invalid_argument("Geometry must be finite");
    if(reference<=0 || lower<=0 || upper<=lower || aperture_deg<=0 || aperture_deg>90)
      throw std::invalid_argument("Invalid shell geometry");
  }
};
inline double exp_difference(double x,double y) {
  return x>=y?std::exp(x)*(-std::expm1(y-x)):std::exp(y)*std::expm1(x-y);
}
// Entries ii, ie, ei, ee. These map both SH source families into both
// regional Legendre families, before multiplication by the angular overlap.
inline std::array<double,4> legendre_radial(const Geometry& g,double n,double nu) {
  const double s=g.shell_log(),t=g.reference/g.lower;
  const double denominator=-std::expm1(-(2*nu+1)*s);
  if(!(nu>0) || !(denominator>0))throw std::invalid_argument("Invalid regional degree");
  const double internal=std::pow(t,n+2),external=std::pow(1/t,n-1);
  return {
    (n+1)/(nu+1)*(-std::expm1(-(n+nu+1)*s))/denominator*internal,
    n/(nu+1)*std::expm1((n-nu)*s)/denominator*external,
    (n+1)/nu*exp_difference(-(2*nu+1)*s,-(n+nu+1)*s)/denominator*internal,
    n/nu*exp_difference((n-nu)*s,-(2*nu+1)*s)/denominator*external};
}
inline double radial_basis(const Geometry& g,int p,double radius) {
  const double s=g.shell_log();
  if(p==0)return 1/std::sqrt(std::expm1(s));
  const double A=2*p*std::acos(-1.0),angle=p*std::acos(-1.0)*std::log(radius/g.lower)/s;
  return std::sqrt(2/(s*(A*A+s*s)))*std::sqrt(g.lower/radius)*
         (A*std::cos(angle)+s*std::sin(angle));
}
// rho = (1/a) integral_a^b SH_radial(r) R_p(r) dr.
// p=0 denotes the constant radial sequence. Both SH radial powers are included.
inline double radial_moment(const Geometry& g,int n,int p,bool external) {
  const double s=g.shell_log(),t=g.reference/g.lower;
  if(p==0) {
    const double r0=1/std::sqrt(std::expm1(s));
    if(external)return r0*std::pow(1/t,n)*std::expm1((n+1)*s)/(n+1);
    return r0*std::pow(t,n+1)*(n==0?s:-std::expm1(-n*s)/n);
  }
  const double pi=std::acos(-1.0),omega=p*pi/s;
  const double q=external?n+.5:-(n+.5);
  const double T=p%2?-std::exp(q*s)-1:std::expm1(q*s);
  const double scale=std::sqrt(2/(s*(4*p*double(p)*pi*pi+s*s)));
  return scale*p*pi*(2*q-1)*T/(q*q+omega*omega)*
         (external?std::pow(1/t,n):std::pow(t,n+1));
}
struct Mode {std::string family;int sequence=0;double degree=0;};
struct BoundaryOperators {
  Geometry geometry;
  int order=0,degree=0,first_degree=0;
  int quadrature_count=0,coincident_overlap_count=0;
  std::vector<double> roots,norms;
  std::vector<Mode> modes;
  Matrix forward,boundary,regional_boundary,gamma_integral,gamma_analytic;
  // x = [SH-internal n_first..N, SH-external n_first..N], one azimuthal channel.
  // r = [R-internal roots, R-external roots, Mehler p=1..P, constant radial mode].
  // B rows = [lower SH moments, upper SH moments, side p=1..P, side p=0].
};
inline BoundaryOperators assemble_boundary(const Geometry& geometry,int m,int N,
    int internal_count,int external_count,int mehler_count,int n_min=1,
    int quadrature_count=0,bool analytical_overlap=true) {
  geometry.validate();
  if(m<0 || N<std::max(m,n_min) || n_min<0 || n_min>1 ||
      internal_count<0 || external_count<0 || mehler_count<0 ||
      internal_count>2048 || external_count>2048 || mehler_count>1024)
    throw std::invalid_argument("Invalid boundary truncations");
  BoundaryOperators op;
  op.geometry=geometry;op.order=m;op.degree=N;op.first_degree=std::max(m,n_min);
  const int L=N-op.first_degree+1,K=std::max(internal_count,external_count);
  const int R=internal_count+external_count+mehler_count+1;
  const int nq=quadrature_count?quadrature_count:std::max(128,2*(N+K)+32);
  if(nq<2)throw std::invalid_argument("Invalid angular quadrature count");
  op.quadrature_count=nq;
  op.forward=Matrix(R,2*L);op.boundary=Matrix(2*L+mehler_count+1,2*L);
  op.regional_boundary=Matrix(op.boundary.rows,R);
  op.gamma_integral=Matrix(K,L);op.gamma_analytic=Matrix(K,L);
  op.roots.resize(K);op.norms.resize(K);
  for(int k=0;k<K;++k) {
    op.roots[k]=Legendre::FindDirichletDegree(m,k+1,geometry.theta());
    if(!std::isfinite(op.roots[k]))throw std::runtime_error("Regional root calculation failed");
  }
  const auto q=gauss_interval(std::cos(geometry.theta()),1,nq);
  Matrix S(nq,L),Rbasis(nq,K),C(L,L),X(L,K);
  for(int point=0;point<nq;++point) {
    const double theta=std::acos(q.x[point]);
    const auto angular=sh_angular(m,N,theta);
    for(int j=0;j<L;++j)S(point,j)=angular.value[j+op.first_degree];
    for(int k=0;k<K;++k) {
      const auto value=Legendre::EvaluateLegacySchmidt(m,op.roots[k],theta);
      if(value.status!=0)throw std::runtime_error("Regional angular evaluation failed");
      Rbasis(point,k)=static_cast<double>(value.value);
    }
  }
  for(int i=0;i<L;++i)for(int j=0;j<L;++j) {
    long double sum=0;for(int point=0;point<nq;++point)sum+=q.w[point]*S(point,i)*S(point,j);
    C(i,j)=static_cast<double>(sum);
  }
  for(int k=0;k<K;++k) {
    long double sum=0;for(int point=0;point<nq;++point)sum+=q.w[point]*Rbasis(point,k)*Rbasis(point,k);
    op.norms[k]=static_cast<double>(sum);
    if(!(op.norms[k]>0))throw std::runtime_error("Regional angular norm is unresolved");
    for(int j=0;j<L;++j) {
      sum=0;for(int point=0;point<nq;++point)sum+=q.w[point]*S(point,j)*Rbasis(point,k);
      X(j,k)=static_cast<double>(sum);
      op.gamma_integral(k,j)=X(j,k)/op.norms[k];
    }
  }
  const auto side=sh_angular(m,N,geometry.theta());
  for(int k=0;k<K;++k) {
    const double nu=op.roots[k];
    const double derivative=Legendre::BoundaryDegreeDerivativeLegacySchmidt(m,nu,geometry.theta());
    for(int j=0;j<L;++j) {
      const int n=j+op.first_degree;
      const double gap=n*double(n+1)-nu*(nu+1);
      double gamma=op.gamma_integral(k,j);
      if(std::fabs(gap)>1e-7*std::max({1.0,n*double(n+1),nu*(nu+1)})) {
        if(!std::isfinite(derivative) || derivative==0)
          throw std::runtime_error("Boundary degree derivative is unresolved");
        gamma=(2*nu+1)/gap*side.value[n]/derivative;
      }else ++op.coincident_overlap_count;
      op.gamma_analytic(k,j)=gamma;
      const auto radial=legendre_radial(geometry,n,nu);
      const double overlap=analytical_overlap?gamma:op.gamma_integral(k,j);
      if(k<internal_count) {
        op.forward(k,j)=radial[0]*overlap;
        op.forward(k,L+j)=radial[1]*overlap;
      }
      if(k<external_count) {
        op.forward(internal_count+k,j)=radial[2]*overlap;
        op.forward(internal_count+k,L+j)=radial[3]*overlap;
      }
    }
  }
  for(int k=0;k<internal_count;++k)op.modes.push_back({"internal",k+1,op.roots[k]});
  for(int k=0;k<external_count;++k)op.modes.push_back({"external",k+1,op.roots[k]});
  const double s=geometry.shell_log(),t=geometry.reference/geometry.lower;
  for(int i=0;i<L;++i) {
    for(int j=0;j<L;++j) {
      const int n=j+op.first_degree;
      op.boundary(i,j)=C(i,j)*(n+1)*std::pow(t,n+2);
      op.boundary(i,L+j)=C(i,j)*(-n)*std::pow(1/t,n-1);
      op.boundary(L+i,j)=C(i,j)*(-(n+1))*std::pow(geometry.reference/geometry.upper,n+2);
      op.boundary(L+i,L+j)=C(i,j)*n*std::pow(geometry.upper/geometry.reference,n-1);
    }
    for(int k=0;k<K;++k) {
      const double nu=op.roots[k];
      if(k<internal_count) {
        op.regional_boundary(i,k)=X(i,k)*(nu+1);
        op.regional_boundary(L+i,k)=X(i,k)*(-(nu+1))*std::exp(-(nu+2)*s);
      }
      if(k<external_count) {
        op.regional_boundary(i,internal_count+k)=X(i,k)*(-nu);
        op.regional_boundary(L+i,internal_count+k)=X(i,k)*nu*std::exp((nu-1)*s);
      }
    }
  }
  for(int p=1;p<=mehler_count+1;++p) {
    const int radial_p=p<=mehler_count?p:0;
    const int row=2*L+p-1,column=internal_count+external_count+p-1;
    double value=0,derivative=0;
    if(radial_p)Mehler::Mehler_1rank(m,radial_p*std::acos(-1.0)/s,
                                   geometry.aperture_deg,geometry.aperture_deg,value,derivative);
    else Constant::Constant_1rank(m,geometry.aperture_deg,geometry.aperture_deg,value,derivative);
    if(!std::isfinite(value) || value==0)throw std::runtime_error("Side angular normalization is unresolved");
    op.regional_boundary(row,column)=value;
    op.modes.push_back({radial_p?"mehler":"constant",radial_p,radial_p?-.5:0});
    for(int j=0;j<L;++j) {
      const int n=j+op.first_degree;
      // V_SH has prefactor R_E, V_R has prefactor a. The extra R_E/a
      // belongs here even for an internal-only SH input.
      op.boundary(row,j)=t*side.value[n]*radial_moment(geometry,n,radial_p,false);
      op.boundary(row,L+j)=t*side.value[n]*radial_moment(geometry,n,radial_p,true);
      op.forward(column,j)=op.boundary(row,j)/value;
      op.forward(column,L+j)=op.boundary(row,L+j)/value;
    }
  }
  finite(op.forward);finite(op.boundary);finite(op.regional_boundary);
  return op;
}
inline double duality_residual(const BoundaryOperators& op,bool spherical_only=false) {
  const auto product=multiply(op.regional_boundary,op.forward);
  std::vector<double> difference,reference;
  const int end=spherical_only?2*(op.degree-op.first_degree+1):op.boundary.rows;
  for(int i=0;i<end;++i)for(int j=0;j<op.boundary.cols;++j) {
    difference.push_back(product(i,j)-op.boundary(i,j));reference.push_back(op.boundary(i,j));
  }
  return norm(difference)/std::max(norm(reference),1e-300);
}
} // namespace SHRSCHA
