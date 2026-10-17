#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace SHRSCHA {
struct Matrix {
  int rows=0, cols=0;
  std::vector<double> data;
  Matrix() = default;
  Matrix(int r,int c):rows(r),cols(c) {
    if(r<0 || c<0) throw std::invalid_argument("Negative matrix dimensions");
    data.assign(static_cast<std::size_t>(r)*c,0);
  }
  double& operator()(int r,int c) {return data.at(static_cast<std::size_t>(r)*cols+c);}
  double operator()(int r,int c) const {return data.at(static_cast<std::size_t>(r)*cols+c);}
};
inline Matrix transpose(const Matrix& a) {
  Matrix b(a.cols,a.rows);
  for(int i=0;i<a.rows;++i) for(int j=0;j<a.cols;++j) b(j,i)=a(i,j);
  return b;
}
inline Matrix multiply(const Matrix& a,const Matrix& b) {
  if(a.cols!=b.rows) throw std::invalid_argument("Incompatible matrix dimensions");
  Matrix c(a.rows,b.cols);
  for(int i=0;i<a.rows;++i) for(int k=0;k<a.cols;++k)
    for(int j=0;j<b.cols;++j) c(i,j)+=a(i,k)*b(k,j);
  return c;
}
inline std::vector<double> multiply(const Matrix& a,const std::vector<double>& x) {
  if(a.cols!=static_cast<int>(x.size())) throw std::invalid_argument("Invalid coefficient dimension");
  std::vector<double> b(a.rows);
  for(int i=0;i<a.rows;++i) {
    long double sum=0;
    for(int j=0;j<a.cols;++j) sum+=static_cast<long double>(a(i,j))*x[j];
    b[i]=static_cast<double>(sum);
  }
  return b;
}
inline double norm(const std::vector<double>& a) {
  long double s=0; for(double v:a) s+=static_cast<long double>(v)*v;
  return static_cast<double>(std::sqrt(s));
}
inline double relative_difference(const std::vector<double>& a,const std::vector<double>& b) {
  if(a.size()!=b.size()) throw std::invalid_argument("Invalid comparison dimensions");
  auto d=a; for(std::size_t i=0;i<a.size();++i)d[i]-=b[i];
  return norm(d)/std::max(norm(a),1e-300);
}
inline void finite(const Matrix& a) {
  for(double v:a.data) if(!std::isfinite(v))
    throw std::runtime_error("Non-finite matrix entry; geometry or degree is numerically unresolved");
}
struct Quadrature {std::vector<double> x,w;};
inline Quadrature gauss_interval(double lo,double hi,int n) {
  if(!std::isfinite(lo) || !std::isfinite(hi) || !(hi>lo) || n<2)
    throw std::invalid_argument("Invalid quadrature interval or count");
  Quadrature q; q.x.resize(n); q.w.resize(n);
  for(int i=0;i<(n+1)/2;++i) {
    long double z=std::cos(std::acos(-1.0L)*(i+.75L)/(n+.5L));
    long double derivative=0;
    bool converged=false;
    for(int step=0;step<100;++step) {
      long double p0=1,p1=z;
      for(int k=2;k<=n;++k) {
        const long double p=((2*k-1)*z*p1-(k-1)*p0)/k;
        p0=p1;p1=p;
      }
      derivative=n*(z*p1-p0)/(z*z-1);
      const long double change=p1/derivative;
      z-=change;
      if(std::fabs(change)<2e-16L) {converged=true;break;}
    }
    if(!converged) throw std::runtime_error("Gauss-Legendre root failed to converge");
    // Re-evaluate the derivative at the final root.
    long double p0=1,p1=z;
    for(int k=2;k<=n;++k) {const long double p=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=p;}
    derivative=n*(z*p1-p0)/(z*z-1);
    const long double weight=(hi-lo)/((1-z*z)*derivative*derivative);
    q.x[i]=static_cast<double>((hi+lo)/2-(hi-lo)*z/2);
    q.x[n-1-i]=static_cast<double>((hi+lo)/2+(hi-lo)*z/2);
    q.w[i]=q.w[n-1-i]=static_cast<double>(weight);
  }
  return q;
}
struct Angular {std::vector<double> value,derivative;};
// Independent real Schmidt semi-normalized integer-degree recurrence.
inline Angular sh_angular(int m,int N,double theta) {
  if(m<0 || N<m || !(theta>0 && theta<std::acos(-1.0)))
    throw std::invalid_argument("Invalid SH angular mode");
  Angular result;result.value.resize(N+1);result.derivative.resize(N+1);
  const double x=std::cos(theta),st=std::sin(theta);
  double diagonal=m?st:1;
  for(int j=2;j<=m;++j)diagonal*=std::sqrt((2.0*j-1)/(2.0*j))*st;
  double previous=diagonal,older=0;
  for(int n=m;n<=N;++n) {
    double p=diagonal;
    if(n==m+1)p=std::sqrt(2.0*m+1)*x*previous;
    else if(n>m+1) {
      const double denominator=n*double(n)-m*double(m);
      p=(2.0*n-1)/std::sqrt(denominator)*x*previous-
        std::sqrt(((n-1)*double(n-1)-m*double(m))/denominator)*older;
    }
    result.value[n]=p;
    result.derivative[n]=(n*x*p-(n>m?std::sqrt(n*double(n)-m*double(m))*previous:0))/st;
    older=previous;previous=p;
  }
  return result;
}
} // namespace SHRSCHA
