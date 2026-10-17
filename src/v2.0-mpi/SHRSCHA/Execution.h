#pragma once
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef SHRSCHA_USE_MPI
#include "../MPIParallel.h"
#endif

namespace SHRSCHA {
inline int execution_rank() {
#ifdef SHRSCHA_USE_MPI
  return MPIParallel::rank();
#else
  return 0;
#endif
}
inline int execution_size() {
#ifdef SHRSCHA_USE_MPI
  return MPIParallel::size();
#else
  return 1;
#endif
}
inline void combine_coefficients(std::vector<double>& values) {
#ifdef SHRSCHA_USE_MPI
  MPIParallel::allreduce_sum(values);
#else
  (void)values;
#endif
}
inline std::string combine_report(const std::string& local) {
#ifdef SHRSCHA_USE_MPI
  if(local.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))throw std::runtime_error("MPI report is too large");
  const int length=static_cast<int>(local.size());
  std::vector<int> counts(execution_rank()==0?execution_size():0),offsets(counts.size());
  if(MPI_Gather(&length,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD)!=MPI_SUCCESS)
    throw std::runtime_error("MPI report length gather failed");
  int total=0;
  for(std::size_t i=0;i<counts.size();++i) {
    if(counts[i]>std::numeric_limits<int>::max()-total)throw std::runtime_error("MPI gathered report is too large");
    offsets[i]=total;total+=counts[i];
  }
  std::string result(static_cast<std::size_t>(total),'\0');
  if(MPI_Gatherv(local.data(),length,MPI_CHAR,(result.empty()?nullptr:&result[0]),counts.data(),offsets.data(),MPI_CHAR,0,MPI_COMM_WORLD)!=MPI_SUCCESS)
    throw std::runtime_error("MPI report gather failed");
  return result;
#else
  return local;
#endif
}
template<typename Function> int guarded_application(Function function) {
#ifdef SHRSCHA_USE_MPI
  std::unique_ptr<MPIParallel::Environment> environment;
#endif
  try {
#ifdef SHRSCHA_USE_MPI
    int count=0;char** arguments=nullptr;
    environment.reset(new MPIParallel::Environment(count,arguments));
    if(execution_rank()!=0)std::cout.setstate(std::ios_base::failbit);
#endif
    return function();
  }catch(const std::exception& error) {
    std::cerr<<"Calculation failed: "<<error.what()<<'\n';
#ifdef SHRSCHA_USE_MPI
    if(MPIParallel::active())MPI_Abort(MPI_COMM_WORLD,1);
#endif
    return 1;
  }
}
} // namespace SHRSCHA
