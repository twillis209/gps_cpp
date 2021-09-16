#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>
#include <gps.hpp>
#include <omp.h>

using namespace gps;
using namespace Eigen;

int main(int argc, const char* argv[]) {
  omp_set_num_threads(2);

  std::cout << omp_get_num_threads() << std::endl;

  size_t n = 1e4;

  boost::mt19937 mt(42u);
  boost::uniform_01<boost::mt19937&> unif(mt);

  std::vector<double> u;
  std::vector<double> v;

  for(int i = 0; i < n; ++i) {
    u.push_back(unif());
    v.push_back(unif());
  }

  std::vector<double> ecdf = bivariateEcdfPar(u,v);

  /*
  ArrayXd toAdd = ArrayXd::Constant(n, 1.);

  ArrayXXd ptr(2, n);

  for(int i = 0; i < n; i++) {
    ptr(0, i) = unif();
    ptr(1, i) = unif();
  }

  ArrayXd ecdf_arr = StOpt::fastCDFOnSample(ptr, toAdd);
  */

  /*
  for(int i = 0; i < n; ++i) {
    u.push_back(unif());
    v.push_back(unif());
  }

  double gps = fastCDFOnSample(u,v);

  std::cout << gps << std::endl;


  std::vector<double> gps_sample = rgps(1, mt, 5, 0, 1e6);

  for(auto i: gps_sample) std::cout << i << std::endl;
  */

  return 0;
}
