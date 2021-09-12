#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>

using namespace Eigen;

int main(int argc, const char* argv[]) {
  size_t n = 4e4;

  boost::mt19937 mt(42u);
  boost::uniform_01<boost::mt19937&> unif(mt);

  ArrayXd toAdd = ArrayXd::Constant(n, 1.);

  ArrayXXd ptr(2, n);

  for(int i = 0; i < n; i++) {
    ptr(0, i) = unif();
    ptr(1, i) = unif();
  }

  ArrayXd ecdf_arr = StOpt::fastCDFOnSample(ptr, toAdd);

  return 0;
}
