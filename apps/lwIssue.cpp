#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <boost/random.hpp>
#include <iostream>

using namespace Eigen;

int main(int argc, const char* argv[]) {
  boost::mt19937 mt(42u);
  boost::uniform_01<boost::mt19937&> unif(mt);

  size_t n = 1e5;

  ArrayXd toAdd = ArrayXd::Constant(n, 1.);

  ArrayXXd ptr(2, n);

  for(int i = 0; i < n; ++i) {
    ptr(0, i) = unif();
    ptr(1, i) = unif();
  }

  ArrayXd cdf = StOpt::fastCDFOnSample(ptr, toAdd);

  return 0;
}
