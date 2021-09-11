#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>
#include <gps.hpp>
#include <rapidcsv.h>

using namespace gps;
using namespace Eigen;

int main(int argc, const char* argv[]) {
  size_t n = 1e2;

  std::vector<double> gps_sample = rgps(n, 5, 0, 1e4, 42u);

  //  for(auto i: gps_sample) std::cout << i << std::endl;

  return 0;
}
