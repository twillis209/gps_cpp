#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>
#include <gps.hpp>
#include <rapidcsv.h>

using namespace gps;
using namespace Eigen;

// TODO causes segfault when used in the for loop to initialise the u and v vectors
double uniformRand(boost::mt19937& generator) {
  static boost::uniform_01<boost::mt19937> dist(generator);
  return dist();
}

int main(int argc, const char* argv[]) {
  size_t n = 5;

  boost::mt19937 mt(43u);

  std::vector<double> gps_sample = rgps(n, 5, 0, 1e4, mt);

  for(auto i: gps_sample) std::cout << i << std::endl;

  return 0;
}
