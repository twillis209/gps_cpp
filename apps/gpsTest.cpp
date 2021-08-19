#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>
#include <gps.hpp>
#include <rapidcsv.h>

using namespace gps;
using namespace Eigen;
using namespace rapidcsv;

// TODO causes segfault when used in the for loop to initialise the u and v vectors
double uniformRand(boost::mt19937& generator) {
  static boost::uniform_01<boost::mt19937> dist(generator);
  return dist();
}

int main() {
  /*
  size_t n = 1e3;
  std::cout << "size_t: " << n << std::endl;
  std::vector<double> u(n);
  std::vector<double> v(n);
  boost::mt19937 generator;

  for(int i = 0; i < n; ++i) {
    u[i] = uniformRand(generator);
    v[i] = uniformRand(generator);
  }

  std::cout << u[0] << ' ' << u[n-1] << std::endl;
  std::cout << v[0] << ' ' << v[n-1] << std::endl;
  */

  Document doc("test/data/1e3_unif.csv");

  std::vector<double> u = doc.GetColumn<double>("u");
  std::vector<double> v = doc.GetColumn<double>("v");

  double gpsTest = gpsStat(u,v);

  std::cout << gpsTest << std::endl;

  return 0;
}
