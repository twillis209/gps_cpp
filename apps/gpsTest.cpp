#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>
#include <gps.hpp>

using namespace Eigen;

// TODO causes segfault when used in the for loop to initialise the u and v vectors
double uniformRand(boost::mt19937& generator) {
  static boost::uniform_01<boost::mt19937> dist(generator);
  return dist();
}

int main() {
//  size_t n = 100000;
//  std::vector<double> u(n);
//  std::vector<double> v(n);
//
//  boost::normal_distribution<double> normalDistrib;
//  boost::variate_generator<boost::mt19937 &, boost::normal_distribution<double> > normalRand(generator, normalDistrib);
//
//  for(int i = 0; i < n; ++i) {
//    u[i] = normalRand();
//    v[i] = normalRand();
//  }
//
//  double gpsTest = gps(u,v);
//
//  std::cout << gpsTest << std::endl;
//

  boost::mt19937 generator;

  std:: cout << uniformRand(generator) << ' ' << uniformRand(generator) << std::endl;

  return 0;
}
