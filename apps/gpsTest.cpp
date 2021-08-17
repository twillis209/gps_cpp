#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>
#include <gps.hpp>

using namespace Eigen;

// TODO causes segfault when used in the for loop to initialise the u and v vectors
double uniformRand() {
  boost::mt19937 generator(42);
  static boost::uniform_01<boost::mt19937> uniformDist(generator);
  return uniformDist();
}

int main() {
//  size_t n = 100000;
//  std::vector<double> u(n);
//  std::vector<double> v(n);
//
//  boost::mt19937 generator;
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

  std::vector u({.1,.2,.3,.4,.5});
  std::vector v({.1,.2,.3,.5,.4});

  int u_max = std::distance(u.begin(), std::max_element(u.begin(), u.end()));
  int v_max = std::distance(v.begin(), std::max_element(v.begin(), v.end()));

  std::cout << "u_max: " << u_max << " v_max: " << v_max << std::endl;

  return 0;
}
