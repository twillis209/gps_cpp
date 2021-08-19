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

//  boost::mt19937 generator;
//
//  std:: cout << uniformRand(generator) << ' ' << uniformRand(generator) << std::endl;

//  Document doc("/home/tw395/rds/hpc-work/gps/test/data/1e6_unif.csv");
//
//  std::vector<double> u = doc.GetColumn<double>("u");
//  std::vector<double> v = doc.GetColumn<double>("v");
//
//  std::cout << u[0] << ' ' << u[1] << std::endl;
//  std::cout << v[0] << ' ' << v[1] << std::endl;

  std::vector u({.1, .2, .3, .4, .5});
  std::vector v({.1, .2, .3, .5, .4});

  double gpsResult = gpsStat(u,v);

  std::cout << gpsResult << std::endl;

//  double gpsResult = gps(u,v);
//
//  std::cout << gpsResult << std::endl;

  return 0;
}
