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
  size_t n = 9e4;
  std::cout << "size_t: " << n << std::endl;
  std::vector<double> u(n);
  std::vector<double> v(n);
  /*
  boost::mt19937 generator;

  for(int i = 0; i < n; ++i) {
    u[i] = uniformRand(generator);
    v[i] = uniformRand(generator);
  }
  */

  /*
  Document doc("/home/tw395/rds/hpc-work/gps/out/1e6_unif.csv");

  std::vector<double> u = doc.GetColumn<double>("u");
  std::vector<double> v = doc.GetColumn<double>("v");
  Document doc("/home/tw395/rds/hpc-work/gps/out/1e6_unif.csv");


  doc.SetColumn<double>("gps", gpsTest);

  doc.Save("/home/tw395/rds/hpc-work/gps/out/1e6_unif.csv");
  */

  /*
  double gpsTest = gpsStat(u,v);

  std::cout << gpsTest << std::endl;
  */

  /*
  boost::mt19937 seed(5u);
  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<>> random_n(seed, boost::exponential_distribution<>(1)) ;

  for(int i = 0; i < 10; ++i) std::cout << random_n() << std::endl;
  */

  std::vector<double> rexp_sample = mix_rexp(100, 5, 0.01, false, 42u);

  for(auto i: rexp_sample) std::cout << i << std::endl;

  return 0;
}
