#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>
#include <gps.hpp>
#include <rapidcsv.h>
#include <fastCDFOnSample.h>

using namespace gps;
using namespace Eigen;

int main(int argc, const char* argv[]) {
  boost::mt19937 mt(42u);

  size_t n = 1e6;

  /*
  std::vector<double> u;
  std::vector<double> v;

  for(int i = 0; i < n; ++i) {
    u.push_back(unif());
    v.push_back(unif());
  }

  double gps = gpsStat(u, v);

  std::cout << gps << std::endl;
  */

  /*
    boost::uniform_01<boost::mt19937&> unif(mt);
  std::vector<double> expSample = mix_rexp(2e4, mt, 5, 0, true);
  std::vector<double> gpsSample = permuteAndSampleGps(std::vector<double>(expSample.begin(), expSample.begin() + 1e4), std::vector<double>(expSample.begin() + 1e4, expSample.end()), 10);

  for(auto i: gpsSample) std::cout << i << std::endl;
  */

  /*
  size_t n = 1e2;

  std::vector<double> gps_sample = rgps(n, 5, 0, 1e4, 42u);

  //  for(auto i: gps_sample) std::cout << i << std::endl;


  std::vector<double> gps_sample = rgps(5, mt, 5, 0, 1e4);

  for(auto i: gps_sample) std::cout << i << std::endl;


  std::vector<int> vec = {0, 1, 2, 3, 4};

  std::vector<int> shuffledVec = shuffleVec(vec, mt);


  for(auto i: shuffledVec) std::cout << i << std::endl;

  std::cout << std::endl;

  for(auto i: vec) std::cout << i << std::endl;
  */

  return 0;
}
