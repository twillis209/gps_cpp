#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>
#include <gps.hpp>
#include <omp.h>
#include <rapidcsv.h>

using namespace gps;
using namespace rapidcsv;
using namespace Eigen;

int main(int argc, const char* argv[]) {
  /*
  Document data("u_v_unif_n_1999556.csv");

  std::vector<double> u = data.GetColumn<double>("u");
  std::vector<double> v = data.GetColumn<double>("v");

  double gpsResult = gpsStat(u,v);

  std::cout << gpsResult << std::endl;

  std::vector<double> permSample = permuteAndSampleGps(u, v, 3);

  for(auto i : permSample) std::cout << i << std::endl;

  std::vector<double> vNoDup = perturbDuplicates(v);

  std::vector<double> u({.01042, .01042, .01059, .01059});

  std::vector<double> uNoDup = perturbDuplicates(u);

  std::cout.precision(10);

  std::sort(uNoDup.begin(), uNoDup.end());
  std::sort(vNoDup.begin(), vNoDup.end());
  */

  Document data("pruned_pid_sys-scl.csv");

  std::vector<double> u = data.GetColumn<double>("P.A");
  std::vector<double> v = data.GetColumn<double>("P.B");

  //std::cout.precision(20);

  std::vector<double> uNoDup = perturbDuplicates(u);
  std::vector<double> vNoDup = perturbDuplicates(v);

  double gpsResult = gpsStat(uNoDup, vNoDup);

  std::vector<std::vector<double>> gpsPermutations;

  #pragma omp parallel for
  for(int i = 0; i < 100; ++i) {
    gpsPermutations.push_back(permuteAndSampleGps(uNoDup, vNoDup, 100));
  }

  std::vector<double> gpsResults;

  for(int i = 0; i < 4; ++i){
    for(int j = 0; j < 10; ++j) {
      gpsResults.push_back(gpsPermutations[i][j]);
    }
  }

  for(auto i : gpsResults) std::cout << i << std::endl;

  //  for(auto i : gpsPermutations) std::cout << i << std::endl;

  //  for(int i = 0; i < v.size(); ++i) std::cout << v[i] << ',' << vNoDup[i] << std::endl;

  /*
  std::map<double, int> freqMap;

  for(size_t i = 0; i < u.size(); ++i){
    freqMap[u[i]]++;

    if(freqMap[u[i]] > 1) {
      //values[i] = values[i] + (freqMap[values[i]] * std::numeric_limits<double>::epsilon());
      u[i] = u[i] + (freqMap[u[i]] * 1e-4);
    }
  }
  .6507
  .4134

  boost::mt19937 mt(42u);

  boost::uniform_01<boost::mt19937&> unif(mt);

  size_t n = 1e6;

  std::vector<double> u;

  for(int i = 0; i < n; ++i) {
    u.push_back(unif());
  }

  */


  /*
  std::map<double, int> freqMap;

  for(size_t i = 0; i < n; ++i){
    freqMap[noDup[i]]++;
  }

  for(size_t i = 0; i < n; ++i){
    if(freqMap[noDup[i]] > 1) std::cout << u[i] << ',' << noDup[i] << std::endl;
  }
  */

  return 0;
}
