#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <iostream>
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>
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

  Document data("pruned_pid_pso.csv");

  std::vector<double> u = data.GetColumn<double>("P.A");
  std::vector<double> v = data.GetColumn<double>("P.B");

  std::cout.precision(20);

  std::vector<double> uNoDup = perturbDuplicates(u);
  std::vector<double> vNoDup = perturbDuplicates(v);

  std::sort(uNoDup.begin(), uNoDup.end());

  std::cout << (std::adjacent_find(uNoDup.begin(), uNoDup.end()) == uNoDup.end()) << std::endl;

  std::cout << (std::adjacent_find(vNoDup.begin(), vNoDup.end()) == vNoDup.end()) << std::endl;

  std::vector<double> vNoDupCopy = vNoDup;

  std::sort(vNoDupCopy.begin(), vNoDupCopy.end());

  vNoDupCopy[vNoDupCopy.size()-1] = 0.0;

  std::cout << gpsStat(vNoDup, vNoDupCopy) << std::endl;

  /*
  for(int i = 0; i < u.size(); ++i) std::cout << u[i] << ',' << uNoDup[i] << std::endl;
  */

    return 0;
}
