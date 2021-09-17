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
  //omp_set_dynamic(0);
  //omp_set_num_threads(2);

  //std::cout << omp_get_num_threads() << std::endl;

  Document data("u_v_unif_n_1999556.csv");

  std::vector<double> u = data.GetColumn<double>("u");
  std::vector<double> v = data.GetColumn<double>("v");

  double gpsResult = gpsStat(u,v, true);

  std::cout << gpsResult << std::endl;
  /*

  std::vector<double> ecdf = bivariateEcdfPar(data.GetColumn<double>("u"), data.GetColumn<double>("v"));

  for(int i = 0; i < 1e3; ++i) std::cout << ecdf[i] << std::endl;
  */

  /*
  boost::mt19937 mt(42u);
  boost::uniform_01<boost::mt19937&> unif(mt);
  */

  /*
  for(int i = 0; i < n; ++i) {
    //    u.push_back(i*1e-5);
    //    v.push_back(i*1e-5);
    double u_s = unif();
    double v_s = unif();
    u.push_back(u_s);
    v.push_back(v_s);
    std::cout << u_s << ',' << v_s << std::endl;
  }
  */
  std::vector<double> ecdf = bivariateEcdfLW(u,v);

  //std::cout << ecdf[(n-2)] << std::endl;

  /*
  std::sort(u.begin(), u.end());
  std::sort(v.begin(), v.end());

  bool isUDuplicated = std::adjacent_find(u.begin(), u.end()) != u.end();
  bool isVDuplicated = std::adjacent_find(v.begin(), v.end()) != v.end();

  std::cout << isUDuplicated << ' ' << isVDuplicated << std::endl;
  */
  /*
  ArrayXd toAdd = ArrayXd::Constant(n, 1.);

  ArrayXXd ptr(2, n);

  for(int i = 0; i < n; i++) {
    ptr(0, i) = unif();
    ptr(1, i) = unif();
  }

  ArrayXd ecdf_arr = StOpt::fastCDFOnSample(ptr, toAdd);
  */

  /*
  for(int i = 0; i < n; ++i) {
    u.push_back(unif());
    v.push_back(unif());
  }

  double gps = fastCDFOnSample(u,v);

  std::cout << gps << std::endl;


  std::vector<double> gps_sample = rgps(1, mt, 5, 0, 1e6);

  for(auto i: gps_sample) std::cout << i << std::endl;
  */

  return 0;
}
