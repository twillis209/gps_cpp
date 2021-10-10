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
  Document data("pruned_pid_sys-scl.csv");

  std::vector<double> u = data.GetColumn<double>("P.A");
  std::vector<double> v = data.GetColumn<double>("P.B");

  std::vector<double> uNoDup = perturbDuplicates(u);
  std::vector<double> vNoDup = perturbDuplicates(v);

  std::vector<double> uEcdf = ecdf(uNoDup);
  std::vector<double> vEcdf = ecdf(vNoDup);

  std::vector<double> bivariateEcdf = bivariateEcdfLW(uNoDup, vNoDup);

  std::cout << "u,v,cdf_u,cdf_v,cdf_uv,denom,maximand" << std::endl;

  size_t n = uNoDup.size();

  for(int i = 0; i < n; ++i) {
    double cdf_uv = bivariateEcdf[i];
    double cdf_u = uEcdf[i];
    double cdf_v = vEcdf[i];
    double denom = sqrt(cdf_u*cdf_v - pow(cdf_u, 2.)*pow(cdf_v, 2.));

    double numerator = abs(cdf_uv - cdf_u*cdf_v);

    double maximand = sqrt((double) n / log((double) n))*numerator/denom;

    std::cout << uNoDup[i] << "," << vNoDup[i] << "," << cdf_u << "," << cdf_v << "," << cdf_uv << "," << denom << "," << maximand << std::endl;

  }

    return 0;
}
