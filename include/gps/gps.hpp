#pragma once

#include <vector>
#include <map>
#include <boost/random.hpp>

namespace gps {
  std::vector<double> ecdf(std::vector<double> reference);

  double gpsStat(std::vector<double> u, std::vector<double> v, std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)> bivariateEcdf);

  std::vector<double> bivariateEcdfLW(const std::vector<double>& u,
                       const std::vector<double>& v);

  std::vector<double> bivariateEcdfPar(const std::vector<double>& u, const std::vector<double>& v);

  std::vector<double> permuteAndSampleGps(std::vector<double> u, std::vector<double> v, size_t n, std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)> bivariateEcdf);

  std::vector<double> perturbDuplicates(std::vector<double> values);

  std::vector<double> perturbDuplicates_addEpsilon(std::vector<double> values, double multiple);

  std::map<double, int> returnFreqMap(std::vector<double> values);

  std::vector<double> simplePerturbation(std::vector<double> values);
}
