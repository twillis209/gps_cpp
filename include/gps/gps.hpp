#pragma once

#include <vector>
#include <map>
#include <boost/random.hpp>

namespace gps {
  std::vector<double> ecdf(std::vector<double> reference);

  double gpsStat(std::vector<double> u, std::vector<double> v);

  std::vector<double> bivariateEcdfLW(const std::vector<double>& u,
                       const std::vector<double>& v);

  std::vector<double> bivariateEcdfPar(const std::vector<double>& u, const std::vector<double>& v);

  std::vector<double> mix_rexp(size_t n, boost::mt19937& mt, double altRate, double altWeight, bool pvalScale);

  std::vector<double> rgps(size_t n, boost::mt19937& mt, double altRate, double altWeight, size_t noOfSnps);

  std::vector<double> permuteAndSampleGps(std::vector<double> u, std::vector<double> v, size_t n);

  std::vector<double> perturbDuplicates(std::vector<double> values);

  std::map<double, int> returnFreqMap(std::vector<double> values);
}
