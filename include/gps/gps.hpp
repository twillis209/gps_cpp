#pragma once

#include <vector>
#include <boost/random.hpp>

namespace gps {

  double rgps(size_t n);

  double gpsStat(std::vector<double> u,
                       std::vector<double> v);

  std::vector<double> bivariateEcdfLW(const std::vector<double>& u,
                       const std::vector<double>& v);

  std::vector<double> mix_rexp(size_t n, double altRate, double altWeight, bool pvalScale, boost::mt19937 mt);

  std::vector<double> rgps(size_t n, double altRate, double altWeight, size_t noOfSnps, unsigned int seed);

}
