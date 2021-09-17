#pragma once

#include <vector>
#include <boost/random.hpp>

namespace gps {

  double gpsStat(std::vector<double> u,
                 std::vector<double> v, bool lw);

  std::vector<double> bivariateEcdfLW(const std::vector<double>& u,
                       const std::vector<double>& v);

  std::vector<double> bivariateEcdfPar(const std::vector<double>& u, const std::vector<double>& v);

  std::vector<double> mix_rexp(size_t n, boost::mt19937& mt, double altRate, double altWeight, bool pvalScale);

  std::vector<double> rgps(size_t n, boost::mt19937& mt, double altRate, double altWeight, size_t noOfSnps, bool lw);

}
