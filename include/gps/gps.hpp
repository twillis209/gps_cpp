#pragma once

#include <vector>

namespace gps {

double gpsStat(std::vector<double> u,
                       std::vector<double> v);

std::vector<double> bivariateEcdfLW(const std::vector<double>& u,
                       const std::vector<double>& v);

void bivariateEcdfLW_no_return(const std::vector<double>& u,
                                      const std::vector<double>& v);
}
