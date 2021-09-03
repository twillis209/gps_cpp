#pragma once

#include <vector>

namespace gps {

double gpsStat(std::vector<double> u,
                       std::vector<double> v);

std::vector<double> bivariateEcdfLW(const std::vector<double>& u,
                       const std::vector<double>& v);

std::vector<double> mix_rexp(int n, double altRate = 5, double altWeight = 0.01, bool pvalScale);
