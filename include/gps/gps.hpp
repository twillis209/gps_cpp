#pragma once

#include <vector>

double gps(std::vector<double> u,
                       std::vector<double> v);

std::vector<double> bivariateEcdfLW(const std::vector<double>& u,
                       const std::vector<double>& v);
