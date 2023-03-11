#pragma once

#include <vector>

namespace PPEcdf {
  std::vector<double> bivariatePPEcdf(std::vector<double> u, std::vector<double> v);
  std::vector<double> bivariatePPEcdfOrbTree(std::vector<double> u, std::vector<double> v);
  template <typename T> std::vector<size_t> idxSort(const std::vector<T> &v);
  template <typename T> std::vector<T> reindex(const std::vector<T>& v, const std::vector<typename std::vector<T>::size_type>& idx);
}
