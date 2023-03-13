#pragma once

#include <vector>
#include <cstddef>

namespace PPEcdf {
  std::vector<double> bivariatePPEcdf(const std::vector<double>& u, const std::vector<double>& v);
  template <typename T> std::vector<std::size_t> idxSort(const std::vector<T> &v);
  template <typename T> std::vector<T> reindex(const std::vector<T>& v, const std::vector<typename std::vector<T>::size_type>& idx);
}
