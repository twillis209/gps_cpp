#include <PPEcdf.hpp>
#include <map>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>

namespace PPEcdf {

  std::vector<double> bivariatePPEcdf(std::vector<double> u, std::vector<double> v) {
    assert(u.size()  == v.size());

    size_t n = u.size();

    std::vector<size_t> idx = idxSort(u);

    std::vector<double> u_sorted = reindex(u, idx);

    std::vector<double> v_sorted = reindex(v, idx);

    std::multiset<double> v_set;

    std::vector<double> ecdf(n);

    for(size_t i = 0; i < n; i++) {
      auto it = v_set.insert(v_sorted[i]);

      size_t m = 0;

      while(it != v_set.begin()) {

      }

    }
  }

  // From https://stackoverflow.com/a/12399290
  std::vector<size_t> idxSort(const std::vector<double>& v) {

    std::vector<size_t> idx(v.size());

    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx; 
  }

  std::vector<double> reindex(const std::vector<double>& v, const std::vector<size_t>& idx) {
    assert(v.size()  == idx.size());

    std::vector<double> u(v.size());

    for(size_t i = 0; i < idx.size(); i++) {
      u[i] = v[idx[i]];
    }

    return u;
  }
}
