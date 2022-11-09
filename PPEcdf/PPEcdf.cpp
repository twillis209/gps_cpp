#include <PPEcdf.hpp>
#include <map>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ranked_index.hpp>
#include <boost/multi_index/member.hpp>

using namespace boost::multi_index;

namespace PPEcdf {

  std::vector<double> bivariatePPEcdf(std::vector<double> u, std::vector<double> v) {
    assert(u.size()  == v.size());

    size_t n = u.size();

    std::map<std::pair<double, double>, size_t> posMap;

    // TODO duplicate keys? Perhaps could handle with vector
    for(size_t i = 0; i < n; i++) {
      posMap[{u[i], v[i]}] = i;
    }

    std::vector<size_t> idx = idxSort(u);

    std::vector<double> u_sorted = reindex(u, idx);

    std::vector<double> v_sorted = reindex(v, idx);

    std::vector<double> ecdf(n);

    typedef multi_index_container<double, indexed_by<ranked_non_unique<identity<double>>>> double_multiset;

    double_multiset v_set;

    for(size_t i = 0; i < n; i++) {
      double_multiset::iterator it = v_set.insert(v_sorted[i]).first;

      double_multiset::size_type m = v_set.rank(it);

      ecdf[posMap[{u_sorted[i], v_sorted[i]}]] = (double) (m+1) / n;

    }

    return ecdf;
  }

  // From https://stackoverflow.com/a/12399290
  std::vector<size_t> idxSort(const std::vector<double>& v) {

    std::vector<size_t> idx(v.size());

    std::iota(idx.begin(), idx.end(), 0);

    // TODO Is this still stable?
    std::stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] <= v[i2];});

    return idx; 
  }

  std::vector<double> reindex(const std::vector<double>& v, const std::vector<size_t>& idx) {
    assert(v.size()  == idx.size());

    std::vector<double> u(v.size());

    for(size_t i = 0; i < idx.size(); i++) {
      u[idx[i]] = v[i];
    }

    return u;
  }
}
