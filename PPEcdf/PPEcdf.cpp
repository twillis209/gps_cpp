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
using namespace std;

namespace PPEcdf {

  vector<double> bivariatePPEcdf(const vector<double>& u, const vector<double>& v) {
    if(u.size() != v.size()) {
      throw invalid_argument("Size of u and v differs.");
    }

    size_t n = u.size();

    map<pair<double, double>, size_t> maxIxMap;

    vector<vector<double>::size_type> u_idx = idxSort(u);

    vector<double> u_sorted = reindex(u, u_idx);

    vector<double> v_sorted = reindex(v, u_idx);

    vector<double> ecdf(n, 0.0);

    typedef multi_index_container<double, indexed_by<ranked_non_unique<identity<double>>>> double_multiset;

    double_multiset v_set;

    for(size_t i = 0; i < n; i++) {
      auto it = v_set.insert(v_sorted[i]).first;

      double_multiset::size_type m = v_set.rank(it);

      maxIxMap[make_pair(u_sorted[i], v_sorted[i])] = m+1;
    }


    for(size_t i = 0; i < n; i++) {
      ecdf[u_idx[i]] = (double) maxIxMap[make_pair(u_sorted[i], v_sorted[i])] / n;
    }

    return ecdf;
  }

  // From https://stackoverflow.com/a/12399290
  template <typename T> vector<size_t> idxSort(const vector<T> &v) {

    vector<size_t> idx(v.size());

    iota(idx.begin(), idx.end(), 0);

    stable_sort(idx.begin(), idx.end(),
         [&v](size_t left, size_t right) {return v[left] < v[right];});

    return idx;
}

  template <typename T> vector<T> reindex(const vector<T>& v, const vector<typename vector<T>::size_type>& idx) {
    assert(v.size() == idx.size());

    vector<T> reindexed(v.size());

    /*
    for(typename vector<T>::size_type i = 0; i < idx.size(); i++) {
      reindexed[i] = v[idx[i]];
    }
    */
    for(typename vector<T>::size_type i = 0; i < idx.size(); i++) {
      reindexed[i] = v[idx[i]];
    }

    return reindexed;
  }
}
