#include <PPEcdf.hpp>
#include <orbtree.h>
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

  vector<double> bivariatePPEcdf(vector<double> u, vector<double> v) {
    if(u.size() != v.size()) {
      throw invalid_argument("Size of u and v differs.");
    }

    size_t n = u.size();

    map<pair<double, double>, size_t> posMap;

    for(size_t i = 0; i < n; i++) {
      posMap[{u[i], v[i]}] = i;
    }

    vector<size_t> idx = idxSort(u);

    vector<double> u_sorted = reindex(u, idx);

    vector<double> v_sorted = reindex(v, idx);

    vector<double> ecdf(n, 0.0);

    typedef multi_index_container<double, indexed_by<ranked_non_unique<identity<double>>>> double_multiset;

    double_multiset v_set;

    for(size_t i = 0; i < n; i++) {
      double_multiset::iterator it = v_set.insert(v_sorted[i]).first;

      double_multiset::size_type m = v_set.rank(it);

      size_t ix = posMap[{u_sorted[i], v_sorted[i]}];

      ecdf[ix] = (double) (m+1) / n;
    }

    return ecdf;
  }

  vector<double> bivariatePPEcdfOrbTree(vector<double> u, vector<double> v) {
    if(u.size() != v.size()) {
      throw invalid_argument("Size of u and v differs.");
    }

    size_t n = u.size();

    map<pair<double, double>, size_t> posMap;

    for(size_t i = 0; i < n; i++) {
      posMap[{u[i], v[i]}] = i;
    }

    vector<size_t> idx = idxSort(u);

    vector<double> u_sorted = reindex(u, idx);

    vector<double> v_sorted = reindex(v, idx);

    vector<double> ecdf(n, 0.0);

    orbtree::rankmultiset<double, uint32_t> multiset;

    for(size_t i = 0; i < n; i++) {
      // TODO really want to be able to insert and get back rank
      multiset.insert(v_sorted[i]);

      uint32_t m = multiset.get_sum(v_sorted[i]);

      size_t ix = posMap[{u_sorted[i], v_sorted[i]}];

      ecdf[ix] = (double) (m+1) / n;
    }

    return ecdf;
  }

  // TODO maybe I'm calling it incorrectly? But the other implementation didn't work either.
  // From https://stackoverflow.com/a/12399290
  template <typename T> vector<size_t> idxSort(const vector<T> &v) {

    vector<size_t> idx(v.size());

    iota(idx.begin(), idx.end(), 0);

    stable_sort(idx.begin(), idx.end(),
         [&v](size_t left, size_t right) {return v[left] < v[right];});

    return idx;
}

  template <typename T> vector<T> reindex(const vector<T>& v, const vector<size_t>& idx) {
    assert(v.size() == idx.size());

    vector<double> reindexed(v.size());

    for(size_t i = 0; i < idx.size(); i++) {
      reindexed[i] = v[idx[i]];
    }

    return reindexed;
  }
}
