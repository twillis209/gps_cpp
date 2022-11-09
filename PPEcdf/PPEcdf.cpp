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

  vector<double> bivariatePPEcdf(vector<double> u, vector<double> v) {
    assert(u.size()  == v.size());

    size_t n = u.size();

    map<pair<double, double>, size_t> posMap;

    // TODO duplicate keys? Perhaps could handle with vector
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

      cout << u_sorted[i] << " " << v_sorted[i] << " " << " " << m << " " << ecdf[i] << " " << endl;

      //ecdf[posMap[{u_sorted[i], v_sorted[i]}]] = (double) (m+1) / n;

    }

    cout << endl;

    //for(size_t i = 0; i < n; i++) {
    //  cout << u[i] << " " << v[i] << " " << ecdf[i] << " " << endl;
    //}

    return ecdf;
  }

 
  // From https://stackoverflow.com/a/12399290
  template<typename T>
  vector<size_t> idxSort(vector<T> v) {

    vector<pair<T, size_t>> v_pair(v.size());

    for(size_t i = 0; i < v.size(); i++) {
      v_pair[i] = make_pair(v[i], i);
    }

    stable_sort(v_pair.begin(), v_pair.end());

    vector<size_t> idx(v.size());

    for(size_t i = 0; i < v.size(); i++) {
      idx[i] = v_pair[i].second;
    }

    for(auto x: v_pair) cout << x.first << " " << x.second << endl;

    return idx; 
  }

  vector<double> reindex(const vector<double>& v, const vector<size_t>& idx) {
    assert(v.size() == idx.size());

    vector<double> reindexed(v.size());

    for(size_t i = 0; i < idx.size(); i++) {
      reindexed[idx[i]] = v[i];
    }

    return reindexed;
  }
}
