#include <PPEcdf.hpp>
#include <orbtree.h>
#include <map>
#include <vector>
#include <iostream>
#include <numeric>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ranked_index.hpp>
#include <boost/multi_index/member.hpp>

using namespace std;
using namespace PPEcdf;
using namespace orbtree;
using namespace boost::multi_index;

int main(int argc, const char* argv[]) {

  //vector u({.1, .2, .3, .5, .4, .5});
  //vector v({.5, .4, .3, .1, .2, .1});
  // ecdf should be 0.1666667 0.1666667 0.1666667 0.3333333 0.1666667 0.3333333
  vector u({.1, .2, .3, .4, .5});
  vector v({.1, .2, .3, .5, .4});
  // ecdf should be .167, .333, .5, .667, 1.0, 1.0

  size_t n = u.size();

  // TODO need to handle the case wherein we have duplicate double-double pairs (surely quite rare, though?)
  // TODO could index by <double, double, size_t> to take into account duplicates, with the size_t noting the no. of duplicates so far + 1?

  map<pair<double, double>, size_t> maxIxMap;

  //vector<vector<double>::size_type> u_idx = idxSort(u);
  vector<vector<double>::size_type> u_idx = idxSort(u);

  vector<double> u_sorted = reindex(u, u_idx);

  vector<double> v_sorted = reindex(v, u_idx);

  //vector<vector<double>::size_type> v_idx(n, 0);
  vector<double> v_idx(n, 0);

  iota(v_idx.begin(), v_idx.end(), 0);

  // TODO Need to check that this in fact reindexes v_idx
  v_idx = reindex(v_idx, u_idx);

  // TODO this cannot be the way this needs to be done
  vector<size_t> v_idx_size_t(v_idx.begin(), v_idx.end());

  /*
  cout << "u_idx" << endl;
  for(auto value : u_idx) cout << value << endl;

  cout << "v_idx_size_t" << endl;
  for(auto value : v_idx_size_t) cout << value << endl;
  */

  vector<double> ecdf(n, 0.0);

  typedef multi_index_container<double, indexed_by<ranked_non_unique<identity<double>>>> double_multiset;

  double_multiset v_set;

  for(size_t i = 0; i < n; i++) {
    cout << "i: " << i << endl;

    cout << "u_sorted[i]: " << u_sorted[i] << " v_sorted[i]: " << v_sorted[i] << endl;

    auto it = v_set.insert(v_sorted[i]).first;

    double_multiset::size_type m = v_set.rank(it);

    maxIxMap[make_pair(u_sorted[i], v_sorted[i])] = m+1;
  }

  cout << endl;

  for(size_t i = 0; i < n; i++) {
    cout << "u_sorted[i]: " << u_sorted[i] << " v_sorted[i]: " << v_sorted[i] << endl;
    cout << "u_idx[i]: " << u_idx[i] << " v_idx[i]: " << v_idx[i] << endl;

    cout << "v_sorted[i]: " << maxIxMap[make_pair(u_sorted[i], v_sorted[i])] << endl;

    ecdf[u_idx[i]] = (double) maxIxMap[make_pair(u_sorted[i], v_sorted[i])] / n;
    //ecdf[v_idx_size_t[i]] = (double) max(m+freqMap[v_sorted[i]], i+(freqMap[u[i]]-1)) / n;
  }

  cout << endl;

  for(auto v : ecdf) cout << v << endl;

  return 0;
}
