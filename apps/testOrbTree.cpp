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
  vector u({.1, .2, .3, .4, .5, .5});
  vector v({.1, .2, .3, .4, .5, .5});

  size_t n = u.size();

  // TODO need to handle the case wherein we have duplicate double-double pairs (surely quite rare, though?)
  // TODO could index by <double, double, size_t> to take into account duplicates, with the size_t noting the no. of duplicates so far + 1?

  map<pair<double, double>, size_t> posMap;

  for(size_t i = 0; i < n; i++) {
    posMap[{u[i], v[i]}] = i;
  }

  vector<size_t> u_idx = idxSort(u);

  vector<double> u_sorted = reindex(u, u_idx);

  vector<double> v_sorted = reindex(v, u_idx);

  vector<size_t> v_idx(n);

  iota(v_idx.begin(), v_idx.end(), 0);

  // TODO expecting this to give index back into original v vector
  //v_idx = reindex(v_idx, u_idx);

  vector<double> ecdf(n, 0.0);

  map<double, size_t> freqMap;

  for(size_t i = 0; i < v_sorted.size(); ++i) {
    freqMap[v_sorted[i]]++;
  }

  typedef multi_index_container<double, indexed_by<ranked_non_unique<identity<double>>>> double_multiset;

  double_multiset v_set;

  for(size_t i = 0; i < n; i++) {
    cout << "i: " << i << endl;

    cout << "v_sorted[i]: " << v_sorted[i] << endl;

    v_set.insert(v_sorted[i]);
  }

  for(size_t i = 0; i < n; i++) {
    cout << "i: " << i << endl;

    cout << "v_sorted[i]: " << v_sorted[i] << endl;

    cout << "find(v_sorted[i]): " << endl;

    auto it = v_set.find(v_sorted[i]);

    double_multiset::size_type m = v_set.rank(it);

    cout << "rank: " << m << endl;

    size_t ix = posMap[{u_sorted[i], v_sorted[i]}];

    cout << "ix: " << ix << endl << endl;

    // TODO: can probably get the freqMap info from the multiset
    ecdf[ix] = (double) (m+freqMap[v_sorted[i]]) / n;

    /*
    cout << "m: " << m << endl;

    size_t ix = posMap[{u_sorted[i], v_sorted[i]}];

    cout << "ix: " << ix << endl << endl;

    ecdf[ix] = (double) (m+1) / n;
    */
  }

  auto it = v_set.lower_bound(.5);
  auto end = v_set.upper_bound(.5);

  //for(; it != end; ++it)
  //  cout << it << endl;
  /*
  for(size_t i = 0; i < n; i++) {
    cout << "i: " << i << endl;

    cout << "v_sorted[i]: " << v_sorted[i] << endl;

    double_multiset::iterator it = v_set.insert(v_sorted[i]).first;

    double_multiset::size_type m = v_set.rank(it);

    cout << "m: " << m << endl;

    size_t ix = posMap[{u_sorted[i], v_sorted[i]}];

    cout << "ix: " << ix << endl << endl;

    ecdf[ix] = (double) (m+1) / n;
  }
  */


  /*
  rankmultiset<double, uint32_t> multiset;

  for(size_t i = 0; i < n; i++) {
    // TODO really want to be able to insert and get back rank
    cout << "i: " << i << endl;
    cout << "v_sorted[i]: " << v_sorted[i] << endl;

    multiset.insert(v_sorted[i]);

    // counting those strictly less than
    uint32_t m = multiset.get_sum(v_sorted[i]);

    cout << "m: " << m << endl;

    size_t ix = posMap[{u_sorted[i], v_sorted[i]}];

    cout << "ix: " << ix << endl << endl;

    ecdf[ix] = (double) (m+1) / n;
  }
  */

  for(auto x: ecdf) cout << x << endl;

  rankmultimap<double, double, uint32_t> multimap;

  /*
  for(size_t i = 0; i < n; i++) {
    // TODO really want to be able to insert and get back rank
    multiset.insert(v_sorted[i]);

    uint32_t m = multiset.get_sum(v_sorted[i]);

    size_t ix = posMap[{u_sorted[i], v_sorted[i]}];

    ecdf[ix] = (double) (m+1) / n;
  }
  */

  return 0;
}
