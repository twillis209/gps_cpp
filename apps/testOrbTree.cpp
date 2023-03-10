#include <PPEcdf.hpp>
#include <orbtree.h>
#include <map>
#include <vector>
#include <iostream>
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
  vector v({.1, .2, .3, .4, .4, .5});

  size_t n = u.size();

  // TODO need to handle the case wherein we have duplicate double-double pairs (surely quite rare, though?)
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
    cout << "i: " << i << endl;

    cout << "v_sorted[i]: " << v_sorted[i] << endl;

    double_multiset::iterator it = v_set.insert(v_sorted[i]).first;

    double_multiset::size_type m = v_set.rank(it);

    cout << "m: " << m << endl;

    size_t ix = posMap[{u_sorted[i], v_sorted[i]}];

    cout << "ix: " << ix << endl << endl;

    ecdf[ix] = (double) (m+1) / n;
  }

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
