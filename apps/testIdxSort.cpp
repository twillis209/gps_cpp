#include <PPEcdf.hpp>
#include <iostream>

using namespace PPEcdf;
using namespace std;

int main(int argc, const char* argv[]) {
  vector<double> v({0.4, 0.1, 0.3, 0.2});
  vector<size_t> idx({3, 0, 2, 1});

  vector<size_t> sortedIdx = idxSort(v);

  for(auto x : sortedIdx) cout << x << endl;

  return 0;
}
