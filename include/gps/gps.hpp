#pragma once

#include <vector>
#include <map>
#include <functional>
#include <optional>
#include <cstddef>

using namespace std;

namespace gps {
  vector<double> ecdf(const vector<double>& reference);

  double gpsStat(const vector<double>& u, const vector<double>& v, function<vector<double>(const vector<double>&, const vector<double>&)> bivariateEcdf, function<double (const double&, const double&, const double&)> weightFunction);

  double meanStat(vector<double> u, vector<double> v, function<vector<double>(const vector<double>&, const vector<double>&)> bivariateEcdf, function<double (const double&, const double&, const double&)> weightFunction);

  vector<double> bivariateEcdfPar(const vector<double>& u, const vector<double>& v);

  vector<double> permuteAndSampleStat(vector<double> u,
                                          vector<double> v,
                                          size_t n,
                                      function<double (vector<double>, vector<double>, function<vector<double>(const vector<double>&, const vector<double>&)>, function<double (const double&, const double&, const double&)>)> statistic,
                                          function<vector<double>(const vector<double>&, const vector<double>&)> bivariateEcdf,
                                      function<double (const double&, const double&, const double&)> weightFunction, const optional<unsigned int>& seed = nullopt);

  map<double, int> returnFreqMap(vector<double> values);

  double gpsWeight(double cdf_u, double cdf_v, double cdf_uv);

  double pseudoADWeight(double cdf_u, double cdf_v, double cdf_uv);

  double squareNumerator(double cdf_u, double cdf_v, double cdf_uv);

  double normNumerator(double cdf_u, double cdf_v, double cdf_uv);
}
