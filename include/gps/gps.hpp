#pragma once

#include <vector>
#include <map>
#include <boost/random.hpp>

using namespace std;

namespace gps {
  vector<double> ecdf(vector<double> reference);

  double gpsStat(vector<double> u, vector<double> v, function<vector<double>(const vector<double>&, const vector<double>&)> bivariateEcdf, function<double (const double&, const double&, const double&)> weightFunction);

  double meanStat(vector<double> u, vector<double> v, function<vector<double>(const vector<double>&, const vector<double>&)> bivariateEcdf, function<double (const double&, const double&, const double&)> weightFunction);

  vector<double> bivariateEcdfLW(const vector<double>& u,
                       const vector<double>& v);

  vector<double> bivariateEcdfPar(const vector<double>& u, const vector<double>& v);

  vector<double> permuteAndSampleGps(vector<double> u,
                                          vector<double> v,
                                          size_t n,
                                          function<vector<double>(const vector<double>&, const vector<double>&)> bivariateEcdf,
                                          function<double (const double&, const double&, const double&)> weightFunction);

  vector<double> permuteAndSampleMeanStat(
                                               vector<double> u,
                                               vector<double> v,
                                               size_t n,
                                               function<double (vector<double>, vector<double>, function<vector<double>(const vector<double>&, const vector<double>&)>, function<double (const double&, const double&, const double&)>)> stat,
                                               function<vector<double>(const vector<double>&, const vector<double>&)> bivariateEcdf,
                                               function<double (const double&, const double&, const double&)> weightFunction);

  vector<double> perturbDuplicates(vector<double> values);

  vector<double> perturbDuplicates_addEpsilon(vector<double> values, double multiple);

  map<double, int> returnFreqMap(vector<double> values);

  vector<double> simplePerturbation(vector<double> values);

  map<double, int> returnFreqMap(vector<double> values);

  double gpsWeight(double cdf_u, double cdf_v, double cdf_uv);

  double pseudoADWeight(double cdf_u, double cdf_v, double cdf_uv);
}
