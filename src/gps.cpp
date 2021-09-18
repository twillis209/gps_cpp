#include <gps.hpp>
#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <boost/math/distributions/empirical_cumulative_distribution_function.hpp>
#include <boost/random.hpp>
#include <math.h>
#include <map>

using namespace Eigen;
using boost::math::empirical_cumulative_distribution_function;

namespace gps {

  std::vector<double> ecdf(std::vector<double> reference) {
    std::vector<double> refCopy = reference;

    size_t n = reference.size();

    std::sort(refCopy.begin(), refCopy.end());

    std::vector<double> estEcdf;

    for(int i = 0; i < n; ++i) {
      estEcdf.push_back((double) (std::upper_bound(refCopy.begin(), refCopy.end(), reference[i])-refCopy.begin())/n);
    }

    return estEcdf;
  }

  double gpsStat(std::vector<double> u, std::vector<double> v, bool lw = false) {
  if(u.size() != v.size()) {
    throw std::invalid_argument("Size of u and v differs.");
  }

  size_t n = u.size();

if(std::distance(u.begin(), std::max_element(u.begin(), u.end())) == std::distance(
v.begin(), std::max_element(v.begin(), v.end()))) {
    throw std::invalid_argument("Indices of largest elements of u and v coincide. GPS statistic is undefined in this case.");
  }

  std::vector<double> uCopy = u;
  std::vector<double> vCopy = v;

  std::sort(uCopy.begin(), uCopy.end());
  std::sort(vCopy.begin(), vCopy.end());

  std::vector<double> uEcdf;
  std::vector<double> vEcdf;

  for(int i = 0; i < n; ++i) {
    uEcdf.push_back((double) (std::upper_bound(uCopy.begin(), uCopy.end(), u[i])-uCopy.begin())/n);
    vEcdf.push_back((double) (std::upper_bound(vCopy.begin(), vCopy.end(), v[i])-vCopy.begin())/n);
  }

  std::vector<double> bivariateEcdf;

  /*
    sort vector, count number of duplicates, then perturb instances of same value requisite number of times

    save duplicates in map, then go through unsorted vector and perturb with a multiple of eps based on remaining count as long as count > 1 (don't need to perturb all duplicates, just n - 1)

    TODO could just omit duplicate values?

    Need to perturb by a relevant value; add a fixed constant will have an effect which differs based on scale


  std::map<double, int> uFreqMap;
  std::map<double, int> vFreqMap;

  for(size_t i = 0; i < n; ++i){
    uFreqMap[u[i]]++;
    vFreqMap[v[i]]++;
  }

  for(size_t i = 0; i < n; ++i) {
    if(uFreqMap[u[i]] > 1) {
      u[i] = u[i] + (uFreqMap[u[i]] * std::numeric_limits<double>::epsilon());
      uFreqMap[u[i]]--;
    }

    if(vFreqMap[v[i]] > 1) {
      v[i] = v[i] + (vFreqMap[v[i]] * std::numeric_limits<double>::epsilon());
      vFreqMap[v[i]]--;
    }
  }

  */

  if(lw) {
    bivariateEcdf = bivariateEcdfLW(u,v);
  } else {
    bivariateEcdf = bivariateEcdfPar(u,v);
  }

  double max = -std::numeric_limits<double>::max();

  // Referencing v[i] and u[i] leads to a segfault as I have moved them

  for(int i = 0; i < n; ++i) {
    double cdf_uv = bivariateEcdf[i];
    double cdf_u = uEcdf[i];
    double cdf_v = vEcdf[i];
    double denom = sqrt(cdf_u*cdf_v - pow(cdf_u, 2.)*pow(cdf_v, 2.));

    double numerator = abs(cdf_uv - cdf_u*cdf_v);

    double maximand = sqrt((double) n / log((double) n))*numerator/denom;

    if(maximand > max) max = maximand;
  }

  return max;
}

std::vector<double> bivariateEcdfLW(const std::vector<double>& u, const std::vector<double>& v) {
  if(u.size() != v.size()) {
    throw std::invalid_argument("Size of u and v differs.");
  }

  size_t n = u.size();

  ArrayXd toAdd = ArrayXd::Constant(n, 1.);

  // NB: XXd for two-dimensional and double entries
  ArrayXXd ptr(2, n);

  for(int i = 0; i < n; i++) {
    ptr(0, i) = u[i];
    ptr(1, i) = v[i];
  }

  ArrayXd ecdf_arr = StOpt::fastCDFOnSample(ptr, toAdd);

  return std::vector<double>(ecdf_arr.data(), ecdf_arr.data() + ecdf_arr.size());
}

  std::vector<double> bivariateEcdfPar(const std::vector<double>& u, const std::vector<double>& v) {
    if(u.size() != v.size()) {
      throw std::invalid_argument("Size of u and v differs.");
    }

    size_t n = u.size();

    std::vector<double> ecdf;

    ecdf.reserve(n);

    #pragma openmp parallel for
    for(size_t i = 0; i < n; ++i) {
        int count = 0;

        for(size_t j = 0; j < n; ++j) {
          if(u[j] <= u[i] && v[j] <= v[i]) {
            ++count;
          }
        }

        ecdf[i] = (double) count/n;
      }

    return ecdf;
  }

  std::vector<double> mix_rexp(size_t n, boost::mt19937& mt, double altRate = 5, double altWeight = 0.01, bool pvalScale = false) {
  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<>> exp_1(mt, boost::exponential_distribution<>(1));

  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<>> exp_alt(mt, boost::exponential_distribution<>(altRate));

  boost::uniform_01<boost::mt19937&> unif(mt);

  std::vector<double> sample;

  sample.reserve(n);

  for(int i = 0; i < n; ++i) {
    double variate = (unif() <= altWeight ? exp_alt() : exp_1());
    sample.push_back(variate);
  }

  if(pvalScale) std::transform(sample.begin(), sample.end(), sample.begin(), [](double d) -> double { return exp(-d);});

  return sample;
}

  std::vector<double> rgps(size_t n, boost::mt19937& mt, double altRate = 5, double altWeight = 0.01, size_t noOfSnps = 1e4, bool lw = false) {

  std::vector<double> gps_sample;

  gps_sample.reserve(n);

  for(int i = 0; i < n; ++i) {
    std::vector<double> expSample = mix_rexp(2*noOfSnps, mt, altRate, altWeight, true);

    gps_sample.push_back(
                         gpsStat(
                                 std::vector<double>(expSample.begin(), expSample.begin() + noOfSnps),
                                 std::vector<double>(expSample.begin() + noOfSnps, expSample.end()), lw)
                         );
  }

  return gps_sample;
}

}
