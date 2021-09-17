#include <gps.hpp>
#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <boost/math/distributions/empirical_cumulative_distribution_function.hpp>
#include <boost/random.hpp>
#include <math.h>

using namespace Eigen;
using boost::math::empirical_cumulative_distribution_function;

namespace gps {

// TODO can probably do away with some of the copying here
  double gpsStat(std::vector<double> u, std::vector<double> v, bool lw = false) {
  if(u.size() != v.size()) {
    throw std::invalid_argument("Size of u and v differs.");
  }

  size_t n = u.size();

  if(std::distance(u.begin(), std::max_element(u.begin(), u.end())) == std::distance(v.begin(), std::max_element(v.begin(), v.end()))) {
    throw std::invalid_argument("Indices of largest elements of u and v coincide. GPS statistic is undefined in this case.");
  }

  std::vector<double> u_copy = u;
  std::vector<double> v_copy = v;
  std::vector<double> bivariate_ecdf;
  if(lw) {
    bivariate_ecdf = bivariateEcdfLW(u,v);
  } else {
    bivariate_ecdf = bivariateEcdfPar(u,v);
  }

  // TODO hopefully this doesn't modify u or v by using std::move
  // TODO it does have to sort the data; is that done in-place?
  auto ecdf_u = empirical_cumulative_distribution_function(std::move(u));
  auto ecdf_v = empirical_cumulative_distribution_function(std::move(v));

  double max = -std::numeric_limits<double>::max();

  // Referencing v[i] and u[i] leads to a segfault as I have moved them

  for(int i = 0; i < n; ++i) {
    double cdf_uv = bivariate_ecdf[i];
    double cdf_u = ecdf_u(u_copy[i]);
    double cdf_v = ecdf_v(v_copy[i]);

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

        ecdf[i] = (double) count;
      }

    std::transform(ecdf.begin(), ecdf.end(), ecdf.begin(), [n](double d) -> double { return d/(double) n;});

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
