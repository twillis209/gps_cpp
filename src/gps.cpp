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
double gpsStat(std::vector<double> u, std::vector<double> v) {
  if(u.size() != v.size()) {
    throw std::invalid_argument("Size of u and v differs.");
  }

  size_t n = u.size();

  if(std::distance(u.begin(), std::max_element(u.begin(), u.end())) == std::distance(v.begin(), std::max_element(v.begin(), v.end()))) {
    throw std::invalid_argument("Indices of largest elements of u and v coincide. GPS statistic is undefined in this case.");
  }

  std::vector<double> u_copy = u;
  std::vector<double> v_copy = v;

  std::vector<double> bivariate_ecdf = bivariateEcdfLW(u,v);
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

  std::vector<double> mix_rexp(size_t n, double altRate = 5, double altWeight = 0.01, bool pvalScale = false, unsigned int seed = 42u) {
  boost::mt19937 mt(seed);

  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<>> exp_1(mt, boost::exponential_distribution<>(1));

  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<>> exp_alt(mt, boost::exponential_distribution<>(altRate));

  boost::uniform_01<boost::mt19937> unif(mt);

  std::vector<double> sample;

  sample.reserve(n);

  for(int i = 0; i < n; ++i) {
    double variate = (unif() <= altWeight ? exp_alt() : exp_1());
    sample.push_back(variate);
  }

  if(pvalScale) std::transform(sample.begin(), sample.end(), sample.begin(), [](double d) -> double { return exp(-d);});

  return sample;
}

  std::vector<double> rgps(size_t n, double altRate = 5, double altWeight = 0.01, size_t noOfSnps = 1e4, unsigned int seed = 42u) {

  boost::mt19937 mt(seed);

  boost::exponential_distribution<> expDist(1);
  boost::exponential_distribution<> expAlt(altRate);

  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<>> exp_1(mt, expDist);
  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<>> exp_alt(mt, expAlt);

  //boost::uniform_01<boost::mt19937> unif(mt);
  boost::uniform_01<boost::mt19937&> unif(mt);

  std::vector<double> gps_sample;

  gps_sample.reserve(n);

  for(int i = 0; i < n; ++i) {
      double gps = NULL;

      int j = 1;

      while(j < 6 && gps == NULL) {

        std::vector<double> exp_sample;

        exp_sample.reserve(2*noOfSnps);

        for(int k = 0; k < 2*noOfSnps; ++k) {
          double variate = (unif() <= altWeight ? exp_alt() : exp_1());
          exp_sample.push_back(variate);
        }

        std::transform(exp_sample.begin(), exp_sample.end(), exp_sample.begin(), [](double d) -> double { return exp(-d);});

        if(std::distance(exp_sample.begin(), std::max_element(exp_sample.begin(), exp_sample.begin() + noOfSnps)) != std::distance(exp_sample.begin() + noOfSnps, std::max_element(exp_sample.begin() + noOfSnps, exp_sample.end()))) {
          gps = gpsStat(std::vector<double>(exp_sample.begin(), exp_sample.begin() + noOfSnps), std::vector<double>(exp_sample.begin() + noOfSnps, exp_sample.end()));
        }

        ++j;
    }

      if(gps == NULL) {
        throw std::runtime_error("Failed to generate GPS sample realisation after 5 attempts");
      } else {
        gps_sample.push_back(gps);
      }

  }

  return gps_sample;
}

}
