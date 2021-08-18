#include <gps.hpp>
#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <boost/math/distributions/empirical_cumulative_distribution_function.hpp>
#include <math.h>

using namespace Eigen;
using boost::math::empirical_cumulative_distribution_function;

// TODO can probably do away with some of the copying ehre
double gps(std::vector<double> u, std::vector<double> v) {
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

  ArrayXXd ptr(2, n);

  for(int i = 0; i < n; i++) {
    ptr(0, i) = u[i];
    ptr(1, i) = v[i];
  }

  ArrayXd ecdf_arr = StOpt::fastCDFOnSample(ptr, toAdd);

  return std::vector<double>(ecdf_arr.data(), ecdf_arr.data() + ecdf_arr.size());
}
