#include <gps.hpp>
#include <Eigen/Dense>
#include <fastCDFOnSample.h>
#include <boost/math/distributions/empirical_cumulative_distribution_function.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <math.h>
#include <map>
#include <omp.h>
#include <numeric>

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

  double gpsStat(std::vector<double> u, std::vector<double> v, std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)> bivariateEcdf) {
    if(u.size() != v.size()) {
      throw std::invalid_argument("Size of u and v differs.");
    }

    size_t n = u.size();

    if(std::distance(u.begin(), std::max_element(u.begin(), u.end())) == std::distance(
    v.begin(), std::max_element(v.begin(), v.end()))) {
        throw std::invalid_argument("Indices of largest elements of u and v coincide. GPS statistic is undefined in this case.");
      }

    std::vector<double> uEcdf = ecdf(u);
    std::vector<double> vEcdf = ecdf(v);

    std::vector<double> bivariateEcdfVec = bivariateEcdf(u, v);

    double max = -std::numeric_limits<double>::max();

    for(int i = 0; i < n; ++i) {
      double cdf_uv = bivariateEcdfVec[i];
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

    #pragma omp parallel for
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

  // TODO can make this faster by fitting the univariate ecdfs once
  std::vector<double> permuteAndSampleGps(std::vector<double> u, std::vector<double> v, size_t n, std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&)> bivariateEcdf) {
    std::vector<double> sample;

    size_t i = 0;

    while(i < n) {
      boost::range::random_shuffle(v);

      if(std::distance(u.begin(), std::max_element(u.begin(), u.end())) != std::distance(v.begin(), std::max_element(v.begin(), v.end()))) {
        sample.push_back(gpsStat(u, v, bivariateEcdf));
        ++i;
      }

    }

    return sample;
  }

  std::vector<double> perturbDuplicates(std::vector<double> values) {
    std::map<double, int> freqMap;

    for(size_t i = 0; i < values.size(); ++i){
      freqMap[values[i]]++;

      if(freqMap[values[i]] > 1) {
        for(int j = 1; j < freqMap[values[i]]; ++j) {
          values[i] = nextafter(values[i], 1.0);
        }
      }
    }

    return values;
  }


  std::vector<double> perturbDuplicates_addEpsilon(std::vector<double> values, double multiple) {
    std::map<double, int> freqMap;

    for(size_t i = 0; i < values.size(); ++i) {
      freqMap[values[i]]++;

      if(freqMap[values[i]] > 1) {
        double candidate_replacement = values[i] + (double) (freqMap[values[i]])* multiple*std::numeric_limits<double>::epsilon();
        //double candidate_replacement = std::nextafter(values[i], 1.0);

        while((0.5 * (candidate_replacement + values[i])) == values[i]) {
          std::cout << "Replacement is equal, incrementing" << std::endl;
          //candidate_replacement = std::nextafter(values[i], 1.0);
          candidate_replacement += (freqMap[values[i]])*multiple*std::numeric_limits<double>::epsilon();
        }

        values[i] = candidate_replacement;
      }
    }

    return values;
  }

  std::map<double, int> returnFreqMap(std::vector<double> values) {
    std::map<double, int> freqMap;

    for(size_t i = 0; i < values.size(); ++i){
      freqMap[values[i]]++;
    }

    return freqMap;
  }

  std::vector<double> simplePerturbation(std::vector<double> values) {
    https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
    std::vector<size_t> indices(values.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j){ return values[i] < values[j];} );

    //std::cout << "values.size(): " << values.size() << std::endl;

    std::cout.precision(20);

    for(size_t i = 0; i < (indices.size()-1); ++i) {
      //std::cout << "Incrementing " << i << " and " << i+1 << std::endl;
      //std::cout << "Incrementing " << indices[i] << " and " << indices[i+1] << std::endl;
      while((0.5 * (values[indices[i]] + values[indices[i+1]])) == values[indices[i+1]]) {
        //std::cout << "Incrementing " << values[indices[i]] << " and " << values[indices[i+1]] << std::endl;
        values[indices[i+1]] += 1000.0*std::numeric_limits<double>::epsilon();
      }
    }

    return values;
  }

}
