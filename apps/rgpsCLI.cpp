#include <iostream>
#include <gps.hpp>
#include <boost/program_options.hpp>
#include <boost/random.hpp>

namespace po = boost::program_options;
using namespace po;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  int cores;
  int draws;
  size_t noOfSnps;
  double altWeight;
  double altRate;
  unsigned int seed;

  desc.add_options()
    ("help", "Print help message")
    ("cores,c", po::value<int>(&cores), "No. of cores")
    ("draws,n", po::value<int>(&draws), "No. of draws per core")
    ("noOfSnps,p", po::value<size_t>(&noOfSnps)->default_value(1e4), "No. of SNPs")
    ("altWeight,w", po::value<double>(&altWeight)->default_value(0), "Mixing weight for non-standard exp. component")
    ("altRate,r", po::value<double>(&altRate)->default_value(5), "Rate parameter for non-standard exp. component")
    ("seed,s", po::value<unsigned int>(&seed)->default_value(42u), "Seed for PRNG")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("cores")) {
    boost::mt19937 mt(seed);

    boost::uniform_int<> dist(0, 255);

    boost::variate_generator<boost::mt19937&, boost::uniform_int<>> seedGen(mt, dist);

    unsigned int seeds[cores];

    for(size_t i = 0; i < cores; ++i) {
      seeds[i] = seedGen();
    }


    std::vector<std::vector<double>> samplesPerCore;

    #pragma omp parallel for
    for(size_t i = 0; i < cores; ++i) {
      boost::mt19937 subMt(seeds[i]);

      samplesPerCore.push_back(gps::rgps(draws, subMt, altRate, altWeight, noOfSnps));
    }

    std::vector<double> gpsSample;

    for(int i = 0; i < cores; ++i) {
      for(int j = 0; j < draws; ++j) {
        gpsSample.push_back(samplesPerCore[i][j]);
      }
    }

    for(auto i : gpsSample) std::cout << i << std::endl;

    } else {
      std::cout << desc << std::endl;
    }

  return 0;
}
