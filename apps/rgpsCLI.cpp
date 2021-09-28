#include <iostream>
#include <gps.hpp>
#include <boost/program_options.hpp>
#include <boost/random.hpp>

namespace po = boost::program_options;
using namespace po;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  size_t n;
  size_t noOfSnps;
  double altWeight;
  double altRate;
  unsigned int seed;

  desc.add_options()
    ("help", "Print help message")
    ("size,n", po::value<size_t>(&n), "No. of variates")
    ("noOfSnps,p", po::value<size_t>(&noOfSnps)->default_value(1e4), "No. of SNPs")
    ("altWeight,w", po::value<double>(&altWeight)->default_value(0), "Mixing weight for non-standard exp. component")
    ("altRate,r", po::value<double>(&altRate)->default_value(5), "Rate parameter for non-standard exp. component")
    ("seed,s", po::value<unsigned int>(&seed)->default_value(42u), "Seed for PRNG")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // TODO to parallelise, could generate uniform ints from the seed to, in turn, seed threads like in the permute CLI
  if(vm.count("size")) {
    boost::mt19937 mt(seed);


    std::vector<double> gpsSample = gps::rgps(n, mt, altRate, altWeight, noOfSnps);

    for(auto i : gpsSample) std::cout << i << std::endl;
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
