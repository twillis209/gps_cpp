#include <iostream>
#include <gps.hpp>
#include <boost/program_options.hpp>
#include <rapidcsv.h>

namespace po = boost::program_options;
using namespace po;
using namespace gps;
using namespace rapidcsv;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  int cores;
  int drawsPerCore;
  std::string inputFile;
  std::string colLabelA;
  std::string colLabelB;

  desc.add_options()
    ("help", "Print help message")

    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("cores,c", po::value<int>(&cores), "No. of cores")
    ("draws,n", po::value<int>(&drawsPerCore), "No. of draws per core")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile);

    std::vector<double> u = data.GetColumn<double>(colLabelA);
    std::vector<double> v = data.GetColumn<double>(colLabelB);

    std::vector<double> uNoDup = perturbDuplicates(u);
    std::vector<double> vNoDup = perturbDuplicates(v);

    std::vector<std::vector<double>> gpsPermutations;

    #pragma omp parallel for
    for(int i = 0; i < cores; ++i) {
      gpsPermutations.push_back(permuteAndSampleGps(uNoDup, vNoDup, drawsPerCore, &bivariateEcdfLW));
    }

    std::vector<double> gpsResults;

    for(int i = 0; i < cores; ++i) {
      for(int j = 0; j < drawsPerCore; ++j) {
        gpsResults.push_back(gpsPermutations[i][j]);
      }
    }

    for(auto i : gpsResults) std::cout << i << std::endl;

  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
