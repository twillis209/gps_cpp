#include <iostream>
#include <gps.hpp>
#include <boost/program_options.hpp>
#include <rapidcsv.h>
#include <omp.h>

namespace po = boost::program_options;
using namespace po;
using namespace gps;
using namespace rapidcsv;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  std::string inputFile;
  std::string outputFile;
  std::string columnA;
  std::string columnB;
  int cores;
  int drawsPerCore;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
    ("columnA,a", po::value<std::string>(&columnA), "Label of column A")
    ("columnB,b", po::value<std::string>(&columnB), "Label of column B")
    ("cores,c", po::value<int>(&cores), "No. of cores")
    ("draws,n", po::value<int>(&drawsPerCore), "No. of draws per core")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document input(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = input.GetColumn<double>(columnA);
    std::vector<double> v = input.GetColumn<double>(columnB);

    std::vector<double> uNoDup = perturbDuplicates(u);
    std::vector<double> vNoDup = perturbDuplicates(v);

    std::vector<std::vector<double>> gpsPermutations;

    // TODO remove me
    #pragma omp parallel
    {
      std::cout << omp_get_num_threads() << std::endl;
    }


    #pragma omp parallel for
    for(int k = 0; k < cores; ++k) {
      gpsPermutations.push_back(permuteAndSampleGps(uNoDup, vNoDup, drawsPerCore));
    }

    std::stringstream stringOutput;

    stringOutput << "GPS" << std::endl;

    for(int i = 0; i < cores; ++i) {
      for(int j = 0; j < drawsPerCore; ++j) {
        stringOutput << gpsPermutations[i][j] << std::endl;
      }
    }

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
