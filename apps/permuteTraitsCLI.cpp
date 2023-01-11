#include <iostream>
#include <gps.hpp>
#include <boost/program_options.hpp>
#include <rapidcsv.h>
#include <omp.h>
#include <PPEcdf.hpp>

namespace po = boost::program_options;
using namespace po;
using namespace gps;
using namespace rapidcsv;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  std::string inputFile;
  std::string outputFile;
  std::string colLabelA;
  std::string colLabelB;
  int perturbN = 0;
  double epsilonMultiple = 2.0;
  int cores;
  int draws;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("perturbN,p", po::value<int>(&perturbN), "No. of perturbation iterations")
    ("epsilonMultiple,e", po::value<double>(&epsilonMultiple), "Multiple of epsilon to use in perturbation procedure")
    ("cores,c", po::value<int>(&cores), "No. of cores")
    ("draws,n", po::value<int>(&draws), "No. of draws")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document input(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u;
    std::vector<double> v;

    try {
      u = input.GetColumn<double>(colLabelA);
    } catch(std::out_of_range stod){
      for(size_t i = 0; i < input.GetRowCount(); ++i){
        try {
          u.push_back(input.GetCell<double>(colLabelA, i));
        } catch(std::out_of_range stod) {
          u.push_back(1.0);
          //logOutput << traitA << "\t" << i+1 << std::endl;
        }
      }
    }

    try {
      v = input.GetColumn<double>(colLabelB);
    } catch(std::out_of_range stod){
      for(size_t i = 0; i < input.GetRowCount(); ++i){
        try {
          v.push_back(input.GetCell<double>(colLabelB, i));
        } catch(std::out_of_range stod) {
          v.push_back(1.0);
          //logOutput << traitB << "\t" << i+1 << std::endl;
        }
      }
    }

    for(size_t i = 0; i < perturbN; ++i) {
      u = perturbDuplicates_addEpsilon(u, epsilonMultiple);
      v = perturbDuplicates_addEpsilon(v, epsilonMultiple);
    }

    std::vector<std::vector<double>> gpsPermutations;

    int drawsPerCore = draws / cores;

    int remainingDraws = draws - drawsPerCore*cores;

    omp_set_num_threads(cores);

    #pragma omp parallel for
    for(int k = 0; k < cores; ++k) {
      gpsPermutations.push_back(permuteAndSampleGps(u, v, drawsPerCore, &PPEcdf::bivariatePPEcdf));
    }

    if(remainingDraws > 0) {
      gpsPermutations.push_back(permuteAndSampleGps(u, v, remainingDraws, &PPEcdf::bivariatePPEcdf));
    }

    std::stringstream stringOutput;

    stringOutput << "GPS" << std::endl;

    for(int i = 0; i < gpsPermutations.size(); ++i) {
      for(int j = 0; j < gpsPermutations[i].size(); ++j) {
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
