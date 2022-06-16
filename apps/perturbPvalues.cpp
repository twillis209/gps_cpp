#include <iostream>
#include <map>
#include <gps.hpp>
#include <rapidcsv.h>
#include <random>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace po;
using namespace gps;
using namespace rapidcsv;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  std::string inputFile;
  std::string colLabelA;
  std::string colLabelB;
  std::string traitA;
  std::string traitB;
  double epsilonMultiple = 2.0;
  int perturbN = 100;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("traitA,c", po::value<std::string>(&traitA), "Trait A")
    ("traitB,d", po::value<std::string>(&traitB), "Trait B")
    ("epsilonMultiple,e", po::value<double>(&epsilonMultiple), "Multiple of epsilon to use in perturbation procedure")
    ("perturbN,p", po::value<int>(&perturbN), "No. of perturbation iterations")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = data.GetColumn<double>(colLabelA);
    std::vector<double> v = data.GetColumn<double>(colLabelB);

    try {
      u = data.GetColumn<double>(colLabelA);
    } catch(std::out_of_range stod){
      for(size_t i = 0; i < data.GetRowCount(); ++i){
        try {
          u.push_back(data.GetCell<double>(colLabelA, i));
        } catch(std::out_of_range stod) {
          u.push_back(1.0);
        }
      }
    }

    try {
      v = data.GetColumn<double>(colLabelB);
    } catch(std::out_of_range stod){
      for(size_t i = 0; i < data.GetRowCount(); ++i){
        try {
          v.push_back(data.GetCell<double>(colLabelB, i));
        } catch(std::out_of_range stod) {
          v.push_back(1.0);
        }
      }
    }

    for(size_t i = 0; i < perturbN; ++i) {
      u = perturbDuplicates_addEpsilon(u, epsilonMultiple);
      v = perturbDuplicates_addEpsilon(v, epsilonMultiple);
    }

    for(size_t i = 0; i < u.size(); ++i) {
      if(u[i] > 1.0 || v[i] > 1.0) {
        u.erase(u.end()-(i+1));
        v.erase(v.end()-(i+1));
      }
    }

    std::cout.precision(20);

    std::cout << traitA << "\t" << traitB << std::endl;

    for(size_t i = 0; i < u.size(); ++i) {
      std::cout << u[i] << "\t" << v[i] << std::endl;
    }
  } else {
    std::cout << desc << std::endl;
  }

  return 0;
}
