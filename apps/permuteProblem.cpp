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
  std::string aOutputFile;
  std::string bOutputFile;
  std::string colLabelA;
  std::string colLabelB;
  std::string perturbFn;
  double epsilonMultiple = 2.0;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("epsilonMultiple,e", po::value<double>(&epsilonMultiple), "Multiple of epsilon to use in perturbation procedure")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = data.GetColumn<double>(colLabelA);
    std::vector<double> v = data.GetColumn<double>(colLabelB);

    std::map<double, int> uFreqMap = returnFreqMap(u);
    std::map<double, int> vFreqMap = returnFreqMap(v);

    std::vector<double> uNoDup = perturbDuplicates_addEpsilon(u, epsilonMultiple);
    std::vector<double> vNoDup = perturbDuplicates_addEpsilon(v, epsilonMultiple);

    for(size_t i = 0; i < 99; ++i) {
      uNoDup = perturbDuplicates_addEpsilon(uNoDup, epsilonMultiple);
      vNoDup = perturbDuplicates_addEpsilon(vNoDup, epsilonMultiple);
    }

    std::map<double, int> uNoDupFreqMap = returnFreqMap(uNoDup);
    std::map<double, int> vNoDupFreqMap = returnFreqMap(vNoDup);

    std::cout.precision(20);

    for(size_t i = 0; i < u.size(); ++i) {
      std::cout << u[i] << "\t" << uFreqMap[u[i]] << "\t" << uNoDup[i] << "\t" << uNoDupFreqMap[uNoDup[i]] << "\t" << v[i] << "\t" << vFreqMap[u[i]] << "\t" << vNoDup[i] << "\t" << vNoDupFreqMap[vNoDup[i]] << std::endl;
    }
  } else {
    std::cout << desc << std::endl;
  }

  return 0;
}
