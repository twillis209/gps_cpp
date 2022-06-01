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

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("perturbFn,f", po::value<std::string>(&perturbFn), "Perturbation function")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::map<std::string, std::function<std::vector<double>(std::vector<double>)>> fnMap =
    {{"nextafter", perturbDuplicates},
     {"addEpsilon", perturbDuplicates_addEpsilon}
    };

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = data.GetColumn<double>(colLabelA);
    std::vector<double> v = data.GetColumn<double>(colLabelB);

    std::map<double, int> uFreqMap = returnFreqMap(u);
    std::map<double, int> vFreqMap = returnFreqMap(v);

    std::vector<double> uNoDup;
    std::vector<double> vNoDup;

    for(size_t i = 0; i < 99; ++i) {
      u = fnMap[perturbFn](u);
      v = fnMap[perturbFn](v);
    }

    uNoDup = u;
    vNoDup = v;

    /*
    for(size_t i = 0; i < u.size(); ++i) {
      if(u[i] != 1.0 && v[i] != 1.0) {
        uNoDup.push_back(u[i]);
        vNoDup.push_back(v[i]);
      }
    }
    */

    std::map<double, int> uNoDupFreqMap = returnFreqMap(uNoDup);
    std::map<double, int> vNoDupFreqMap = returnFreqMap(vNoDup);

    std::cout.precision(20);

    for(size_t i = 0; i < uNoDup.size(); ++i) {
      std::cout << u[i] << "\t" << uFreqMap[u[i]] << "\t" << uNoDup[i] << "\t" << uNoDupFreqMap[uNoDup[i]] << "\t" << v[i] << "\t" << vFreqMap[u[i]] << "\t" << vNoDup[i] << "\t" << vNoDupFreqMap[vNoDup[i]] << std::endl;
    }
  } else {
    std::cout << desc << std::endl;
  }

  return 0;
}
