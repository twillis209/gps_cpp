#include <iostream>
#include <gps.hpp>
#include <rapidcsv.h>
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

  desc.add_options()
    ("help", "Print help message")

    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("aOutputFile,c", po::value<std::string>(&aOutputFile), "Path to A output file")
    ("bOutputFile,d", po::value<std::string>(&bOutputFile), "Path to B output file")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = data.GetColumn<double>(colLabelA);
    std::vector<double> v = data.GetColumn<double>(colLabelB);

    std::vector<double> uNoDup = perturbDuplicates(u);
    std::vector<double> vNoDup = perturbDuplicates(v);

    std::map<double, int> uNoDupFreqMap = returnFreqMap(uNoDup);
    std::map<double, int> vNoDupFreqMap = returnFreqMap(vNoDup);

    std::stringstream aStringOutput;

    for(size_t i = 0; i < uNoDup.size(); ++i) {
      aStringOutput << uNoDup[i] << "\t" << uNoDupFreqMap[uNoDup[i]] << std::endl;
    }

    Document aOutput(aStringOutput, LabelParams(), SeparatorParams('\t'));
    aOutput.Save(aOutputFile);

    std::stringstream bStringOutput;

    for(size_t i = 0; i < vNoDup.size(); ++i) {
      bStringOutput << vNoDup[i] << "\t" << vNoDupFreqMap[vNoDup[i]] << std::endl;
    }

    Document bOutput(bStringOutput, LabelParams(), SeparatorParams('\t'));
    bOutput.Save(bOutputFile);
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
