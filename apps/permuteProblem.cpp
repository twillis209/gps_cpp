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
  std::string colLabelA;
  std::string colLabelB;

  desc.add_options()
    ("help", "Print help message")

    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
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

    std::cout.precision(35);

    std::cout << "Trait_A\tPert_trait_A\tTrait_B\tPert_trait_B" << std::endl;

    for(size_t i = 0; i < uNoDup.size(); ++i) {
      std::cout << u[i] << "\t" << uNoDup[i] << "\t" << v[i] << "\t" << vNoDup[i] << std::endl;
    }

  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
