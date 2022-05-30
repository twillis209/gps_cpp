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

  std::string inputFile;
  std::string outputFile;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    size_t noOfColumns = data.GetColumnCount();

    std::stringstream stringOutput;

    stringOutput << "Trait_A\tTrait_B\tGPS" << std::endl;

    std::vector<double> u = data.GetColumn<double>(4);
    std::vector<double> uNoDup = perturbDuplicates(u);

    // Ignore variant column
    for(size_t i = 6; i < 27; ++i) {
      //std::cout << "pid " << " i " << i << std::endl;
        std::vector<double> v = data.GetColumn<double>(i);
        std::vector<double> vNoDup = perturbDuplicates(v);

        double gps = gpsStat(uNoDup, vNoDup, &bivariateEcdfLW);

        stringOutput << "pid" << '\t' << data.GetColumnName(i) << '\t' << gps << std::endl;
     }

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
