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

    stringOutput << "Trait A\tTrait B\tGPS" << std::endl;

    // Ignore variant column
    for(size_t i = 1; i < noOfColumns-1; ++i) {
      for(size_t j = i+1; j < noOfColumns; ++j) {
        std::vector<double> u = data.GetColumn<double>(i);
        std::vector<double> v = data.GetColumn<double>(j);

        std::vector<double> uNoDup = perturbDuplicates(u);
        std::vector<double> vNoDup = perturbDuplicates(v);

        double gps = gpsStat(uNoDup, vNoDup);

        stringOutput << data.GetColumnName(i) << '\t' << data.GetColumnName(j) << '\t' << gps << std::endl;
      }
    }

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
