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
  std::string traitA;
  std::string traitB;
  std::string outputFile;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("traitA,a", po::value<std::string>(&inputFile), "Label of trait A column")
    ("traitB,b", po::value<std::string>(&inputFile), "Label of trait B column")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    std::cout << "Before data" << std::endl;

    // TODO Problem is here
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::cout << "After data" << std::endl;

    size_t traitAIdx = data.GetColumnIdx(traitA);
    size_t traitBIdx = data.GetColumnIdx(traitB);

    std::stringstream stringOutput;

    stringOutput << "Trait_A\tTrait_B\tGPS" << std::endl;

    std::vector<double> u = data.GetColumn<double>(traitA);
    std::vector<double> v = data.GetColumn<double>(traitB);

    std::vector<double> uNoDup = perturbDuplicates(u);
    std::vector<double> vNoDup = perturbDuplicates(v);

    double gps = gpsStat(uNoDup, vNoDup);

    stringOutput << traitA << '\t' << traitB << '\t' << gps << std::endl;

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
