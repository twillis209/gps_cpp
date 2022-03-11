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
  std::string traitA;
  std::string traitB;
  std::string colLabelA;
  std::string colLabelB;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
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

    std::vector<double> bivariateEcdf = bivariateEcdfLW(uNoDup, vNoDup);

    std::vector<double> uEcdf = ecdf(uNoDup);
    std::vector<double> vEcdf = ecdf(vNoDup);

    std::stringstream stringOutput;

    stringOutput << "u\tv\tF(u)\tF(v)\tF(u,v)" << std::endl;

    for(int i = 0; i < uNoDup.size(); i++) {
      stringOutput << uNoDup[i] << "\t" << vNoDup[i] << "\t" << uEcdf[i] << "\t" << vEcdf[i] << "\t" << bivariateEcdf[i] << std::endl;

    }

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
