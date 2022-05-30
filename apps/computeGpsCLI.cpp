#include <iostream>
#include <gps.hpp>
#include <boost/program_options.hpp>
#include <rapidcsv.h>
#include <omp.h>

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
  bool lwFlag = false;
  bool perturbFlag = false;
  int cores = 1;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("traitA,c", po::value<std::string>(&traitA), "Trait A")
    ("traitB,d", po::value<std::string>(&traitB), "Trait B")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
    ("lwFlag,l", po::bool_switch(&lwFlag), "Flag to use the fast bivariate ecdf algorithm from Langrene and Warin")
    ("lwFlag,l", po::bool_switch(&lwFlag), "Flag to use the fast bivariate ecdf algorithm from Langrene and Warin")
    ("perturbFlag,p", po::bool_switch(&perturbFlag), "Flag to use the fast bivariate ecdf algorithm from Langrene and Warin")
    ("cores,c", po::value<int>(&cores), "No. of cores")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = data.GetColumn<double>(colLabelA);
    std::vector<double> v = data.GetColumn<double>(colLabelB);

    if(perturbFlag) {
      u = perturbDuplicates(u);
      v = perturbDuplicates(v);
    }

    omp_set_num_threads(cores);

    double gps;

    if(lwFlag) {
      gps = gpsStat(u, v, &bivariateEcdfLW);
    } else {
      gps = gpsStat(u, v, &bivariateEcdfPar);
    }

    std::stringstream stringOutput;

    stringOutput << "Trait_A\tTrait_B\tGPS" << std::endl;

    stringOutput << traitA << "\t" << traitB << "\t" << gps << std::endl;

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
