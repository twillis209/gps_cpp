#include <iostream>
#include <gps.hpp>
#include <boost/program_options.hpp>
#include <rapidcsv.h>
#include <omp.h>
#include <math.h>

namespace po = boost::program_options;
using namespace po;
using namespace gps;
using namespace rapidcsv;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  std::string inputFile;
  std::string outputFile;
  std::string colLabelA;
  std::string colLabelB;
  int cores = 1;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
    ("cores,n", po::value<int>(&cores), "No. of cores")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = data.GetColumn<double>(colLabelA);
    std::vector<double> v = data.GetColumn<double>(colLabelB);

    omp_set_num_threads(cores);

    std::vector<double> bivariateEcdf = bivariateEcdfPar(u, v);

    std::vector<double> uEcdf = ecdf(u);
    std::vector<double> vEcdf = ecdf(v);

    std::stringstream stringOutput;

    stringOutput << "u\tv\tF_u\tF_v\tF_uv\tnum\tdenom\tmaximand" << std::endl;

    for(int i = 0; i < u.size(); i++) {
      double denom = sqrt(uEcdf[i]*vEcdf[i] - pow(uEcdf[i], 2.)*pow(vEcdf[i], 2.));

      double numerator = abs(bivariateEcdf[i] - uEcdf[i]*vEcdf[i]);

      double maximand = sqrt((double) u.size() / log((double) u.size()))*numerator/denom;

      stringOutput << u[i] << "\t" << v[i] << "\t" << uEcdf[i] << "\t" << vEcdf[i] << "\t" << bivariateEcdf[i] << "\t" << numerator << "\t" << denom << "\t" << maximand << std::endl;
    }

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
