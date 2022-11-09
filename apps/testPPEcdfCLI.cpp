#include <PPEcdf.hpp>
#include <rapidcsv.h>
#include <gps.hpp>
#include <boost/program_options.hpp>
#include <omp.h>
#include <iostream>
#include <iomanip>

namespace po = boost::program_options;
using namespace po;
using namespace PPEcdf;
using namespace rapidcsv;
using namespace gps;


int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  std::string inputFile;
  std::string outputFile;
  std::string colLabelA;
  std::string colLabelB;
  int perturbN = 0;
  double epsilonMultiple = 2.0;
  int cores = 1;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
    ("perturbN,p", po::value<int>(&perturbN), "No. of perturbation iterations")
    ("epsilonMultiple,e", po::value<double>(&epsilonMultiple), "Multiple of epsilon to use in perturbation procedure")
    ("cores,n", po::value<int>(&cores), "No. of cores")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams(','));

    std::vector<double> u;
    std::vector<double> v;

    try {
      u = data.GetColumn<double>(colLabelA);
    } catch(std::out_of_range stod){
      for(size_t i = 0; i < data.GetRowCount(); ++i){
        try {
          u.push_back(data.GetCell<double>(colLabelA, i));
        } catch(std::out_of_range stod) {
          u.push_back(1.0);
        }
      }
    }

    try {
      v = data.GetColumn<double>(colLabelB);
    } catch(std::out_of_range stod){
      for(size_t i = 0; i < data.GetRowCount(); ++i){
        try {
          v.push_back(data.GetCell<double>(colLabelB, i));
        } catch(std::out_of_range stod) {
          v.push_back(1.0);
        }
      }
    }

    omp_set_num_threads(cores);

    std::vector<double> ecdf = bivariatePPEcdf(u, v);
    
    std::vector<double> ex_ecdf = bivariateEcdfPar(u, v);

    std::stringstream stringOutput;

    for(size_t i = 0; i < u.size(); i++) {
      stringOutput << u[i] << "\t" << v[i] << "\t" << ecdf[i] << "\t" << ex_ecdf[i] << std::endl;
    }

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);

    } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
