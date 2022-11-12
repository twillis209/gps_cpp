#include <iostream>
#include <gps.hpp>
#include <PPEcdf.hpp>
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
  std::string ecdfArg = "naive";
  int cores = 1;
  int perturbN = 0;
  double epsilonMultiple = 2.0;
  bool deduplicateFlag = false;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("ecdfArg,f", po::value<std::string>(&ecdfArg), "Specifies ecdf algorithm: \"naive\", \"pp\", or \"lw\"")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
    ("cores,n", po::value<int>(&cores), "No. of cores")
    ("perturbN,p", po::value<int>(&perturbN), "No. of perturbation iterations")
    ("epsilonMultiple,e", po::value<double>(&epsilonMultiple), "Multiple of epsilon to use in perturbation procedure")
    ("deduplicate,d", po::bool_switch(&deduplicateFlag), "Remove duplicate values")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = data.GetColumn<double>(colLabelA);
    std::vector<double> v = data.GetColumn<double>(colLabelB);

    if(perturbN > 0) {

      std::cout << "Perturbing..." << std::endl;

      for(size_t i = 0; i < perturbN; ++i) {
        u = perturbDuplicates_addEpsilon(u, epsilonMultiple);
        v = perturbDuplicates_addEpsilon(v, epsilonMultiple);
      }
    }

    std::map<double, int> freqMapU = returnFreqMap(u);
    std::map<double, int> freqMapV = returnFreqMap(v);

    int n = u.size();

    if(deduplicateFlag) {
      std::cout << "Length of u vector before deletion: " << n << std::endl;

      // Delete any values we couldn't perturb away from being duplicates
      for(size_t i = 0; i < n; ++i) {

        //std::cout << i << std::endl;

        if(freqMapU[u[i]] > 1 || freqMapV[v[i]] > 1 || u[i] > 1.0 || v[i] > 1.0) {
          std::cout << "Erasing " << i << " " << std::endl;
          u.erase(u.begin()+i);
          v.erase(v.begin()+i);
        }
      }

      std::cout << "Length of u vector after deletion: " << u.size() << std::endl;
    }

    omp_set_num_threads(cores);

    std::vector<double> bivariateEcdf;

    if(ecdfArg == "naive") {
      bivariateEcdf = bivariateEcdfPar(u, v);
    } else if(ecdfArg == "pp") {
      bivariateEcdf = PPEcdf::bivariatePPEcdf(u, v);
    } else if(ecdfArg == "lw") {
      bivariateEcdf = bivariateEcdfLW(u, v);
    } else {
      std::cout << "Invalid ecdfArg argument, using naive algorithm" << std::endl;
      bivariateEcdf = bivariateEcdfPar(u, v);
    }

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
