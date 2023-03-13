#include <iostream>
#include <gps.hpp>
#include <PPEcdf.hpp>
#include <boost/program_options.hpp>
#include <rapidcsv.h>
#include <omp.h>
#include <math.h>
#include <cstddef>

namespace po = boost::program_options;
using namespace po;
using namespace gps;
using namespace rapidcsv;
using namespace std;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  string inputFile;
  string outputFile;
  string colLabelA;
  string colLabelB;
  string ecdfArg = "naive";
  int cores = 1;
  double epsilonMultiple = 2.0;
  bool deduplicateFlag = false;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<string>(&colLabelB), "Label of column B")
    ("ecdfArg,f", po::value<string>(&ecdfArg), "Specifies ecdf algorithm: \"naive\" or \"pp\"")
    ("outputFile,o", po::value<string>(&outputFile), "Path to output file")
    ("cores,n", po::value<int>(&cores), "No. of cores")
    ("deduplicate,d", po::bool_switch(&deduplicateFlag), "Remove duplicate values")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    vector<double> u = data.GetColumn<double>(colLabelA);
    vector<double> v = data.GetColumn<double>(colLabelB);

    map<double, int> freqMapU = returnFreqMap(u);
    map<double, int> freqMapV = returnFreqMap(v);

    int n = u.size();

    omp_set_num_threads(cores);

    vector<double> bivariateEcdf;

    if(ecdfArg == "naive") {
      bivariateEcdf = bivariateEcdfPar(u, v);
    } else if(ecdfArg == "pp") {
      bivariateEcdf = PPEcdf::bivariatePPEcdf(u, v);
    } else {
      cout << "Invalid ecdfArg argument, using naive algorithm" << endl;
      bivariateEcdf = bivariateEcdfPar(u, v);
    }

    vector<double> uEcdf = ecdf(u);
    vector<double> vEcdf = ecdf(v);

    stringstream stringOutput;

    stringOutput << "u\tv\tF_u\tF_v\tF_uv\tnum\tdenom\tmaximand" << endl;

    for(int i = 0; i < u.size(); i++) {
      double denom = sqrt(uEcdf[i]*vEcdf[i] - pow(uEcdf[i], 2.)*pow(vEcdf[i], 2.));

      double numerator = abs(bivariateEcdf[i] - uEcdf[i]*vEcdf[i]);

      double maximand = sqrt((double) u.size() / log((double) u.size()))*numerator/denom;

      stringOutput << u[i] << "\t" << v[i] << "\t" << uEcdf[i] << "\t" << vEcdf[i] << "\t" << bivariateEcdf[i] << "\t" << numerator << "\t" << denom << "\t" << maximand << endl;
    }

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);
  } else {
      cout << desc << endl;
  }

  return 0;
}
