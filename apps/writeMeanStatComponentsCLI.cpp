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
using namespace std;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  string inputFile;
  string outputFile;
  string colLabelA;
  string colLabelB;
  string ecdfArg = "naive";
  string weight = "gps";
  int cores = 1;
  int perturbN = 0;
  double epsilonMultiple = 2.0;
  bool deduplicateFlag = false;

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<string>(&colLabelB), "Label of column B")
    ("ecdfArg,f", po::value<string>(&ecdfArg), "Specifies ecdf algorithm: \"naive\", \"pp\", or \"lw\"")
    ("weight,w", po::value<string>(&weight), "Weight function to use")
    ("outputFile,o", po::value<string>(&outputFile), "Path to output file")
    ("cores,n", po::value<int>(&cores), "No. of cores")
    ("perturbN,p", po::value<int>(&perturbN), "No. of perturbation iterations")
    ("epsilonMultiple,e", po::value<double>(&epsilonMultiple), "Multiple of epsilon to use in perturbation procedure")
    ("deduplicate,d", po::bool_switch(&deduplicateFlag), "Remove duplicate values")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {

    function<vector<double>(const vector<double>&, const vector<double>&)> ecdfFun;

    if(ecdfArg == "naive") {
      ecdfFun = bivariateEcdfPar;
    } else if(ecdfArg == "pp") {
      ecdfFun = PPEcdf::bivariatePPEcdf;
    } else if(ecdfArg == "lw") {
      ecdfFun = bivariateEcdfLW;
    } else {
      cout << "Invalid ecdf argument, using naive algorithm" << endl;
      ecdfFun = bivariateEcdfPar;
    }

    function<double (double, double, double)> weightFun;

    if(weight == "gps") {
      weightFun = gpsWeight;
    } else if(weight == "pseudoAD") {
      weightFun = pseudoADWeight;
    } else {
      throw invalid_argument("Unrecognised weight argument");
    }

    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    vector<double> u = data.GetColumn<double>(colLabelA);
    vector<double> v = data.GetColumn<double>(colLabelB);

    if(perturbN > 0) {

      cout << "Perturbing..." << endl;

      for(size_t i = 0; i < perturbN; ++i) {
        u = perturbDuplicates_addEpsilon(u, epsilonMultiple);
        v = perturbDuplicates_addEpsilon(v, epsilonMultiple);
      }
    }

    map<double, int> freqMapU = returnFreqMap(u);
    map<double, int> freqMapV = returnFreqMap(v);

    int n = u.size();

    if(deduplicateFlag) {
      cout << "Length of u vector before deletion: " << n << endl;

      // Delete any values we couldn't perturb away from being duplicates
      for(size_t i = 0; i < n; ++i) {

        //cout << i << endl;

        if(freqMapU[u[i]] > 1 || freqMapV[v[i]] > 1 || u[i] > 1.0 || v[i] > 1.0) {
          cout << "Erasing " << i << " " << endl;
          u.erase(u.begin()+i);
          v.erase(v.begin()+i);
        }
      }

      cout << "Length of u vector after deletion: " << u.size() << endl;
    }

    omp_set_num_threads(cores);

    vector<double> uEcdf = ecdf(u);
    vector<double> vEcdf = ecdf(v);
    vector<double> uvEcdf = ecdfFun(u, v);

    stringstream stringOutput;

    stringOutput << "u\tv\tF_u\tF_v\tF_uv\tnum\tweight" << endl;

    for(int i = 0; i < u.size(); i++) {
      double weight = weightFun(uEcdf[i], vEcdf[i], uvEcdf[i]);

      double numerator = pow(uvEcdf[i] - uEcdf[i]*vEcdf[i], 2.);

      stringOutput << u[i] << "\t" << v[i] << "\t" << uEcdf[i] << "\t" << vEcdf[i] << "\t" << uvEcdf[i] << "\t" << numerator << "\t" << weight << endl;
    }

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);
  } else {
      cout << desc << endl;
  }

  return 0;
}
