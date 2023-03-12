#include <iostream>
#include <iomanip>
#include <gps.hpp>
#include <PPEcdf.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <rapidcsv.h>
#include <omp.h>

namespace po = boost::program_options;
using namespace po;
using namespace gps;
using namespace rapidcsv;
using namespace std;

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  string inputFile;
  string outputFile;
  string logFile;
  string timingFile;
  string traitA;
  string traitB;
  string colLabelA;
  string colLabelB;
  string ecdf = "naive";
  double epsilonMultiple = 2.0;
  int cores = 1;
  bool deduplicateFlag = false;
  string statistic = "gps";
  string weight = "gps";

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<string>(&colLabelB), "Label of column B")
    ("traitA,c", po::value<string>(&traitA), "Trait A")
    ("traitB,d", po::value<string>(&traitB), "Trait B")
    ("outputFile,o", po::value<string>(&outputFile), "Path to output file")
    ("timingFile,j", po::value<string>(&timingFile), "Path to timing file")
    ("logFile,g", po::value<string>(&logFile), "Path to log file")
    ("ecdf,f", po::value<string>(&ecdf), "Specifies ecdf algorithm: \"naive\" or \"pp\"")
    ("cores,n", po::value<int>(&cores), "No. of cores")
    ("statistic,s", po::value<string>(&statistic), "Statistic to compute")
    ("weight,w", po::value<string>(&weight), "Weight function to use")
    ;

  if(statistic != "gps" && statistic != "mean") {
    throw invalid_argument("Unrecognised statistic argument");
  }

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {

    function<vector<double>(const vector<double>&, const vector<double>&)> ecdfFun;

    if(ecdf == "naive") {
      ecdfFun = bivariateEcdfPar;
    } else if(ecdf == "pp") {
      ecdfFun = PPEcdf::bivariatePPEcdf;
    } else {
      cout << "Invalid ecdf argument, using naive algorithm" << endl;
      ecdfFun = bivariateEcdfPar;
    }

    function<double (vector<double>, vector<double>, function<vector<double>(const vector<double>&, const vector<double>&)>, function<double (const double&, const double&, const double&)>)> statFun;

    if(statistic == "gps") {
      statFun = gpsStat;
    } else if(statistic == "mean") {
      statFun = meanStat;
    } else {
      throw invalid_argument("Unrecognised statistic argument");
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

    stringstream logOutput;

    logOutput << "Trait\tRow" << endl;

    vector<double> u;
    vector<double> v;

    try {
      u = data.GetColumn<double>(colLabelA);
    } catch(out_of_range stod){
      for(size_t i = 0; i < data.GetRowCount(); ++i){
        try {
          u.push_back(data.GetCell<double>(colLabelA, i));
        } catch(out_of_range stod) {
          u.push_back(1.0);
          logOutput << traitA << "\t" << i+1 << endl;
        }
      }
    }

    try {
      v = data.GetColumn<double>(colLabelB);
    } catch(out_of_range stod){
      for(size_t i = 0; i < data.GetRowCount(); ++i){
        try {
          v.push_back(data.GetCell<double>(colLabelB, i));
        } catch(out_of_range stod) {
          v.push_back(1.0);
          logOutput << traitB << "\t" << i+1 << endl;
        }
      }
    }

    double gps;

    cout.precision(20);

    map<double, int> freqMapU = returnFreqMap(u);
    map<double, int> freqMapV = returnFreqMap(v);

    int n = u.size();

    omp_set_num_threads(cores);

    cout << "Computing the GPS statistic..." << endl;

    boost::timer::cpu_timer timer;

    gps = statFun(u, v, ecdfFun, weightFun);

    if(vm.count("timingFile")) {
      timer.stop();
      ofstream timingOut(timingFile);
      timingOut << timer.format();
      timingOut.close();
    }

    stringstream stringOutput;

    stringOutput << "Trait_A\tTrait_B\tGPS" << endl;

    stringOutput << traitA << "\t" << traitB << "\t" << gps << endl;

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);

    if(!logFile.empty()) {
      Document log(logOutput, LabelParams(), SeparatorParams('\t'));
      log.Save(logFile);
    }
  } else {
      cout << desc << endl;
  }

  return 0;
}
