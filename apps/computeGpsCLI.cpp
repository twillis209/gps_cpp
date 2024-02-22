#include <iostream>
#include <iomanip>
#include <gps.hpp>
#include <PPEcdf.hpp>
#include <CLI/CLI.hpp>
#include <boost/timer/timer.hpp>
#include <rapidcsv.h>
#include <omp.h>
#include <cstddef>

using namespace CLI;
using namespace gps;
using namespace rapidcsv;
using namespace std;

int main(int argc, const char* argv[]) {
  App app{"Compute GPS statistic"};

  string inputFile;
  string outputFile;
  string logFile;
  string timingFile;
  string traitA;
  string traitB;
  string colLabelA;
  string colLabelB;
  string ecdf = "pp";
  int cores = 1;
  string statistic = "gps";
  string weight = "gps";

  app.add_option("-i,--inputFile", inputFile, "Path to input file");
  app.add_option("-a,--colLabelA", colLabelA, "Label of column A");
  app.add_option("-b,--colLabelB", colLabelB, "Label of column B");
  app.add_option("-c,--traitA", traitA, "Trait A");
  app.add_option("-d,--traitB", traitB, "Trait B");
  app.add_option("-o,--outputFile", outputFile, "Path to output file");
  app.add_option("-j,--timingFile", timingFile, "Path to timing file");
  app.add_option("-g,--logFile", logFile, "Path to log file");
  app.add_option("-f,--ecdf", ecdf, "Specifies ecdf algorithm: 'naive' or 'pp'");
  app.add_option("-n,--cores", cores, "No. of cores");
  app.add_option("-s,--statistic", statistic, "Statistic to compute");
  app.add_option("-w,--weight", weight, "Weight function to use");

  CLI11_PARSE(app, argc, argv);

  if(statistic != "gps" && statistic != "mean") {
    throw invalid_argument("Unrecognised statistic argument");
  }


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

    int n = u.size();

    omp_set_num_threads(cores);

    cout << "Computing the GPS statistic..." << endl;

    boost::timer::cpu_timer timer;

    gps = statFun(u, v, ecdfFun, weightFun);

    if(!timingFile.empty()) {
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

  return 0;
}
