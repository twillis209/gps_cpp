#include <iostream>
#include <stdexcept>
#include <gps.hpp>
#include <CLI/CLI.hpp>
#include <rapidcsv.h>
#include <omp.h>
#include <PPEcdf.hpp>
#include <cstddef>

using namespace CLI;
using namespace gps;
using namespace rapidcsv;
using namespace std;

int main(int argc, const char* argv[]) {
  App app{"Generate permuted realisations of the GPS statistic"};

  string inputFile;
  string outputFile;
  string colLabelA;
  string colLabelB;
  int cores = 1;
  int draws;
  string statistic = "gps";
  string weight = "gps";

  app.add_option("-i,--inputFile", inputFile, "Path to input file");
  app.add_option("-o,--outputFile", outputFile, "Path to output file");
  app.add_option("-a,--colLabelA", colLabelA, "Label of column A");
  app.add_option("-b,--colLabelB", colLabelB, "Label of column B");
  app.add_option("-n,--cores", cores, "No. of cores");
  app.add_option("-d,--draws", cores, "No. of GPS realisations to generate");

  CLI11_PARSE(app, argc, argv);

  if(statistic != "gps" && statistic != "mean") {
    throw invalid_argument("Unrecognised statistic argument");
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

  Document input(inputFile, LabelParams(), SeparatorParams('\t'));

  vector<double> u;
  vector<double> v;

  try {
    u = input.GetColumn<double>(colLabelA);
  } catch(out_of_range stod){
    for(size_t i = 0; i < input.GetRowCount(); ++i){
      try {
        u.push_back(input.GetCell<double>(colLabelA, i));
      } catch(out_of_range stod) {
        u.push_back(1.0);
      }
    }
  }

  try {
    v = input.GetColumn<double>(colLabelB);
  } catch(out_of_range stod){
    for(size_t i = 0; i < input.GetRowCount(); ++i){
      try {
        v.push_back(input.GetCell<double>(colLabelB, i));
      } catch(out_of_range stod) {
        v.push_back(1.0);
      }
    }
  }

  vector<vector<double>> gpsPermutations;

  int drawsPerCore = draws / cores;

  int remainingDraws = draws - drawsPerCore*cores;

  omp_set_num_threads(cores);

  #pragma omp parallel for
  for(int k = 0; k < cores; ++k) {
      gpsPermutations.push_back(permuteAndSampleStat(u, v, drawsPerCore, statFun, PPEcdf::bivariatePPEcdf, weightFun));
  }

  if(remainingDraws > 0) {
    gpsPermutations.push_back(permuteAndSampleStat(u, v, remainingDraws, statFun, PPEcdf::bivariatePPEcdf, weightFun));
  }

  stringstream stringOutput;

  stringOutput << "GPS" << endl;

  for(int i = 0; i < gpsPermutations.size(); ++i) {
    for(int j = 0; j < gpsPermutations[i].size(); ++j) {
      stringOutput << gpsPermutations[i][j] << endl;
    }
  }

  Document output(stringOutput, LabelParams(), SeparatorParams('\t'));

  output.Save(outputFile);

  return 0;
}
