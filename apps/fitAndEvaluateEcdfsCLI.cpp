#include <iostream>
#include <gps.hpp>
#include <PPEcdf.hpp>
#include <CLI/CLI.hpp>
#include <rapidcsv.h>
#include <omp.h>
#include <math.h>
#include <cstddef>

using namespace CLI;
using namespace gps;
using namespace rapidcsv;
using namespace std;
 
int main(int argc, const char* argv[]) {
  App app{"Compute and evaluate ecdf"};

  string inputFile;
  string outputFile;
  string colLabelA;
  string colLabelB;
  string ecdfArg = "naive";
  int cores = 1;

  app.add_option("-i,--inputFile", inputFile, "Path to input file")->required()->check(CLI::ExistingFile);
  app.add_option("-o,--outputFile", outputFile, "Path to output file")->required();
  app.add_option("-a,--colLabelA", colLabelA, "Label of column A")->required();
  app.add_option("-b,--colLabelB", colLabelB, "Label of column B")->required();
  app.add_option("-f,--ecdf", ecdfArg, "Specifies ecdf algorithm: 'naive' or 'pp'");
  app.add_option("-n,--cores", cores, "No. of cores");

  CLI11_PARSE(app, argc, argv);

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

  return 0;
}
