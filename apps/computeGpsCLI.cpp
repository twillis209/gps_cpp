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

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  std::string inputFile;
  std::string outputFile;
  std::string logFile;
  std::string perturbedFile;
  std::string timingFile;
  std::string traitA;
  std::string traitB;
  std::string colLabelA;
  std::string colLabelB;
  std::string ecdf = "naive";
  int perturbN = 0;
  double epsilonMultiple = 2.0;
  int cores = 1;
  bool deduplicateFlag = false;
  std::string statistic = "gps";

  desc.add_options()
    ("help", "Print help message")
    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("traitA,c", po::value<std::string>(&traitA), "Trait A")
    ("traitB,d", po::value<std::string>(&traitB), "Trait B")
    ("outputFile,o", po::value<std::string>(&outputFile), "Path to output file")
    ("timingFile,j", po::value<std::string>(&timingFile), "Path to timing file")
    ("logFile,g", po::value<std::string>(&logFile), "Path to log file")
    ("perturbedFile,t", po::value<std::string>(&perturbedFile), "Path to file containing perturbed u and v vectors")
    ("ecdf,f", po::value<std::string>(&ecdf), "Specifies ecdf algorithm: \"naive\", \"pp\", or \"lw\"")
    ("perturbN,p", po::value<int>(&perturbN), "No. of perturbation iterations")
    ("epsilonMultiple,e", po::value<double>(&epsilonMultiple), "Multiple of epsilon to use in perturbation procedure")
    ("cores,n", po::value<int>(&cores), "No. of cores")
    ;

  if(statistic != "gps" && statistic != "mean") {
    throw std::invalid_argument("Unrecognised statistic argument");
  }

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::stringstream logOutput;

    logOutput << "Trait\tRow" << std::endl;

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
          logOutput << traitA << "\t" << i+1 << std::endl;
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
          logOutput << traitB << "\t" << i+1 << std::endl;
        }
      }
    }

    double gps;

    std::cout.precision(20);

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
          if(freqMapU[u[i]] > 1 || freqMapV[v[i]] > 1 || u[i] > 1.0 || v[i] > 1.0) {
            u.erase(u.end()-(i+1));
            v.erase(v.end()-(i+1));
          }
        }

        std::cout << "Length of u vector after deletion: " << u.size() << std::endl;
      }

      if(!perturbedFile.empty()) {

        std::stringstream perturbedOutput;

        perturbedOutput << traitA << "\t" << traitB << std::endl;

        perturbedOutput << std::setprecision(20);

        for(size_t i = 0; i < u.size(); ++i) {
          perturbedOutput << u[i] << "\t" << v[i] << std::endl;
        }

        Document perturbedOutputDoc(perturbedOutput, LabelParams(), SeparatorParams('\t'));
        perturbedOutputDoc.Save(perturbedFile);
      }


    omp_set_num_threads(cores);

    std::cout << "Computing the GPS statistic..." << std::endl;

    boost::timer::cpu_timer timer;

    if(statistic == "gps") {
      if(ecdf == "naive") {
        gps = gpsStat(u, v, &bivariateEcdfPar);
      } else if(ecdf == "pp") {
        gps = gpsStat(u, v, &PPEcdf::bivariatePPEcdf);
      } else if(ecdf == "lw") {
        gps = gpsStat(u, v, &bivariateEcdfLW);
      } else {
        std::cout << "Invalid ecdf argument, using naive algorithm" << std::endl;
        gps = gpsStat(u, v, &bivariateEcdfPar);
      }
    } else if(statistic == "mean") {
      if(ecdf == "naive") {
        gps = meanStat(u, v, &bivariateEcdfPar, &gpsWeight);
      } else if(ecdf == "pp") {
        gps = meanStat(u, v, &PPEcdf::bivariatePPEcdf, &gpsWeight);
      } else if(ecdf == "lw") {
        gps = meanStat(u, v, &bivariateEcdfLW, &gpsWeight);
      } else {
        std::cout << "Invalid ecdf argument, using naive algorithm" << std::endl;
        gps = meanStat(u, v, &bivariateEcdfPar, &gpsWeight);
      }
    }

    if(vm.count("timingFile")) {
      timer.stop();
      std::ofstream timingOut(timingFile);
      timingOut << timer.format();
      timingOut.close();
    }

    std::stringstream stringOutput;

    stringOutput << "Trait_A\tTrait_B\tGPS" << std::endl;

    stringOutput << traitA << "\t" << traitB << "\t" << gps << std::endl;

    Document output(stringOutput, LabelParams(), SeparatorParams('\t'));
    output.Save(outputFile);

    if(!logFile.empty()) {
      Document log(logOutput, LabelParams(), SeparatorParams('\t'));
      log.Save(logFile);
    }
  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
