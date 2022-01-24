#include <iostream>
#include <map>
#include <gps.hpp>
#include <rapidcsv.h>
#include <random>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace po;
using namespace gps;
using namespace rapidcsv;

std::vector<double> perturbDuplicates_nextafter(std::vector<double> values) {
  std::map<double, int> freqMap;

  for(size_t i = 0; i < values.size(); ++i){
    freqMap[values[i]]++;

    if(freqMap[values[i]] > 1) {
      for(int j = 1; j < freqMap[values[i]]; ++j) {
        values[i] = nextafter(values[i], DBL_MAX);
      }
    }
  }

  return values;
}

std::vector<double> perturbDuplicates_addEpsilon(std::vector<double> values) {
  std::map<double, int> freqMap;

  for(size_t i = 0; i < values.size(); ++i){
    freqMap[values[i]]++;

    if(freqMap[values[i]] > 1) {
      // TODO found that this works by accident, but it increments by more than epsilon
      for(int j = 1; j < freqMap[values[i]]; ++j) {
        values[i] = values[i] + std::numeric_limits<double>::epsilon();
      }
    }
  }

  return values;
}

std::vector<double> perturbDuplicates_addVariable(std::vector<double> values, double addend) {
  std::map<double, int> freqMap;

  for(size_t i = 0; i < values.size(); ++i){
    freqMap[values[i]]++;

    if(freqMap[values[i]] > 1) {
      // TODO found that this works by accident, but it increments by more than epsilon
      for(int j = 1; j < freqMap[values[i]]; ++j) {
        values[i] += addend;
      }
    }
  }

  return values;
}

std::vector<double> perturbDuplicates_addUnif(std::vector<double> values) {
  std::mt19937 gen(42);
  std::uniform_real_distribution unif(0.0, 1e-10);

  std::map<double, int> freqMap;

  for(size_t i = 0; i < values.size(); ++i){
    freqMap[values[i]]++;

    if(freqMap[values[i]] > 1) {
      // TODO found that this works by accident, but it increments by more than epsilon
      for(int j = 1; j < freqMap[values[i]]; ++j) {
        values[i] += unif(gen);
      }
    }
  }

  return values;
}

int main(int argc, const char* argv[]) {
  po::options_description desc("Allowed options");

  std::string inputFile;
  std::string outputFile;
  std::string colLabel;

  desc.add_options()
    ("help", "Print help message")

    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabel,a", po::value<std::string>(&colLabel), "Label of column")
    ("outputFile,b", po::value<std::string>(&outputFile), "Path to output file")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = data.GetColumn<double>(colLabel);

    std::map<double, int> uFreqMap = returnFreqMap(u);

    std::vector<double> uNoDup = perturbDuplicates_addEpsilon(u);

    std::map<double, int> uNoDupFreqMap = returnFreqMap(uNoDup);

    std::stringstream stringOutput;

    for(size_t i = 0; i < uNoDup.size(); ++i) {
      stringOutput << u[i] << "\t" << uFreqMap[u[i]] << "\t" << uNoDup[i] << "\t" << uNoDupFreqMap[uNoDup[i]] << std::endl;
    }

    Document outputDoc(stringOutput, LabelParams(), SeparatorParams('\t'));
    outputDoc.Save(outputFile);

  } else {
      std::cout << desc << std::endl;
  }

  return 0;
}
