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
  std::string aOutputFile;
  std::string bOutputFile;
  std::string colLabelA;
  std::string colLabelB;

  desc.add_options()
    ("help", "Print help message")

    ("inputFile,i", po::value<std::string>(&inputFile), "Path to input file")
    ("colLabelA,a", po::value<std::string>(&colLabelA), "Label of column A")
    ("colLabelB,b", po::value<std::string>(&colLabelB), "Label of column B")
    ("aOutputFile,c", po::value<std::string>(&aOutputFile), "Path to A output file")
    ("bOutputFile,d", po::value<std::string>(&bOutputFile), "Path to B output file")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if(vm.count("inputFile")) {
    Document data(inputFile, LabelParams(), SeparatorParams('\t'));

    std::vector<double> u = data.GetColumn<double>(colLabelA);
    std::vector<double> v = data.GetColumn<double>(colLabelB);

    std::map<double, int> uFreqMap = returnFreqMap(u);
    std::map<double, int> vFreqMap = returnFreqMap(v);

    std::vector<double> uNoDup = perturbDuplicates_addUnif(u);
    std::vector<double> vNoDup = perturbDuplicates_addUnif(v);

    std::map<double, int> uNoDupFreqMap = returnFreqMap(uNoDup);
    std::map<double, int> vNoDupFreqMap = returnFreqMap(vNoDup);

    std::stringstream aStringOutput;

    for(size_t i = 0; i < uNoDup.size(); ++i) {
      aStringOutput << u[i] << "\t" << uFreqMap[u[i]] << "\t" << uNoDup[i] << "\t" << uNoDupFreqMap[uNoDup[i]] << std::endl;
    }

    Document aOutput(aStringOutput, LabelParams(), SeparatorParams('\t'));
    aOutput.Save(aOutputFile);

    std::stringstream bStringOutput;

    for(size_t i = 0; i < vNoDup.size(); ++i) {
      bStringOutput << v[i] << "\t" << vFreqMap[v[i]] << "\t" << vNoDup[i] << "\t" << vNoDupFreqMap[vNoDup[i]] << std::endl;
    }

    Document bOutput(bStringOutput, LabelParams(), SeparatorParams('\t'));
    bOutput.Save(bOutputFile);
  } else {
      std::cout << desc << std::endl;

      // TODO remove me
      double testDoubleA = 0.1;
      double testDoubleB = nextafter(testDoubleA, 1.0);
      double testDoubleC = nextafter(testDoubleA, DBL_MAX);

      std::cout.precision(20);

      std::cout << testDoubleA << std::endl;
      std::cout << testDoubleB << std::endl;
      std::cout << testDoubleC << std::endl;
      std::cout << (testDoubleA == testDoubleB) << std::endl;
  }

  return 0;
}
