# README

`gps_cpp` is a C++ library which provides an implementation of the genome-wide pairwise-association signal sharing test (the 'GPS test') invented by Li et al. (DOI: 10.1038/nm.3933).

This project's directory structure was borrowed from a sample project provided by Henry Schreiner [here](https://gitlab.com/CLIUtils/modern-cmake/-/tree/master/examples/extended-project).

## CLI programs

We compute the GPS test statistic for p-values from a pair of GWAS using the `computeGpsCLI` application. We generate null realisations of the GPS test statistic with the `permuteTraitsCLI` application. You can see the use of these programs in the `snakemake` pipeline for our publication 'Accurate detection of shared genetic architecture from GWAS summary statistics in the small-sample context' [here](https://github.com/twillis209/gps_paper_pipeline).

### `computeGpsCLI`

This program expects a tab-separated, uncompressed file with two columns of p-values with labels corresponding to those passed as the `colLabelA` and `colLabelB` command-line arguments.

The essential command-line arguments are as follows:
- `--inputFile/-i`: path to input file
- `--outputFile/-o`: path to output file
- `--colLabelA/-a`: label of first p-value column
- `--colLabelB/-b`: label of second p-value column
- `--traitA/-c`: label of first trait for output file
- `--traitB/-d`: label of second trait for output file

In addition the following are optional arguments:
- `--timingFile/-j`: path to file to contain GPS running time for the purpose of evaluating different ecdf algorithms
- `--logFile/-g`: path to file in which to log input file values which could not be read
- `--ecdf/-f`: string determining ecdf algorithm to use with the GPS (legacy argument)
- `--cores/-n`: number of cores. Only speeds things up when the naive ecdf algorithm is used.

### `permuteTraitsCLI`

This program expects a tab-separated, uncompressed file with two columns of p-values with labels corresponding to those passed as the `colLabelA` and `colLabelB` command-line arguments.

The work of running the permutations is mapped over the cores provided (with the number of cores specified by the `--cores` argument).

Note that the Perisic and Posse ecdf algorithm is used for this program as it is the fastest of the three we considered.

The essential command-line arguments are as follows:
- `--inputFile/-i`: path to input file
- `--outputFile/-o`: path to output file
- `--colLabelA/-a`: label of first p-value column
- `--colLabelB/-b`: label of second p-value column
- `--draws/-n`: the number of permutations to generate

In addition the following are optional arguments:
- `--cores/-n`: number of cores. 

## Computing p-values

The C++ code in this repository suffices to compute the GPS statistic, but not to obtain a p-value for it. We've included a simple CLI R program in the `R` directory which we use to compute a GPS test p-value using a null realisations of the GPS test statistic. The script fits a generalised extreme value distribution (GEVD) to null realisations of the GPS statistic and reports a p-value plus the GEVD parameter estimates and their standard errors.

## Build process

`gps_cpp` is built with CMake. With `gps_cpp` as your working directory:

```
mkdir build
cd build
cmake ..
make
```

### Dependencies

`gps_cpp` depends on the Boost library. CMake will look for this as part of the build process. Using old versions of Boost on my local cluster, it seems `gps_cpp` can be built with a version as old as 1.59.0 (provided it was built with `gcc` 5.4.0, not 4.8.5). It should work with newer versions, too. Note that Boost can be a bit of a pain to install if you're not used to this sort of thing. Using a package manager is probably easiest; `gps_cpp` should be compatible with all versions of the Debian `libboost-all-dev` package listed [here](packages.debian.org/search&keywords=libboost-all-dev).

`gps_cpp` also depends on the `rapidcsv` and `Catch2` libraries, but these should be downloaded and built as part of the build process. See the `CMakeLists.txt` files for more details.

### Unit tests

`Catch2` tests can be run from the top-level `gps_cpp` directory with `./build/test/testGps`; it's necessary to run them from here as several depend on test data files in `gps_cpp/test/data` directory.

## License

This code is licensed under the Lesser GPL. Earlier versions licensed components from the `StOpt` library, also. These have since been removed as we've moved to a different algorithm for computing the bivariate ecdf.
