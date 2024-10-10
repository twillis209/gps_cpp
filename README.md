# README

`gps_cpp` is a C++ library which provides an implementation of the genome-wide pairwise-association signal sharing test (the 'GPS test') invented by [Li et al.](https://doi.org/10.1038/nm.3933)

This project's directory structure was borrowed from a sample project provided by Henry Schreiner [here](https://gitlab.com/CLIUtils/modern-cmake/-/tree/master/examples/extended-project).

## CLI programs

We compute the GPS test statistic for p-values from a pair of GWAS using the `computeGpsCLI` application. We generate null realisations of the GPS test statistic with the `permuteTraitsCLI` application. You can see the use of these programs in the `snakemake` pipeline for our [publication](https://doi.org/10.1371/journal.pgen.1010852) '_Accurate detection of shared genetic architecture from GWAS summary statistics in the small-sample context_' [here](https://github.com/twillis209/gps_paper_pipeline) and in use to the end of actually discovering something (rather than merely evaluating the test's performance) in the pipeline for the [paper](https://doi.org/10.1016/j.clim.2024.110356) '_Leveraging pleiotropy identifies common-variant associations with selective IgA deficiency_' [here](https://github.com/twillis209/igad_paper_pipeline).

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
- `--logFile/-g`: path to file in which to log input file values which could not be read

### `permuteTraitsCLI`

This program expects a tab-separated, uncompressed file with two columns of p-values with labels corresponding to those passed as the `colLabelA` and `colLabelB` command-line arguments.

The work of running the permutations is mapped over the cores provided (with the number of cores specified by the `--cores` argument).

Note that the Perisic and Posse ecdf algorithm is used for this program as it is the fastest of the three we considered.

The essential command-line arguments are as follows:
- `--inputFile/-i`: path to input file
- `--outputFile/-o`: path to output file
- `--colLabelA/-a`: label of first p-value column
- `--colLabelB/-b`: label of second p-value column
- `--draws/-d`: the number of permutations to generate

In addition the following are optional arguments:
- `--cores/-n`: number of cores

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

`gps_cpp` depends on Boost, specifically the `multi_index` library. Earlier I required that this be installed on the user's machine prior to the `gps_cpp` build, but I now (October '24) use the `FetchContent` feature of `cmake` to download Boost to get `multi_index` (this is a header-only library). This should save you quite the headache in installing Boost if you're not familiar with this sort of thing (and even if you are!)

`gps_cpp` also depends on the `rapidcsv`, `Catch2`, and `CLI11` libraries, but these should be downloaded and built as part of the build process (again using `FetchContent`). See the `CMakeLists.txt` files for more details.

### Unit tests

`Catch2` tests can be run from the top-level `gps_cpp` directory with `./build/test/testGps`; it's necessary to run them from this working directory as several tests depend on data files in `gps_cpp/test/data` directory.

## Docker

A `docker` image containing the program can be found [here](https://hub.docker.com/r/twillis209/gps-cpp) and can be obtained with

```
docker pull twillis209/gps-cpp:latest
```

## License

This code is licensed under the Lesser GPL. Earlier versions licensed components from the `StOpt` library, also. These have since been removed as we've moved to a different algorithm for computing the bivariate ecdf.
