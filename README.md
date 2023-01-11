# README

`gps_cpp` is a C++ library which provides an implementation of the genome-wide pairwise-association signal sharing test (the 'GPS test') invented by Li et al. (DOI: 10.1038/nm.3933).

There are two implementations of the GPS test in this repo. The first depends on the very fast multivariate ecdf algorithm published by Langrene and Warin (arXiv:2005.03246) and their implementation of it in the [`StOpt` library](https://gitlab.com/stochastic-control/StOpt). The algorithm is O(n log n) in the bivariate case. I've incorporated a component of the `StOpt` library in the `StOptCDF` directory; this code is licensed under the GNU Lesser General Public License. We used the Langrene and Warin ecdf algorithm in the code supporting our paper on the GPS test. This implemention of the test can still be found on the `gps_paper` branch.

My more recent implementation of the GPS test depends on the bivariate ecdf algorithm of Perisic and Posse from 2005 (https://doi.org/10.1198/106186005X69440). This algorithm also runs in O(n log n) time, but I've found it to be significantly faster than Langrene and Warin's (this is perhaps understandable given the latter treats the general multivariate case and the former only the bivariate case). The version of the GPS test on the `master` branch currently uses the Perisic and Posse algorithm. We implemented it in the `PPEcdf` library.

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
- `--perturbedFile/-t`: path to file in which to log the perturbed input data
- `--ecdf/-f`: string determining ecdf algorithm to use with the GPS. Defaults to the naive algorithm which does not require perturbation.
- `--perturbN/-p`: determines the number of iterations of the perturbation procedure. Defaults to 0.
- `--epsilonMultiple/-e`: determines the multiple of `std::numeric_limits<>::epsilon` to add to values in the perturation procedure
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
- `--perturbN/-p`: determines the number of iterations of the perturbation procedure. Defaults to 0.
- `--epsilonMultiple/-e`: determines the multiple of `std::numeric_limits<>::epsilon` to add to values in the perturation procedure
- `--cores/-n`: number of cores. 

## The need to perturb the p-values

The divide-and-conquer algorithm of Langrene and Warin depends on the data points (in our use case, p-values) being distinct. However, GWAS summary statistics are usually sufficiently numerous and/or imprecisely reported that there are duplicate p-values amongst them. In order to make use of the fast ecdf algorithm, we 'perturb' duplicate p-values in order to create a set of unique data points. I've found it necessary to do the same with the Perisic and Posse algorithm for reasons that aren't yet clear; at the moment I put it down to some deficiency of my implementation.

The perturbation procedure does not appear to affect the value of the GPS test statistic expressed to four or five significant figures *in most cases*. You can check the disparity by running `computeGpsCLI` with the command-line argument `-p 0` (no perturbations) and with the argument `-f naive`, which triggers use of the naive bivariate ecdf algorithm. If using this algorithm, you should specify a number of cores with the `-c` flag as it is rather slow. `-f pp` and `-f lw` will trigger use of the Perisic and Posse, and Langrene and Warin ecdf algorithms, respectively.

See the `perturbDuplicates_addEpsilon` function for the implementation of this deduplication approach. In short, we 'spread out' duplicate values by adding multiples of epsilon to all but the first value. To fully remove duplicates, we have to iterate this procedure (e.g. see the `perturbN` argument to the `computeGpsCLI` and `permuteTraitsCLI` scripts which specifies the number of iterations). For the Langrene and Warin implementation, we run 200 iterations of the perturbation procedure to spread out the values sufficiently; this was empirically determined, is particular to the biased selection of GWAS we've looked at so far, and incorporates a generous safety factor. For the Perisic and Posse implementation, one iteration seems to suffice.

We've yet to implement this perturbation approach in a fashion which succeeds with all data sets, but we're working on it. In the meantime, it helps to supply p-values to the highest precision possible; we do this by recalculating p-values from the effect estimates and standard errors commonly supplied in GWAS summary statistics.

In the event the procedure fails to properly deduplicate the values when using the Langrene and Warin algorithm, you'll most likely encounter the following error:

```
computeGpsCLI: opt/include/eigen3/Eigen/src/Core/DenseCoeffsBase.h:365: Eigen::DenseCoeffsBase<Derived, 1>::Scalar& Eigen::DenseCoeffsBase<Derived, 1>::operator()(Eigen::Index, Eigen::Index) [with Derived = Eigen::Array<int, -1, -1>; Eigen::DenseCoeffsBase<Derived, 1>::Scalar = int; Eigen::Index = long int]: Assertion `row >= 0 && row < rows() && col >= 0 && col < cols()' failed.
Aborted
```

We've found it's not sufficient to supply a data set with wholly distinct values in each dimension: the values need to differ by a certain quantum large enough for Langrene and Warin's divide-and-conquer algorithm to split them into sets.

When using the Perisic and Posse algorithm, I compare the GPS statistics I obtain with those computed with the naive algorithm as a sanity check.

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

`gps_cpp` depends on the Boost and Eigen3 libraries. CMake will look for these as part of the build process. In addition, it also depends on the `rapidcsv` and `Catch2` libraries, but these should be downloaded and built as part of the build process. See the `CMakeLists` files for more details.

### Unit tests

`Catch2` tests can be run from the top-level `gps_cpp` directory with `./build/test/testGps`; it's necessary to run them from here as several depend on test data files in `gps_cpp/test/data` directory.

I've incorporated the unit tests from `StOpt` along with the fast ecdf algorithm's implementation taken from there, too. These use the Boost unit testing framework rather than `Catch2`, so at the moment these aren't run together with the `Catch2` tests. Instead they can be run with:

```
./build/test/testFastCDF
./build/test/testFastCDFOnSample
```

## License

This code is licensed under the Lesser GPL. Note that the `StOpt` components we've used have been licensed under the Lesser GPL, also.
