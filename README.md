# README

`gps_cpp` is a C++ library which provides an implementation of the genome-wide pairwise-association signal sharing test (the 'GPS test') invented by Li et al. (DOI: 10.1038/nm.3933).

My implementation of the test depends on the very fast multivariate ecdf algorithm published by Langrene and Warin (arXiv:2005.03246) and their implementation of it in the [`StOpt` library](https://gitlab.com/stochastic-control/StOpt). The permutation testing approach we would take with the GPS test would be completely infeasible without this algorithm. I've incorporated a component of the `StOpt` library in the `StOptCDF` directory; this code is licensed under the GNU Lesser General Public License.

This project's directory structure was borrowed from a sample project provided by Henry Schreiner [here](https://gitlab.com/CLIUtils/modern-cmake/-/tree/master/examples/extended-project).

## CLI programs

We compute the GPS test statistic for p-values from a pair of GWAS using the `computeGpsCLI` application. We generate null realisations of the GPS test statistic with the `permuteTraitsCLI` application. You can see the use of these programs in the `snakemake` pipeline for our publication 'Accurate detection of shared genetic architecture from GWAS summary statistics in the small-sample context' at https://github.com/twillis209/gps_paper_pipeline

## The need to perturb the p-values

The divide-and-conquer algorithm of Langrene and Warin depends on the data points (in our use case, p-values) being distinct. However, GWAS summary statistics are usually sufficiently numerous and/or imprecisely reported that there are duplicate p-values amongst them. In order to make use of the fast ecdf algorithm, we 'perturb' duplicate p-values in order to create a set of unique data points. This does not appear to affect the value of the GPS test statistic expressed to four or five significant figures. 

See the `perturbDuplicates_addEpsilon` function for the implementation of this deduplication approach. In short, we 'spread out' duplicate values by adding multiples of epsilon to all but the first value. To fully remove duplicates, we have to iterate this procedure (e.g. see the `perturbN` argument to the `computeGpsCLI` and `permuteTraitsCLI` scripts which specifies the number of iterations). 

We've yet to implement this perturbation approach in a fashion which succeeds with all data sets, but we're working on it. In the meantime, it helps to supply p-values to the highest precision possible; we do this by recalculating p-values from the effect estimates and standard errors commonly supplied in GWAS summary statistics.

In the event the procedure fails to properly deduplicate the values, you'll most likely encounter the following error:

```
computeGpsCLI: opt/include/eigen3/Eigen/src/Core/DenseCoeffsBase.h:365: Eigen::DenseCoeffsBase<Derived, 1>::Scalar& Eigen::DenseCoeffsBase<Derived, 1>::operator()(Eigen::Index, Eigen::Index) [with Derived = Eigen::Array<int, -1, -1>; Eigen::DenseCoeffsBase<Derived, 1>::Scalar = int; Eigen::Index = long int]: Assertion `row >= 0 && row < rows() && col >= 0 && col < cols()' failed.
Aborted
```

We've found it's not sufficient to supply a data set with wholly distinct values in each dimension: the values need to differ by a certain quantum large enough for the divide-and-conquer algorithm to split them into sets.

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
