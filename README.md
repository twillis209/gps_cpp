# README

`gps_cpp` is a C++ library which provides an implementation of the genome-wide pairwise-association signal sharing test (the 'GPS test') invented by Li et al. (DOI: 10.1038/nm.3933).

My implementation of the test depends on the very fast multivariate ecdf algorithm published by Langrene and Warin (arXiv:2005.03246) and their implementation of it in the [`StOpt` library](https://gitlab.com/stochastic-control/StOpt). The permutation testing approach we would took with the GPS test would be completely infeasible without this algorithm. I've incorporated a component of the `StOpt`library in the `StOptCDF` directory; this code is licensed under the GNU Lesser General Public License.

This project's directory structure was borrowed from a sample project provided by Henry Schreiner [here](https://gitlab.com/CLIUtils/modern-cmake/-/tree/master/examples/extended-project).

## CLI programs

We compute the GPS test statistic for p-values from a pair of GWAS using the `computeGpsCLI` application.

## The need to perturb the p-values

The divide-and-conquer algorithm of Langrene and Warin depends on the data points (in our use case, p-values) being distinct. However, GWAS summary statistics are usually sufficiently numerous and/or imprecisely reported that there are many duplicate values amongst them. In order to make use of the fast ecdf algorithm, we 'perturb' duplicate p-values in order to create a set of unique data points. This does not appear to affect the value of the GPS test statistic expressed to four or five significant figures. See the `perturbDuplicates_addEpsilon` function for the implementation of this deduplication approach.
