#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <PPEcdf.hpp>
#include <gps.hpp>
#include <iostream>
#include <rapidcsv.h>

using namespace rapidcsv;
using namespace gps;
using namespace Catch;
using namespace PPEcdf;
using namespace std;

TEST_CASE( "idxSort sorts simple data", "[idxSort]" ) {
  vector<double> v({0.5, 0.4, 0.3, 0.2, 0.1});
  vector<size_t> idx({4, 3, 2, 1, 0});

  vector<size_t> sortedIdx = idxSort(v);

  REQUIRE_THAT(
               sortedIdx,
               Matchers::Equals(vector<size_t>{4, 3, 2, 1, 0}));
}

TEST_CASE( "idxSort sorts sorted simple data", "[idxSort]" ) {
  vector<double> v({0.1, 0.2, 0.3, 0.4, 0.5});
  vector<size_t> idx({0, 1, 2, 3, 4});

  vector<size_t> sortedIdx = idxSort(v);

  REQUIRE_THAT(
               sortedIdx,
               Matchers::Equals(idx));
}

TEST_CASE( "idxSort sorts simple duplicated data", "[idxSort]" ) {
  vector<double> v({0.4, 0.4, 0.4});
  vector<size_t> idx({0, 1, 2});

  vector<size_t> sortedIdx = idxSort(v);

  REQUIRE_THAT(
               sortedIdx,
               Matchers::Equals(idx));
}

TEST_CASE( "idxSort sorts triplet of doubles", "[idxSort]" ) {
  vector<double> v({0.2, 0.1, 0.3,});
  vector<size_t> idx({1, 0, 2});

  vector<size_t> sortedIdx = idxSort(v);

  REQUIRE_THAT(
               sortedIdx,
               Matchers::Equals(idx));
}

TEST_CASE( "idxSort sorts simple data with duplicates", "[idxSort]" ) {
  vector<double> v({0.4, 0.1, 0.3, 0.2});
  vector<size_t> idx({1, 3, 2, 0});

  vector<size_t> sortedIdx = idxSort(v);

  vector<double> ridx = reindex(v, sortedIdx);

  REQUIRE_THAT(
               idx,
               Matchers::Equals(sortedIdx));
}

TEST_CASE( "idxSort sorts simple data with all duplicates but one", "[idxSort]" ) {
  vector<double> v({.4, .4, .4, .4, .1});
  vector<size_t> idx({4, 0, 1, 2, 3});

  vector<size_t> sortedIdx = idxSort(v);

  REQUIRE_THAT(
               sortedIdx,
               Matchers::Equals(idx));
}

TEST_CASE( "reindex sorts simple data", "[reindex]" ) {
  vector<double> v({.5, .4, .3, .2, .1});
  vector<size_t> idx({4, 3, 2, 1, 0});
  vector<double> reindexed_v({.1, .2, .3, .4, .5});

  vector<double> reindexedV = reindex(v, idx);

  REQUIRE_THAT(
               reindexedV,
               Matchers::Approx(reindexed_v));
}

// NB: runif with set.seed(42)
TEST_CASE( "reindex sorts two vectors", "[reindex]" ) {
  vector u({0.9148060, 0.9370754, 0.2861395, 0.8304476, 0.6417455, 0.5190959, 0.7365883, 0.1346666, 0.6569923, 0.7050648});
  vector v({0.4577418, 0.7191123, 0.9346722, 0.2554288, 0.4622928, 0.9400145, 0.9782264, 0.1174874, 0.4749971, 0.5603327});

  vector sorted_u({0.1346666, 0.2861395, 0.5190959, 0.6417455, 0.6569923, 0.7050648, 0.7365883, 0.8304476, 0.9148060, 0.9370754});
  vector sorted_v({0.1174874, 0.9346722, 0.9400145, 0.4622928, 0.4749971, 0.5603327, 0.9782264, 0.2554288, 0.4577418, 0.7191123});

  vector<size_t> u_sorted_idx({7, 2, 5, 4, 8, 9, 6, 3, 0, 1});

  vector<size_t> idx = idxSort(u);

  REQUIRE_THAT(
               u_sorted_idx,
               Matchers::Equals(idx));

  vector<double> u_sorted = reindex(u, idx);
  vector<double> v_sorted = reindex(v, idx);

  REQUIRE_THAT(
               sorted_u,
               Matchers::Equals(u_sorted));

  REQUIRE_THAT(
               sorted_v,
               Matchers::Equals(v_sorted));
}


TEST_CASE( "ecdf estimates ecdf correctly on simple data set", "[ecdf]" ) {
  vector u({.1, .2, .3, .4, .5});

  vector<double> ecdfResult = ecdf(u);

  REQUIRE_THAT(
               ecdfResult,
               Matchers::Approx(vector<double>{.2, .4, .6, .8, 1.})
               );
}

TEST_CASE( "ecdf estimates ecdf correctly on unsorted simple data set", "[ecdf]" ) {
  vector u({.1, .2, .5, .4, .3});

  vector<double> ecdfResult = ecdf(u);

  REQUIRE_THAT(
               ecdfResult,
               Matchers::Approx(vector<double>{.2, .4, 1., .8, .6})
               );
}

TEST_CASE( "bivariatePPEcdf computes the ecdf for simple data", "[bivariatePPEcdf]" ) {
  vector<double> v({.5, .4, .3, .2, .1});
  vector<double> u({.5, .4, .3, .2, .1});

  vector<double> ecdf = bivariatePPEcdf(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(vector<double>{1, .8, .6, .4, .2}));
}

TEST_CASE( "bivariatePPEcdfOrbTree computes the ecdf for simple data", "[bivariatePPEcdfOrbTree]" ) {
  std::vector<double> v({.5, .4, .3, .2, .1});
  std::vector<double> u({.5, .4, .3, .2, .1});

  std::vector<double> ecdf = bivariatePPEcdfOrbTree(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(std::vector<double>{1, .8, .6, .4, .2}));
}


TEST_CASE( "bivariatePPEcdf computes the ecdf for less simple data", "[bivariatePPEcdf]" ) {
  vector u({.1, .2, .3, .4, .5});
  vector v({.1, .2, .3, .5, .4});

  vector<double> ecdf = bivariatePPEcdf(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(vector<double>{.2, .4, .6, .8, .8}));
}

TEST_CASE( "bivariatePPEcdfOrbTree computes the ecdf for less simple data", "[bivariatePPEcdfOrbTree]" ) {
  vector u({.1, .2, .3, .4, .5});
  vector v({.1, .2, .3, .5, .4});

  vector<double> ecdf = bivariatePPEcdfOrbTree(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(vector<double>{.2, .4, .6, .8, .8}));
}


TEST_CASE( "bivariatePPEcdf computes the ecdf for less simple data (II)", "[bivariatePPEcdf]" ) {
  vector u({.1, .2, .3, .4, .5});
  vector v({.5, .4, .3, .2, .1});

  vector<double> ecdf = bivariatePPEcdf(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(vector<double>{.2, .2, .2, .2, .2}));
}

TEST_CASE( "bivariatePPEcdfOrbTree computes the ecdf for less simple data (II)", "[bivariatePPEcdfOrbTree]" ) {
  vector u({.1, .2, .3, .4, .5});
  vector v({.5, .4, .3, .2, .1});

  vector<double> ecdf = bivariatePPEcdfOrbTree(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(vector<double>{.2, .2, .2, .2, .2}));
}

TEST_CASE( "bivariatePPEcdf computes the ecdf for less simple data (III)", "[bivariatePPEcdf]" ) {
  vector u({0.9148060, 0.9370754, 0.2861395, 0.8304476, 0.6417455, 0.5190959, 0.7365883, 0.1346666, 0.6569923, 0.7050648});
  vector v({0.4577418, 0.7191123, 0.9346722, 0.2554288, 0.4622928, 0.9400145, 0.9782264, 0.1174874, 0.4749971, 0.5603327});

  vector<double> ecdf = bivariatePPEcdf(u, v);
  vector<double> ecdf_res({0.3, 0.7, 0.2, 0.2, 0.2, 0.3, 0.7, 0.1, 0.3, 0.4});

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(ecdf_res));
}

TEST_CASE( "bivariatePPEcdfOrbTree computes the ecdf for less simple data (III)", "[bivariatePPEcdfOrbTree]" ) {
  vector u({0.9148060, 0.9370754, 0.2861395, 0.8304476, 0.6417455, 0.5190959, 0.7365883, 0.1346666, 0.6569923, 0.7050648});
  vector v({0.4577418, 0.7191123, 0.9346722, 0.2554288, 0.4622928, 0.9400145, 0.9782264, 0.1174874, 0.4749971, 0.5603327});

  vector<double> ecdf = bivariatePPEcdfOrbTree(u, v);
  vector<double> ecdf_res({0.3, 0.7, 0.2, 0.2, 0.2, 0.3, 0.7, 0.1, 0.3, 0.4});

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(ecdf_res));
}


TEST_CASE( "PP bivariate ecdf computes correct ecdf for uniformly distributed data", "[bivariatePPEcdf]" ) {

  Document data("test/data/1e3_unif.csv", LabelParams(), SeparatorParams(','));
  Document exemplar("test/data/1e3_unif_ecdf.csv", LabelParams(), SeparatorParams(','));

  vector<double> u = data.GetColumn<double>("u");
  vector<double> v = data.GetColumn<double>("v");
  vector<double> ecdfExemplar = exemplar.GetColumn<double>("ecdf");

  vector<double> ecdf = bivariatePPEcdf(u, v);

  REQUIRE_THAT(
              ecdf,
              Matchers::Approx(ecdfExemplar)
              );
}

TEST_CASE( "PP bivariate ecdf orbtree computes correct ecdf for uniformly distributed data", "[bivariatePPEcdfOrbTree]" ) {

  Document data("test/data/1e3_unif.csv", LabelParams(), SeparatorParams(','));
  Document exemplar("test/data/1e3_unif_ecdf.csv", LabelParams(), SeparatorParams(','));

  vector<double> u = data.GetColumn<double>("u");
  vector<double> v = data.GetColumn<double>("v");
  vector<double> ecdfExemplar = exemplar.GetColumn<double>("ecdf");

  vector<double> ecdf = bivariatePPEcdfOrbTree(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(ecdfExemplar)
               );
}

TEST_CASE( "bivariatePPEcdfOrbTree computes the ecdf for simple data with duplicate", "[bivariatePPEcdfOrbTree]" ) {
  vector u({.1, .2, .3, .4, .5, .5});
  vector v({.5, .4, .3, .2, .1, .1});

  vector<double> ecdf = bivariatePPEcdfOrbTree(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(vector<double>{0.1666667, 0.1666667, 0.1666667, 0.1666667, 0.3333333, 0.3333333}));
}

TEST_CASE( "gps computes correct statistic on uniformly distributed data", "[gps]" ) {

  Document doc("test/data/1e3_unif.csv");

  vector<double> u = doc.GetColumn<double>("u");
  vector<double> v = doc.GetColumn<double>("v");

  double gpsResult = gpsStat(u,v, &bivariatePPEcdf, &gpsWeight);

  REQUIRE(gpsResult == Detail::Approx(0.87333));
}

TEST_CASE( "gpsStat with bivariatePPEcdfOrbTree computes correct statistic on uniformly distributed data", "[gps]" ) {

  Document doc("test/data/1e3_unif.csv");

  vector<double> u = doc.GetColumn<double>("u");
  vector<double> v = doc.GetColumn<double>("v");

  double gpsResult = gpsStat(u,v, &bivariatePPEcdfOrbTree, &gpsWeight);

  REQUIRE(gpsResult == Detail::Approx(0.87333));
}

TEST_CASE( "gpsStat runs in simple case", "[gpsStat]" ) {
  vector u({.1, .2, .3, .4, .5});
  vector v({.1, .2, .3, .5, .4});

  double gpsResult = gpsStat(u,v, &bivariatePPEcdf, &gpsWeight);

  REQUIRE(gpsResult == Detail::Approx(1.439137));
}


TEST_CASE( "Naive bivariate ecdf runs in simple case", "[ecdf]" ) {
  vector u({.1, .2, .3, .4, .5});
  vector v({.1, .2, .3, .5, .4});

  vector<double> ecdfResult = bivariateEcdfPar(u,v);

  REQUIRE_THAT(
               ecdfResult,
               Matchers::Approx(vector<double>{.2, .4, .6, .8, .8})
               );
}

TEST_CASE( "PP bivariate ecdf throws exception with differently-sized input vectors", "[ecdf]") {

  REQUIRE_THROWS_MATCHES(
                         bivariatePPEcdf(vector<double>{.1, .2}, vector<double>{.1, .3, .05}),
                         invalid_argument,
                         Message("Size of u and v differs.")
                         );
}

TEST_CASE( "Naive bivariate ecdf computes correct ecdf for uniformly distributed data", "[ecdf]" ) {

  Document data("test/data/1e3_unif.csv");
  Document exemplar("test/data/1e3_unif_ecdf.csv");

  vector<double> u = data.GetColumn<double>("u");
  vector<double> v = data.GetColumn<double>("v");
  vector<double> ecdfExemplar = exemplar.GetColumn<double>("ecdf");

  vector<double> ecdf = bivariateEcdfPar(u,v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(ecdfExemplar)
               );

}


TEST_CASE( "Naive bivariate ecdf throws exception with differently-sized input vectors", "[ecdf]") {
  REQUIRE_THROWS_MATCHES(
                         bivariateEcdfPar(vector<double>{.1, .2}, vector<double>{.1, .3, .05}),
                         invalid_argument,
                         Message("Size of u and v differs.")
                         );
}

TEST_CASE( "gps throws exception with differently-sized input vectors", "[ecdf]") {

  REQUIRE_THROWS_MATCHES(
                         gpsStat(vector<double>{.1, .2}, vector<double>{.1, .3, .05}, &bivariatePPEcdf, &gpsWeight),
                         invalid_argument,
                         Message("Size of u and v differs.")
                         );
}

TEST_CASE( "gps throws exception on undefined case", "[ecdf]") {

  REQUIRE_THROWS_MATCHES(
                         gpsStat(vector<double>{.1, .2, .3}, vector<double>{.1, .2, .3}, &bivariatePPEcdf, &gpsWeight),
                         invalid_argument,
                         Message("Indices of largest elements of u and v coincide. GPS statistic is undefined in this case.")
                 );
}

TEST_CASE( "perturbDuplicates distinguishes duplicates in simple case", "[perturbDuplicates]") {

  vector<double> u({.1,.1,.1,.2,.3,.4,.5});

  vector<double> noDup = perturbDuplicates_addEpsilon(u, 1.0);

  REQUIRE(noDup[0] == .1);
  REQUIRE(noDup[0] != noDup[1]);
  REQUIRE(noDup[0] != noDup[2]);
}
