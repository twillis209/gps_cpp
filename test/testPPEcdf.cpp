#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <PPEcdf.hpp>
#include <iostream>
#include <rapidcsv.h>

using namespace rapidcsv;
using namespace PPEcdf;
using namespace Catch;

TEST_CASE( "idxSort sorts simple data", "[idxSort]" ) {
  std::vector<double> v({0.5, 0.4, 0.3, 0.2, 0.1});
  std::vector<size_t> idx({4, 3, 2, 1, 0});

  std::vector<size_t> sortedIdx = idxSort(v);

  REQUIRE_THAT(
               sortedIdx,
               Matchers::Equals(std::vector<size_t>{4, 3, 2, 1, 0}));
}

TEST_CASE( "idxSort sorts simple duplicated data", "[idxSort]" ) {
  std::vector<double> v({0.4, 0.4, 0.4});
  std::vector<size_t> idx({0, 1, 2});

  std::vector<size_t> sortedIdx = idxSort(v);

  REQUIRE_THAT(
               sortedIdx,
               Matchers::Equals(idx));
}

TEST_CASE( "idxSort sorts simple data with duplicates", "[idxSort]" ) {
  std::vector<double> v({0.4, 0.1, 0.3, 0.2});
  std::vector<size_t> idx({3, 0, 2, 1});

  std::vector<size_t> sortedIdx = idxSort(v);

  REQUIRE_THAT(
               sortedIdx,
               Matchers::Equals(idx));
}

//TEST_CASE( "idxSort sorts simple data with all duplicates but one", "[idxSort]" ) {
//  std::vector<double> v({.4, .4, .4, .4, .1});
//  std::vector<size_t> idx({4, 3, 2, 1, 0});
//
//  std::vector<size_t> sortedIdx = argSort(v);
//
//  REQUIRE_THAT(
//               sortedIdx,
//               Matchers::Equals(std::vector<size_t>{4, 3, 2, 1, 0}));
//}
//
//TEST_CASE( "reindex sorts simple data", "[reindex]" ) {
//  std::vector<double> v({.5, .4, .3, .2, .1});
//  std::vector<size_t> idx({1, 2, 4, 3, 0});
//  std::vector<double> reindexed_v({.1, .5, .4, .2, .3});
//
//  std::vector<double> reindexedV = reindex(v, idx);
//
//  REQUIRE_THAT(
//               reindexedV,
//               Matchers::Approx(std::vector<double>{.1, .5, .4, .2, .3}));
//}

/*
// NB: runif with set.seed(42)
TEST_CASE( "reindex sorts two vectors", "[reindex]" ) {
  //std::vector u({0.9148060, 0.9370754, 0.2861395, 0.8304476, 0.6417455, 0.5190959, 0.7365883, 0.1346666, 0.6569923, 0.7050648});
  std::vector u({0.3, 0.1, 0.2});
  std::vector v({0.4577418, 0.7191123, 0.9346722, 0.2554288, 0.4622928, 0.9400145, 0.9782264, 0.1174874, 0.4749971, 0.5603327});

  //std::vector<size_t> u_sorted_idx({8, 9, 1, 7, 3, 2, 6, 0, 4, 5});
  std::vector<size_t> u_sorted_idx({2, 0, 1});

  //std::vector sorted_u({0.1346666, 0.2861395, 0.5190959, 0.6417455, 0.6569923, 0.7050648, 0.7365883, 0.8304476, 0.9148060, 0.9370754});
  std::vector sorted_u({0.1346666, 0.2861395, 0.5190959, 0.6417455, 0.6569923, 0.7050648, 0.7365883, 0.8304476, 0.9148060, 0.9370754});
  std::vector sorted_v({0.1174874, 0.9346722, 0.9400145, 0.4622928, 0.4749971, 0.5603327, 0.9782264, 0.2554288, 0.4577418, 0.7191123});

  std::vector<size_t> idx = argSort(u);

  for(auto x: idx) std::cout << x << std::endl;

  REQUIRE_THAT(
               u_sorted_idx,
               Matchers::Equals(idx));

  /*
  std::vector<double> u_sorted = reindex(u, idx);

  REQUIRE_THAT(
               sorted_u,
               Matchers::Equals(u_sorted));
  */
//  std::vector<double> reindexed_u = reindex(u, idx);
//  std::vector<double> reindexed_v = reindex(v, idx);
//
//  REQUIRE_THAT(
//               reindexed_u,
//               Matchers::Approx(sorted_u));
//}
/*
TEST_CASE( "bivariatePPEcdf computes the ecdf for simple data", "[bivariatePPEcdf]" ) {
  std::vector<double> v({.5, .4, .3, .2, .1});
  std::vector<double> u({.5, .4, .3, .2, .1});

  std::vector<double> ecdf = bivariatePPEcdf(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(std::vector<double>{1, .8, .6, .4, .2}));
}

TEST_CASE( "bivariatePPEcdf computes the ecdf for less simple data", "[bivariatePPEcdf]" ) {
  std::vector u({.1, .2, .3, .4, .5});
  std::vector v({.1, .2, .3, .5, .4});

  std::vector<double> ecdf = bivariatePPEcdf(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(std::vector<double>{.2, .4, .6, .8, .8}));
}

TEST_CASE( "bivariatePPEcdf computes the ecdf for less simple data (II)", "[bivariatePPEcdf]" ) {
  std::vector u({.1, .2, .3, .4, .5});
  std::vector v({.5, .4, .3, .2, .1});

  std::vector<double> ecdf = bivariatePPEcdf(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(std::vector<double>{.2, .2, .2, .2, .2}));
}
*/
/*
TEST_CASE( "bivariatePPEcdf computes the ecdf for less simple data (III)", "[bivariatePPEcdf]" ) {
  std::vector u({0.12029992, 0.56309735, 0.08356643, 0.97654664, 0.45778035, 0.46165514, 0.71599403, 0.99431279, 0.20568423, 0.87618803});

  std::vector v({0.1232234, 0.9478716, 0.5595810, 0.2650775, 0.5656846, 0.3824650, 0.7844777, 0.4306922, 0.5496253, 0.6820351});

  std::vector<double> ecdf = bivariatePPEcdf(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(std::vector<double>{0.1,0.6,0.1,0.2,0.4,0.2,0.6,0.4,0.2,0.6}));
}
*/
/*
TEST_CASE( "PP bivariate ecdf computes correct ecdf for uniformly distributed data", "[bivariatePPEcdf]" ) {

  Document data("test/data/1e3_unif.csv");
  Document exemplar("test/data/1e3_unif_ecdf.csv");

  std::vector<double> u = data.GetColumn<double>("u");
  std::vector<double> v = data.GetColumn<double>("v");
  std::vector<double> ecdfExemplar = exemplar.GetColumn<double>("ecdf");

  std::vector<double> ecdf = bivariatePPEcdf(u, v);

  for(size_t i; i < u.size(); i++) {
    std::cout << u[i] << "\t" << v[i] << "\t" << ecdf[i] << std::endl;
  }

  REQUIRE_THAT(
              ecdf,
              Matchers::Approx(ecdfExemplar)
              );
}
*/
/*
TEST_CASE( "bivariatePPEcdf computes the ecdf for simple data with duplicates", "[bivariatePPEcdf]" ) {
  std::vector u({.1, .2, .3, .4, .4});
  std::vector v({.1, .2, .3, .5, .5});

  std::vector<double> ecdf = bivariatePPEcdf(u, v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(std::vector<double>{.2, .4, .6, 1., 1.}));
}
*/

/*

  REQUIRE_THAT(
  ecdf,
  Matchers::Approx(std::vector<double>{.167, .333, .5, .667, .667, .667}));
TEST_CASE( "PP bivariate ecdf computes correct ecdf for uniformly distributed data", "[bivariatePPEcdf]" ) {

  Document data("test/data/1e3_unif.csv");
  Document exemplar("test/data/1e3_unif_ecdf.csv");

  std::vector<double> u = data.GetColumn<double>("u");
  std::vector<double> v = data.GetColumn<double>("v");
  std::vector<double> ecdfExemplar = exemplar.GetColumn<double>("ecdf");

  std::vector<double> ecdf = bivariatePPEcdf(u, v);

  for(auto x : ecdf) {
    std::cout << x << std::endl;
  }
//  REQUIRE_THAT(
//               ecdf,
//               Matchers::Approx(ecdfExemplar)
//               );
}

TEST_CASE( "ecdf estimates ecdf correctly on unsorted simple data set", "[ecdf]" ) {
  std::vector u({.1, .2, .5, .4, .3});

  std::vector<double> ecdfResult = ecdf(u);

  REQUIRE_THAT(
               ecdfResult,
               Matchers::Approx(std::vector<double>{.2, .4, 1., .8, .6})
               );
}


TEST_CASE( "gpsStat runs in simple case", "[gpsStat]" ) {
  std::vector u({.1, .2, .3, .4, .5});
  std::vector v({.1, .2, .3, .5, .4});

  double gpsResult = gpsStat(u,v, &bivariateEcdfLW);

  REQUIRE(gpsResult == Detail::Approx(1.439137));
}


TEST_CASE( "Naive bivariate ecdf runs in simple case", "[ecdf]" ) {
  std::vector u({.1, .2, .3, .4, .5});
  std::vector v({.1, .2, .3, .5, .4});

  std::vector<double> ecdfResult = bivariateEcdfPar(u,v);

  std::cout << ecdfResult[0] << std::endl;

  REQUIRE_THAT(
               ecdfResult,
               Matchers::Approx(std::vector<double>{.2, .4, .6, .8, .8})
               );
}

TEST_CASE( "L&W bivariate ecdf runs in simple case", "[ecdf]" ) {
  std::vector u({.1, .2, .3, .4, .5});
  std::vector v({.1, .2, .3, .5, .4});

  std::vector<double> ecdfResult = bivariateEcdfLW(u,v);

  REQUIRE_THAT(
               ecdfResult,
               Matchers::Approx(std::vector<double>{.2, .4, .6, .8, .8})
               );
}

TEST_CASE( "L&W bivariate ecdf throws exception with differently-sized input vectors", "[ecdf]") {

  REQUIRE_THROWS_MATCHES(
                         bivariateEcdfLW(std::vector<double>{.1, .2}, std::vector<double>{.1, .3, .05}),
                         std::invalid_argument,
                         Message("Size of u and v differs.")
                         );
}

TEST_CASE( "L&W bivariate ecdf computes correct ecdf for uniformly distributed data", "[ecdf]" ) {

  Document data("test/data/1e3_unif.csv");
  Document exemplar("test/data/1e3_unif_ecdf.csv");

  std::vector<double> u = data.GetColumn<double>("u");
  std::vector<double> v = data.GetColumn<double>("v");
  std::vector<double> ecdfExemplar = exemplar.GetColumn<double>("ecdf");

  std::vector<double> ecdf = bivariateEcdfLW(u,v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(ecdfExemplar)
               );

}

TEST_CASE( "Naive bivariate ecdf computes correct ecdf for uniformly distributed data", "[ecdf]" ) {

  Document data("test/data/1e3_unif.csv");
  Document exemplar("test/data/1e3_unif_ecdf.csv");

  std::vector<double> u = data.GetColumn<double>("u");
  std::vector<double> v = data.GetColumn<double>("v");
  std::vector<double> ecdfExemplar = exemplar.GetColumn<double>("ecdf");

  std::vector<double> ecdf = bivariateEcdfPar(u,v);

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(ecdfExemplar)
               );

}


TEST_CASE( "Naive bivariate ecdf throws exception with differently-sized input vectors", "[ecdf]") {
  REQUIRE_THROWS_MATCHES(
                         bivariateEcdfPar(std::vector<double>{.1, .2}, std::vector<double>{.1, .3, .05}),
                         std::invalid_argument,
                         Message("Size of u and v differs.")
                         );
}

TEST_CASE( "gps throws exception with differently-sized input vectors", "[ecdf]") {

  REQUIRE_THROWS_MATCHES(
                         gpsStat(std::vector<double>{.1, .2}, std::vector<double>{.1, .3, .05}, &bivariateEcdfLW),
                         std::invalid_argument,
                         Message("Size of u and v differs.")
                         );
}

TEST_CASE( "gps throws exception on undefined case", "[ecdf]") {

  REQUIRE_THROWS_MATCHES(
                         gpsStat(std::vector<double>{.1, .2, .3}, std::vector<double>{.1, .2, .3}, &bivariateEcdfLW),
                         std::invalid_argument,
                         Message("Indices of largest elements of u and v coincide. GPS statistic is undefined in this case.")
                 );
}

TEST_CASE( "perturbDuplicates distinguishes duplicates in simple case", "[perturbDuplicates]") {

  std::vector<double> u({.1,.1,.1,.2,.3,.4,.5});

  std::vector<double> noDup = perturbDuplicates_addEpsilon(u, 1.0);

  REQUIRE(noDup[0] == .1);
  REQUIRE(noDup[0] != noDup[1]);
  REQUIRE(noDup[0] != noDup[2]);
}
*/
