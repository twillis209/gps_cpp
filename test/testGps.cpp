#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <gps.hpp>
#include <iostream>
#include <rapidcsv.h>

using namespace rapidcsv;
using namespace gps;
using namespace Catch;

TEST_CASE( "ecdf estimates ecdf correctly on simple data set", "[ecdf]" ) {
  std::vector u({.1, .2, .3, .4, .5});

  std::vector<double> ecdfResult = ecdf(u);

  REQUIRE_THAT(
               ecdfResult,
               Matchers::Approx(std::vector<double>{.2, .4, .6, .8, 1.})
               );
}

TEST_CASE( "ecdf estimates ecdf correctly on unsorted simple data set", "[ecdf]" ) {
  std::vector u({.1, .2, .5, .4, .3});

  std::vector<double> ecdfResult = ecdf(u);

  REQUIRE_THAT(
               ecdfResult,
               Matchers::Approx(std::vector<double>{.2, .4, 1., .8, .6})
               );
}

TEST_CASE( "gps computes correct statistic on uniformly distributed data", "[gps]" ) {

  Document doc("test/data/1e3_unif.csv");

  std::vector<double> u = doc.GetColumn<double>("u");
  std::vector<double> v = doc.GetColumn<double>("v");

  double gpsResult = gpsStat(u,v, &bivariateEcdfLW);

  REQUIRE(gpsResult == Detail::Approx(0.87333));
}

TEST_CASE( "gpsStat runs in simple case", "[gpsStat]" ) {
  std::vector u({.1, .2, .3, .4, .5});
  std::vector v({.1, .2, .3, .5, .4});

  double gpsResult = gpsStat(u,v, &bivariateEcdfLW);

  REQUIRE(gpsResult == Detail::Approx(1.439137));
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
/*
TEST_CASE( "Naive bivariate ecdf computes correct ecdf for uniformly distributed data", "[ecdf]" ) {

  Document data("test/data/1e3_unif.csv");
  Document exemplar("test/data/1e3_unif_ecdf.csv");

  std::vector<double> u = data.GetColumn<double>("u");
  std::vector<double> v = data.GetColumn<double>("v");
  std::vector<double> ecdfExemplar = exemplar.GetColumn<double>("ecdf");

  std::vector<double> ecdf = bivariateEcdfPar(u,v);

  std::cout << ecdf[0] << std::endl;

  REQUIRE_THAT(
               ecdf,
               Matchers::Approx(ecdfExemplar)
               );

}
*/

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

  std::vector<double> noDup = perturbDuplicates(u);

  REQUIRE(noDup[0] == .1);
  REQUIRE(noDup[0] != noDup[1]);
  REQUIRE(noDup[0] != noDup[2]);
}
