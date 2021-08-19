#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <gps.hpp>
#include <iostream>
#include <csv.hpp>

TEST_CASE( "gps computes correct statistic on uniformly distributed data", "[gps]" ) {

  // Probably better libraries if I want to read whole columns, but I have spent too much time on this as it is
  csv::CSVReader reader("test/data/1e6_unif.csv")

  std::vector<double> ecdfResult = bivariateEcdfLW(u,v);

  REQUIRE_THAT(
               ecdfResult,
               Catch::Matchers::Approx(std::vector<double>{.2, .4, .6, .8, .8})
               );
}

TEST_CASE( "L&W bivariate ecdf runs in simple case", "[ecdf]" ) {
  std::vector u({.1, .2, .3, .4, .5});
  std::vector v({.1, .2, .3, .5, .4});

  std::vector<double> ecdfResult = bivariateEcdfLW(u,v);

  REQUIRE_THAT(
               ecdfResult,
               Catch::Matchers::Approx(std::vector<double>{.2, .4, .6, .8, .8})
               );
}

TEST_CASE( "L&W bivariate ecdf throws exception with differently-sized input vectors", "[ecdf]") {

  REQUIRE_THROWS_MATCHES(
                         bivariateEcdfLW(std::vector<double>{.1, .2}, std::vector<double>{.1, .3, .05}),
                         std::invalid_argument,
                         Catch::Message("Size of u and v differs.")
                         );
}

TEST_CASE( "gps throws exception with differently-sized input vectors", "[ecdf]") {

  REQUIRE_THROWS_MATCHES(
                         gps(std::vector<double>{.1, .2}, std::vector<double>{.1, .3, .05}),
                         std::invalid_argument,
                         Catch::Message("Size of u and v differs.")
                         );
}

TEST_CASE( "gps throws exception on undefined case", "[ecdf]") {

  REQUIRE_THROWS_MATCHES(
                         gps(std::vector<double>{.1, .2, .3}, std::vector<double>{.1, .2, .3}),
                         std::invalid_argument,
                         Catch::Message("Indices of largest elements of u and v coincide. GPS statistic is undefined in this case.")
                 );
}
