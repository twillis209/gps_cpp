#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <gps.hpp>
#include <iostream>
#include <rapidcsv.h>

using namespace rapidcsv;
using namespace gps;
using namespace Catch;
/*
TEST_CASE( "Usual error is thrown with duplicates on small instance", "[ecdf]" ) {
  std::vector u({.1, .1, .3, .4, .5});
  std::vector v({.1, .2, .3, .4, .5});

  std::vector<double> ecdfResult = bivariateEcdfLW(u, v);

//  REQUIRE_THAT(
//               ecdfResult,
//               Matchers::Approx(std::vector<double>{.2, .2, .4, .6, .8, 1.})
//               );

}
*/

TEST_CASE( "Usual error is thrown with duplicates on real data set", "[ecdf]" ) {
  Document data("test/data/pid_still_all_perturbed.tsv", LabelParams(), SeparatorParams('\t'));

  std::vector<double> u = data.GetColumn<double>("pid");
  std::vector<double> v = data.GetColumn<double>("still");

  std::cout.precision(20);

  std::vector<double> ecdf = bivariateEcdfLW(u,v);
}
