//
// Created by Piotr Kubala on 23/09/2020.
//

#include <catch2/catch.hpp>

#include "matchers/VectorApproxEqualMatcher.h"

#include "evolution/observables/OnsiteOccupations.h"
#include "core/FockBasisGenerator.h"

TEST_CASE("OnsiteOccupations") {
    auto base = std::shared_ptr<FockBasis>(FockBasisGenerator{}.generate(2, 2));
    OnsiteOccupations onsiteOccupations(base);

    SECTION("header") {
        REQUIRE(onsiteOccupations.getHeader() == std::vector<std::string>{"n_1", "n_2"});
    }

    SECTION("values") {
        onsiteOccupations.calculateForState(arma::cx_vec{1./3, 2./3, 2./3});

        REQUIRE_THAT(onsiteOccupations.getValues(), IsApproxEqual(std::vector<double>{2./3, 4./3}, 1e-8));
    }
}