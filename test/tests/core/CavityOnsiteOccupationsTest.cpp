//
// Created by Piotr Kubala on 07/12/2020.
//

#include <catch2/catch.hpp>

#include "matchers/VectorApproxEqualCatchMatcher.h"

#include "core/observables/CavityOnsiteOccupations.h"
#include "core/FockBasisGenerator.h"

TEST_CASE("CavityOnsiteOccupations") {
    auto basis = FockBasisGenerator{}.generate(2, 2);
    auto interactions = std::make_unique<CavityLongInteraction>(1, 1./6, 0);
    CavityOnsiteOccupations onsiteOccupations(std::move(basis), std::move(interactions));

    SECTION("header") {
        REQUIRE(onsiteOccupations.getHeader() == std::vector<std::string>{"n_1_cos", "n_2_cos"});
    }

    SECTION("values") {
        onsiteOccupations.calculateForState(arma::cx_vec{1./3, 2./3, 2./3});

        REQUIRE_THAT(onsiteOccupations.getValues(), IsApproxEqual(std::vector<double>{2./3, 4./6}, 1e-8));
    }
}