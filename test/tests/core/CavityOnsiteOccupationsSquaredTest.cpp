//
// Created by Piotr Kubala on 07/12/2020.
//

#include <catch2/catch.hpp>

#include "matchers/VectorApproxEqualCatchMatcher.h"

#include "core/observables/CavityOnsiteOccupationsSquared.h"
#include "core/FockBasisGenerator.h"

TEST_CASE("CavityOnsiteOccupationsSquared") {
    auto basis = FockBasisGenerator{}.generate(2, 2);
    auto interactions = std::make_unique<CavityLongInteraction>(1, 1./6, 0);
    CavityOnsiteOccupationsSquared onsiteOccupations2(std::move(basis), std::move(interactions));

    SECTION("header") {
        REQUIRE(onsiteOccupations2.getHeader() == std::vector<std::string>{"n_1N_1_cos", "n_1N_2_cos", "n_2N_2_cos"});
    }

    SECTION("values") {
        onsiteOccupations2.calculateForState(arma::cx_vec{1./3, 2./3, 2./3});

        SymmetricMatrix<double> r(2);
        r(0, 0) = 8./9;
        r(1, 0) = 2./9; r(1, 1) = 5./9;

        REQUIRE_THAT(onsiteOccupations2.getOccupationsSquared(), IsApproxEqual(r, 1e-8));
    }
}