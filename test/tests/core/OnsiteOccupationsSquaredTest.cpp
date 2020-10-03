//
// Created by Piotr Kubala on 23/09/2020.
//

#include <catch2/catch.hpp>

#include "matchers/VectorApproxEqualCatchMatcher.h"

#include "core/observables/OnsiteOccupationsSquared.h"
#include "core/FockBasisGenerator.h"

TEST_CASE("OnsiteOccupationsSquared") {
    auto base = std::shared_ptr<FockBasis>(FockBasisGenerator{}.generate(2, 2));
    OnsiteOccupationsSquared onsiteOccupations2(base);

    SECTION("header") {
        REQUIRE(onsiteOccupations2.getHeader() == std::vector<std::string>{"n_1N_1", "n_1N_2", "n_2N_2"});
    }

    SECTION("values") {
        onsiteOccupations2.calculateForState(arma::cx_vec{1./3, 2./3, 2./3});

        SymmetricMatrix<double> r(2);
        r(0, 0) = 8./9;
        r(1, 0) = 4./9; r(1, 1) = 20./9;

        REQUIRE_THAT(onsiteOccupations2.getOccupationsSquared(), IsApproxEqual(r, 1e-8));
    }
}