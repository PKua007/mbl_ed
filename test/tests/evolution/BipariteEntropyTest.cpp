//
// Created by Piotr Kubala on 23/09/2020.
//

#include <catch2/catch.hpp>

#include "matchers/VectorApproxEqualMatcher.h"

#include "evolution/observables/BipariteEntropy.h"
#include "core/FockBaseGenerator.h"

TEST_CASE("BipariteEntropy") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(2, 4));
    BipariteEntropy bipariteEntropy(fockBase);

    SECTION("header") {
        REQUIRE(bipariteEntropy.getHeader() == std::vector<std::string>{"S"});
    }

    SECTION("values when all subsystems with nonzero entropy") {
        arma::cx_vec vec = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        vec = arma::normalise(vec);

        bipariteEntropy.calculateForState(vec);

        REQUIRE(bipariteEntropy.getValues().size() == 1);
        REQUIRE(bipariteEntropy.getValues().front() == Approx(0.8461643733772487));
    }

    SECTION("values when one of subsystems has zero entropy") {
        arma::cx_vec vec = {1, 2, 0, 0, 5, 0, 0, 8, 9, 10};
        vec = arma::normalise(vec);

        bipariteEntropy.calculateForState(vec);

        REQUIRE(bipariteEntropy.getValues().size() == 1);
        REQUIRE(bipariteEntropy.getValues().front() == Approx(0.3446104320908521));
    }
}