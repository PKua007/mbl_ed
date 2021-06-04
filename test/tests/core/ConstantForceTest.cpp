//
// Created by pkua on 04.06.2021.
//

#include <catch2/catch.hpp>

#include "core/terms/ConstantForce.h"
#include "core/HamiltonianGenerator.h"
#include "core/FockBasisGenerator.h"


TEST_CASE("ConstantForce") {
    SECTION("individual sites") {
        HamiltonianGenerator generator(FockBasisGenerator{}.generate(1, 4), false);
        const auto &fockBase = *generator.getFockBasis();
        ConstantForce constantForce(3);

        REQUIRE(constantForce.calculate(fockBase[0], generator) == Approx(-4.5));
        REQUIRE(constantForce.calculate(fockBase[1], generator) == Approx(-1.5));
        REQUIRE(constantForce.calculate(fockBase[2], generator) == Approx(1.5));
        REQUIRE(constantForce.calculate(fockBase[3], generator) == Approx(4.5));
    }

    SECTION("all sites") {
        HamiltonianGenerator generator(FockBasisGenerator{}.generate(2, 2), false);
        const auto &fockBase = *generator.getFockBasis();
        ConstantForce constantForce(1);

        REQUIRE(constantForce.calculate(fockBase[1], generator) == Approx(0));
    }
}