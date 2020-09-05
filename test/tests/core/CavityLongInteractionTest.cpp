//
// Created by Piotr Kubala on 09/02/2020.
//

#include <catch2/catch.hpp>

#include "core/terms/CavityLongInteraction.h"
#include "core/HamiltonianGenerator.h"
#include "core/FockBaseGenerator.h"

TEST_CASE("CavityLongInteraction: +-1 interactions") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(3, 3), false);
    const auto &fockBase = *generator.getFockBase();
    CavityLongInteraction longInteraction(3, 0.5, 0);

    REQUIRE(longInteraction.calculate(fockBase[0], generator) == Approx(-9));
    REQUIRE(longInteraction.calculate(fockBase[1], generator) == Approx(-1));
    REQUIRE(longInteraction.calculate(fockBase[2], generator) == Approx(-9));
    REQUIRE(longInteraction.calculate(fockBase[3], generator) == Approx(-1));
    REQUIRE(longInteraction.calculate(fockBase[4], generator) == Approx(-1));
    REQUIRE(longInteraction.calculate(fockBase[5], generator) == Approx(-9));
    REQUIRE(longInteraction.calculate(fockBase[6], generator) == Approx(-9));
    REQUIRE(longInteraction.calculate(fockBase[7], generator) == Approx(-1));
    REQUIRE(longInteraction.calculate(fockBase[8], generator) == Approx(-1));
    REQUIRE(longInteraction.calculate(fockBase[9], generator) == Approx(-9));
}

TEST_CASE("CavityLongInteraction: beta=1/12, phi0=pi/6") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(2, 3), false);
    const auto &fockBase = *generator.getFockBase();
    CavityLongInteraction longInteraction(3, 1./12, 800);
    longInteraction.setPhi0(M_PI/6);    // Check setPhi0 as well

    REQUIRE(longInteraction.calculate(fockBase[0], generator) == Approx(-3));
    REQUIRE(longInteraction.calculate(fockBase[1], generator) == Approx(-1-std::sqrt(3)/2));
    REQUIRE(longInteraction.calculate(fockBase[2], generator) == Approx(-3./4));
    REQUIRE(longInteraction.calculate(fockBase[3], generator) == Approx(-1));
    REQUIRE(longInteraction.calculate(fockBase[4], generator) == Approx(-1./4));
    REQUIRE(longInteraction.calculate(fockBase[5], generator) == Approx(0).margin(1e-8));
}