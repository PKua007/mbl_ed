//
// Created by Piotr Kubala on 16/03/2020.
//

#include <catch2/catch.hpp>

#include "core/terms/QuasiperiodicDisorder.h"
#include "core/FockBaseGenerator.h"
#include "core/HamiltonianGenerator.h"

TEST_CASE("QuasiperiodicDisorder") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(2, 3), false);
    const auto &fockBase = *generator.getFockBase();
    QuasiperiodicDisorder quasiperiodicDisorder(3, 1./12, 800);
    quasiperiodicDisorder.setPhi0(M_PI/6);    // Check setPhi0 as well

    REQUIRE(quasiperiodicDisorder.calculate(fockBase[0], generator) == Approx(3*std::sqrt(3)));
    REQUIRE(quasiperiodicDisorder.calculate(fockBase[1], generator) == Approx(1.5*(std::sqrt(3) + 1)));
    REQUIRE(quasiperiodicDisorder.calculate(fockBase[2], generator) == Approx(1.5*std::sqrt(3)));
    REQUIRE(quasiperiodicDisorder.calculate(fockBase[3], generator) == Approx(3));
    REQUIRE(quasiperiodicDisorder.calculate(fockBase[4], generator) == Approx(1.5));
    REQUIRE(quasiperiodicDisorder.calculate(fockBase[5], generator) == Approx(0).margin(1e-8));
}