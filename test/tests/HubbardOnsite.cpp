//
// Created by Piotr Kubala on 09/02/2020.
//

#include <catch2/catch.hpp>

#include "simulation/terms/HubbardOnsite.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/FockBaseGenerator.h"

TEST_CASE("HubbardOnsite") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(4, 2), false);
    const auto &fockBase = generator.getFockBase();
    HubbardOnsite hubbardOnsite(2);

    REQUIRE(hubbardOnsite.calculate(fockBase[0], generator) == 12);
    REQUIRE(hubbardOnsite.calculate(fockBase[1], generator) == 6);
    REQUIRE(hubbardOnsite.calculate(fockBase[2], generator) == 4);
    REQUIRE(hubbardOnsite.calculate(fockBase[3], generator) == 6);
    REQUIRE(hubbardOnsite.calculate(fockBase[4], generator) == 12);
}