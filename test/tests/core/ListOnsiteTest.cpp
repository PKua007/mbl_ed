//
// Created by Piotr Kubala on 20/03/2020.
//

#include <catch2/catch.hpp>

#include "core/terms/ListOnsite.h"
#include "core/FockBasisGenerator.h"
#include "core/HamiltonianGenerator.h"

TEST_CASE("ListOnsite") {
    HamiltonianGenerator generator(FockBasisGenerator{}.generate(2, 3), false);
    const auto &fockBase = *generator.getFockBasis();
    ListOnsite listOnsite({1, 2, 3});

    REQUIRE(listOnsite.calculate(fockBase[0], generator) == 2);
    REQUIRE(listOnsite.calculate(fockBase[1], generator) == 3);
    REQUIRE(listOnsite.calculate(fockBase[2], generator) == 4);
    REQUIRE(listOnsite.calculate(fockBase[3], generator) == 4);
    REQUIRE(listOnsite.calculate(fockBase[4], generator) == 5);
    REQUIRE(listOnsite.calculate(fockBase[5], generator) == 6);
}