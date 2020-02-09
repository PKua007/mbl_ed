//
// Created by Piotr Kubala on 09/02/2020.
//

#include <catch2/catch.hpp>

#include "simulation/terms/HubbardHop.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/FockBaseGenerator.h"

TEST_CASE("HubbardHop") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(4, 4), true);
    HubbardHop hubbardHop(2);

    SECTION("1 hops") {
        REQUIRE(hubbardHop.calculate({}, {}, 0, 1, generator) == -2);
        REQUIRE(hubbardHop.calculate({}, {}, 1, 2, generator) == -2);
        REQUIRE(hubbardHop.calculate({}, {}, 2, 3, generator) == -2);
        REQUIRE(hubbardHop.calculate({}, {}, 3, 0, generator) == -2);
    }

    SECTION("inverted 1 hops") {
        REQUIRE(hubbardHop.calculate({}, {}, 1, 0, generator) == -2);
        REQUIRE(hubbardHop.calculate({}, {}, 2, 1, generator) == -2);
    }
}