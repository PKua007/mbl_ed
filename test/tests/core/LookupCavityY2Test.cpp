//
// Created by Piotr Kubala on 15/02/2020.
//

#include <catch2/catch.hpp>

#include "core/FockBasisGenerator.h"
#include "core/terms/LookupCavityY2.h"
#include "core/HamiltonianGenerator.h"

TEST_CASE("LookupCavityY2: correct") {
    HamiltonianGenerator generator(FockBasisGenerator{}.generate(2, 3), false);
    CavityConstants cavityConstants;
    cavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}});
    cavityConstants.addRealisation(CavityConstants::Realisation{0.5, {{10, 11, 12}, {13, 14, 15}, {16, 17, 18}}});

    SECTION("default realisation = 0") {
        LookupCavityY2 term(3, cavityConstants);

        // vector{0, 0, 2} => vector{0, 1, 1} => vector{1, 0, 1}
        // hop y constants = 6, 3, -U1/K = -1
        REQUIRE(term.calculate({2, 1, {0, 0, 2}, {0, 1, 1}}, {1, 0, {0, 1, 1}, {1, 0, 1}}, generator) == -18);

        // vector{0, 2, 0} => vector{1, 1, 0} => vector{2, 0, 0}
        // hop y constants = 3, 3, -U1/K = -1
        REQUIRE(term.calculate({1, 0, {0, 2, 0}, {1, 1, 0}}, {1, 0, {1, 1, 0}, {2, 0, 0}}, generator) == -9);

        // vector{2, 0, 0} => vector{1, 1, 0} => vector{2, 0, 0}
        // hop y constants = 3, 3, -U1/K = -1
        REQUIRE(term.calculate({0, 1, {2, 0, 0}, {1, 1, 0}}, {1, 0, {1, 1, 0}, {2, 0, 0}}, generator) == -9);
    }

    SECTION("explicit realisation in constructor = 1") {
        LookupCavityY2 term(3, cavityConstants, 1);

        // vector{0, 0, 2} => vector{0, 1, 1} => vector{1, 0, 1}
        // hop y constants = 15, 12, -U1/K = -1
        REQUIRE(term.calculate({2, 1, {0, 0, 2}, {0, 1, 1}}, {1, 0, {0, 1, 1}, {1, 0, 1}}, generator) == -180);
    }

    SECTION("changing realisation to 1") {
        LookupCavityY2 term(3, cavityConstants, 0);

        term.changeRealisation(1);

        REQUIRE(term.calculate({2, 1, {0, 0, 2}, {0, 1, 1}}, {1, 0, {0, 1, 1}, {1, 0, 1}}, generator) == -180);
    }
}

TEST_CASE("LookupCavityY2: errors") {
    HamiltonianGenerator generator(FockBasisGenerator{}.generate(2, 2), false);
    CavityConstants cavityConstants;
    cavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}, {4, 5, 6}}});
    cavityConstants.addRealisation(CavityConstants::Realisation{0.5, {{7, 8, 9}, {10, 11, 12}}});

    SECTION("too big realisation in constructor") {
        REQUIRE_THROWS(LookupCavityY2(2, cavityConstants, 2));
    }

    SECTION("too big realisation in changeRealisation") {
        LookupCavityY2 term(2, cavityConstants, 0);

        REQUIRE_THROWS(term.changeRealisation(2));
    }

    SECTION("too many sites in vector") {
        CavityConstants smallCavityConstants;
        smallCavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}}});
        LookupCavityY2 term(2, smallCavityConstants, 0);

        REQUIRE_THROWS(term.calculate({0, 1, {2, 0}, {1, 1}}, {0, 1, {1, 1}, {0, 2}}, generator));
    }
}