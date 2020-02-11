//
// Created by Piotr Kubala on 11/02/2020.
//

#include <catch2/catch.hpp>

#include "simulation/FockBaseGenerator.h"
#include "simulation/terms/LookupCavityZ2.h"
#include "simulation/HamiltonianGenerator.h"

TEST_CASE("LookupCavityZ2: correct") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(2, 2), false);
    CavityConstants cavityConstants;
    cavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}, {4, 5, 6}}});
    cavityConstants.addRealisation(CavityConstants::Realisation{0.5, {{7, 8, 9}, {10, 11, 12}}});

    SECTION("default realisation = 0") {
        LookupCavityZ2 term(2, cavityConstants);

        // (vector{1, 1} dot wanniers{2, 5}) ^ 2 * factor{-U1=2 / K=2 = -1} = -49
        REQUIRE(term.calculate({1, 1}, generator) == -49);
    }

    SECTION("explicitly given realisations index = 1") {
        LookupCavityZ2 term(2, cavityConstants, 1);

        // (vector{1, 1} dot wanniers{8, 11}) ^ 2 * factor{-U1=2 / K=2 = -1} = -361
        REQUIRE(term.calculate({1, 1}, generator) == -361);
    }

    SECTION("changing realisation index") {
        LookupCavityZ2 term(2, cavityConstants, 0);

        term.changeRealisation(1);

        // (vector{1, 1} dot wanniers{8, 11}) ^ 2 * factor{-U1=2 / K=2 = -1} = -361
        REQUIRE(term.calculate({1, 1}, generator) == -361);
    }
}

TEST_CASE("LookupCavityZ2: errors") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(2, 2), false);
    CavityConstants cavityConstants;
    cavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}, {4, 5, 6}}});
    cavityConstants.addRealisation(CavityConstants::Realisation{0.5, {{7, 8, 9}, {10, 11, 12}}});

    SECTION("too big realisation in constructor") {
        REQUIRE_THROWS(LookupCavityZ2(2, cavityConstants, 2));
    }

    SECTION("too big realisation in changeRealisation") {
        LookupCavityZ2 term(2, cavityConstants, 0);

        REQUIRE_THROWS(term.changeRealisation(2));
    }

    SECTION("too many sites in vector") {
        CavityConstants smallCavityConstants;
        smallCavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}}});
        LookupCavityZ2 term(2, smallCavityConstants, 0);

        REQUIRE_THROWS(term.calculate({1, 1}, generator));
    }
}