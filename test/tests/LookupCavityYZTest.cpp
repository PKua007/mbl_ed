//
// Created by Piotr Kubala on 11/02/2020.
//

#include <catch2/catch.hpp>

#include "matchers/ArmaApproxEqualMatcher.h"

#include "simulation/FockBaseGenerator.h"
#include "simulation/terms/LookupCavityYZ.h"
#include "simulation/HamiltonianGenerator.h"

TEST_CASE("LookupCavityYZ: correct") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(2, 2), false);
    CavityConstants cavityConstants;
    cavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}, {4, 5, 6}}});
    cavityConstants.addRealisation(CavityConstants::Realisation{0.5, {{7, 8, 9}, {10, 11, 12}}});

    SECTION("default realisation = 0") {
        LookupCavityYZ term(2, cavityConstants);

        // vector{1, 1} => vector{0, 2}
        // vector{1, 1} * wanniers{2, 5} = 7
        // vector{0, 2} * wanniers{2, 5} = 10
        // hop y constant = 3, -U1/K = -1
        // -1 * 3 * (7 + 10) = -51
        // + vice versa
        REQUIRE(term.calculate({0, 1, {1, 1}, {0, 2}}, generator) == -51);
        REQUIRE(term.calculate({1, 0, {0, 2}, {1, 1}}, generator) == -51);
    }

    SECTION("explicit realisation in constructor = 1") {
        LookupCavityYZ term(2, cavityConstants, 1);

        // vector{1, 1} => vector{0, 2}
        // vector{1, 1} * wanniers{8, 11} = 19
        // vector{0, 2} * wanniers{8, 11} = 22
        // hop y constant = 9, -U1/K = -1
        // + vice versa
        REQUIRE(term.calculate({0, 1, {1, 1}, {0, 2}}, generator) == -9 * (19 + 22));
        REQUIRE(term.calculate({1, 0, {0, 2}, {1, 1}}, generator) == -9 * (19 + 22));
    }

    SECTION("changing realisation to 1") {
        LookupCavityYZ term(2, cavityConstants, 0);

        term.changeRealisation(1);

        REQUIRE(term.calculate({0, 1, {1, 1}, {0, 2}}, generator) == -9 * (19 + 22));
    }
}

TEST_CASE("LookupCavityYZ: integration with HamiltonianGenerator") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(2, 2), false);
    CavityConstants cavityConstants;
    cavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}, {4, 5, 6}}});
    generator.addHoppingTerm(std::make_unique<LookupCavityYZ>(2, cavityConstants));

    arma::mat result = arma::mat(generator.generate());

    arma::mat expected = {{  0, -33,   0},
                          {-33,   0, -51},
                          {  0, -51,   0}};
    expected *= M_SQRT2;
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("LookupCavityYZ: errors") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(2, 2), false);
    CavityConstants cavityConstants;
    cavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}, {4, 5, 6}}});
    cavityConstants.addRealisation(CavityConstants::Realisation{0.5, {{7, 8, 9}, {10, 11, 12}}});

    SECTION("too big realisation in constructor") {
        REQUIRE_THROWS(LookupCavityYZ(2, cavityConstants, 2));
    }

    SECTION("too big realisation in changeRealisation") {
        LookupCavityYZ term(2, cavityConstants, 0);

        REQUIRE_THROWS(term.changeRealisation(2));
    }

    SECTION("too many sites in vector") {
        CavityConstants smallCavityConstants;
        smallCavityConstants.addRealisation(CavityConstants::Realisation{0.4, {{1, 2, 3}}});
        LookupCavityYZ term(2, smallCavityConstants, 0);

        REQUIRE_THROWS(term.calculate({0, 1, {1, 1}, {0, 2}}, generator));
    }
}