//
// Created by Piotr Kubala on 10/02/2020.
//

#include <catch2/catch.hpp>

#include "core/CavityConstants.h"

TEST_CASE("CavityConstants: basic") {
    CavityConstants cavityConstants;

    REQUIRE(cavityConstants.empty());
    REQUIRE(cavityConstants.size() == 0);
    REQUIRE(cavityConstants.getNumberOfSites() == 0);

    cavityConstants.addRealisation(CavityConstants::Realisation{1, {{1, 2, 3}, {4, 5, 6}}});
    cavityConstants.addRealisation(CavityConstants::Realisation{2, {{7, 8, 9}, {10, 11, 12}}});
    REQUIRE(cavityConstants.size() == 2);
    REQUIRE(!cavityConstants.empty());
    REQUIRE(cavityConstants[0] == CavityConstants::Realisation{1, {{1, 2, 3}, {4, 5, 6}}});
    REQUIRE(cavityConstants[1] == CavityConstants::Realisation{2, {{7, 8, 9}, {10, 11, 12}}});
    REQUIRE(cavityConstants.getNumberOfSites() == 2);
}

TEST_CASE("CavityConstants: error check") {
    CavityConstants cavityConstants;
    cavityConstants.addRealisation(CavityConstants::Realisation{1, {{2, 3, 4}}});

    REQUIRE_THROWS(cavityConstants.addRealisation(CavityConstants::Realisation{2, {{1, 2, 3}, {4, 5, 6}}}));
}