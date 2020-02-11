//
// Created by Piotr Kubala on 10/02/2020.
//

#include <sstream>

#include <catch2/catch.hpp>

#include "frontend/CavityConstantsReader.h"

TEST_CASE("CavityConstantsReader: empty input") {
    std::istringstream in;

    CavityConstants constants = CavityConstantsReader::load(in);

    REQUIRE(constants.empty());
}

TEST_CASE("CavityConstantsReader: correct input") {
    std::stringstream in;
    in << "1 11 12 13" << std::endl;
    in << "1 21 22 23" << std::endl;
    in << "1 31 32 33" << std::endl;
    in << "2 11 12 13" << std::endl;
    in << "2 21 22 23" << std::endl;
    in << "2 31 32 33" << std::endl;

    CavityConstants constants = CavityConstantsReader::load(in);

    REQUIRE(constants.size() == 2);
    REQUIRE(constants[0] == CavityConstants::Realisation{1, {{11, 12, 13}, {21, 22, 23}, {31, 32, 33}}});
    REQUIRE(constants[1] == CavityConstants::Realisation{2, {{11, 12, 13}, {21, 22, 23}, {31, 32, 33}}});
}

TEST_CASE("CavityConstantsReader: last realisation with too few site entries") {
    std::stringstream in;
    in << "1 11 12 13" << std::endl;
    in << "1 21 22 23" << std::endl;
    in << "1 31 32 33" << std::endl;
    in << "2 11 12 13" << std::endl;
    in << "2 21 22 23" << std::endl;

    REQUIRE_THROWS_WITH(CavityConstantsReader::load(in), Catch::Contains("too few"));
}

TEST_CASE("CavityConstantsReader: different number of sites in realisations") {
    std::stringstream in;
    in << "1 11 12 13" << std::endl;
    in << "1 21 22 23" << std::endl;
    in << "1 31 32 33" << std::endl;
    in << "2 11 12 13" << std::endl;
    in << "2 21 22 23" << std::endl;
    in << "3 21 22 23" << std::endl;

    REQUIRE_THROWS_WITH(CavityConstantsReader::load(in), Catch::Contains("Different"));
}