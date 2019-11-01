//
// Created by pkua on 01.11.2019.
//

#include <catch2/catch.hpp>

#include "FockBaseGenerator.h"

TEST_CASE("3 bosons in 3 sites") {
    FockBaseGenerator generator;

    auto base = generator.generate(3, 3);

    REQUIRE(base.size() == 10);
    REQUIRE(base[0] == std::vector<int>{3, 0, 0});
    REQUIRE(base[1] == std::vector<int>{2, 1, 0});
    REQUIRE(base[2] == std::vector<int>{2, 0, 1});
    REQUIRE(base[3] == std::vector<int>{1, 2, 0});
    REQUIRE(base[4] == std::vector<int>{1, 1, 1});
    REQUIRE(base[5] == std::vector<int>{1, 0, 2});
    REQUIRE(base[6] == std::vector<int>{0, 3, 0});
    REQUIRE(base[7] == std::vector<int>{0, 2, 1});
    REQUIRE(base[8] == std::vector<int>{0, 1, 2});
    REQUIRE(base[9] == std::vector<int>{0, 0, 3});
}