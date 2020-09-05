//
// Created by pkua on 01.11.2019.
//

#include <catch2/catch.hpp>

#include "core/FockBaseGenerator.h"

TEST_CASE("FockBaseGenerator: 3 bosons in 3 sites") {
    FockBaseGenerator generator;

    auto base = generator.generate(3, 3);

    REQUIRE(base->size() == 10);
    REQUIRE((*base)[0] == FockBase::Vector{3, 0, 0});
    REQUIRE((*base)[1] == FockBase::Vector{2, 1, 0});
    REQUIRE((*base)[2] == FockBase::Vector{2, 0, 1});
    REQUIRE((*base)[3] == FockBase::Vector{1, 2, 0});
    REQUIRE((*base)[4] == FockBase::Vector{1, 1, 1});
    REQUIRE((*base)[5] == FockBase::Vector{1, 0, 2});
    REQUIRE((*base)[6] == FockBase::Vector{0, 3, 0});
    REQUIRE((*base)[7] == FockBase::Vector{0, 2, 1});
    REQUIRE((*base)[8] == FockBase::Vector{0, 1, 2});
    REQUIRE((*base)[9] == FockBase::Vector{0, 0, 3});
}