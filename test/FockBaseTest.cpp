//
// Created by pkua on 01.11.2019.
//

#include <catch2/catch.hpp>

#include "FockBase.h"

TEST_CASE("FockBase: basic operations") {
    SECTION("empty base") {
        FockBase base;

        REQUIRE(base.size() == 0);
    }

    SECTION("adding and reading elements") {
        FockBase base;

        base.add(FockBase::Vector{1, 2, 3});
        base.add(FockBase::Vector{4, 5, 6});

        REQUIRE(base.size() == 2);
        REQUIRE(base[0] == FockBase::Vector{1, 2, 3});
        REQUIRE(base[1] == FockBase::Vector{4, 5, 6});
    }

    SECTION("adding not equally sized elements shold throw") {
        FockBase base;
        base.add(FockBase::Vector{1, 2, 3});

        REQUIRE_THROWS(base.add(FockBase::Vector{4, 5, 6, 7}));
    }

    SECTION("iterators") {
        FockBase base;

        base.add(FockBase::Vector{1, 2, 3});
        base.add(FockBase::Vector{4, 5, 6});

        REQUIRE(base.end() - base.begin() == 2);
        REQUIRE(*(base.begin() + 0) == FockBase::Vector{1, 2, 3});
        REQUIRE(*(base.begin() + 1) == FockBase::Vector{4, 5, 6});
    }
}

TEST_CASE("FockBase: searching") {
    SECTION("existing search") {
        FockBase base;
        base.add(FockBase::Vector{1, 2, 3});
        base.add(FockBase::Vector{4, 5, 6});
        base.add(FockBase::Vector{1, 2, 2});

        auto index0 = base.findIndex(FockBase::Vector{1, 2, 3});
        auto index1 = base.findIndex(FockBase::Vector{4, 5, 6});
        auto index2 = base.findIndex(FockBase::Vector{1, 2, 2});

        REQUIRE(*index0 == 0);
        REQUIRE(*index1 == 1);
        REQUIRE(*index2 == 2);
    }

    SECTION("non-existing search") {
        FockBase base;
        base.add(FockBase::Vector{1, 2, 3});
        base.add(FockBase::Vector{4, 5, 6});
        base.add(FockBase::Vector{1, 2, 2});

        auto index0 = base.findIndex(FockBase::Vector{1, 1, 3});
        auto index1 = base.findIndex(FockBase::Vector{4, 5, 6, 7});
        auto index2 = base.findIndex(FockBase::Vector{2, 2, 1});

        REQUIRE(index0 == std::nullopt);
        REQUIRE(index1 == std::nullopt);
        REQUIRE(index2 == std::nullopt);
    }
}