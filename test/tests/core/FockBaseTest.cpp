//
// Created by pkua on 01.11.2019.
//

#include <catch2/catch.hpp>

#include "core/FockBasis.h"

TEST_CASE("FockBase: basic operations") {
    SECTION("empty base") {
        FockBasis base;

        REQUIRE(base.size() == 0);
    }

    SECTION("adding and reading elements") {
        FockBasis base;

        base.add(FockBasis::Vector{1, 2, 3});
        base.add(FockBasis::Vector{4, 5, 6});

        REQUIRE(base.size() == 2);
        REQUIRE(base[0] == FockBasis::Vector{1, 2, 3});
        REQUIRE(base[1] == FockBasis::Vector{4, 5, 6});
    }

    SECTION("adding not equally sized elements should throw") {
        FockBasis base;
        base.add(FockBasis::Vector{1, 2, 3});

        REQUIRE_THROWS(base.add(FockBasis::Vector{4, 5, 6, 7}));
    }

    SECTION("iterators") {
        FockBasis base;

        base.add(FockBasis::Vector{1, 2, 3});
        base.add(FockBasis::Vector{4, 5, 6});

        REQUIRE(base.end() - base.begin() == 2);
        REQUIRE(*(base.begin() + 0) == FockBasis::Vector{1, 2, 3});
        REQUIRE(*(base.begin() + 1) == FockBasis::Vector{4, 5, 6});
    }
}

TEST_CASE("FockBase: searching") {
    SECTION("existing search") {
        FockBasis base;
        base.add(FockBasis::Vector{1, 2, 3});
        base.add(FockBasis::Vector{4, 5, 6});
        base.add(FockBasis::Vector{1, 2, 2});

        auto index0 = base.findIndex(FockBasis::Vector{1, 2, 3});
        auto index1 = base.findIndex(FockBasis::Vector{4, 5, 6});
        auto index2 = base.findIndex(FockBasis::Vector{1, 2, 2});

        REQUIRE(*index0 == 0);
        REQUIRE(*index1 == 1);
        REQUIRE(*index2 == 2);
    }

    SECTION("non-existing search") {
        FockBasis base;
        base.add(FockBasis::Vector{1, 2, 3});
        base.add(FockBasis::Vector{4, 5, 6});
        base.add(FockBasis::Vector{1, 2, 2});

        auto index0 = base.findIndex(FockBasis::Vector{1, 1, 3});
        auto index1 = base.findIndex(FockBasis::Vector{4, 5, 6, 7});
        auto index2 = base.findIndex(FockBasis::Vector{2, 2, 1});
        auto index3 = base.findIndex(FockBasis::Vector{1, 2, 3, 0});

        REQUIRE(index0 == std::nullopt);
        REQUIRE(index1 == std::nullopt);
        REQUIRE(index2 == std::nullopt);
        REQUIRE(index3 == std::nullopt);
    }
}

TEST_CASE("FockBase: number of sites") {
    SECTION("get number of sites") {
        FockBasis base;

        base.add(FockBasis::Vector{1, 2, 3});
        base.add(FockBasis::Vector{4, 5, 6});

        REQUIRE(base.getNumberOfSites() == 3);
    }

    SECTION("if no vectors added it should throw") {
        FockBasis base;

        REQUIRE_THROWS(base.getNumberOfSites());
    }
}