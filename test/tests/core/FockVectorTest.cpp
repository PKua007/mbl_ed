//
// Created by Piotr Kubala on 18/07/2020.
//

#include <catch2/catch.hpp>
#include <sstream>

#include "core/FockVector.h"

TEST_CASE("FockVector: string representation") {
    SECTION("ok") {
        FockVector v("1.0.3.5");

        REQUIRE(v == FockVector{1, 0, 3, 5});
    }

    SECTION("wrong") {
        REQUIRE_THROWS(FockVector("-5.5"));
        REQUIRE_THROWS(FockVector("4.a.6"));
    }
}

TEST_CASE("FockVector: tag representation") {
    SECTION("ok") {
        FockVector unif(5, "unif");
        FockVector dw(6, "dw");
        FockVector empty(0, "Mumbo Jumbo you are afk!");

        REQUIRE(unif == FockVector{1, 1, 1, 1, 1});
        REQUIRE(dw == FockVector{2, 0, 2, 0, 2, 0});
        REQUIRE(empty.empty());
    }

    SECTION("wrong") {
        REQUIRE_THROWS(FockVector(5, "dw"));
        REQUIRE_THROWS(FockVector(4, "non existing tag"));
    }
}

TEST_CASE("FockVector: printing") {
    FockVector empty{};
    FockVector v{0, 3, 7, 0};

    std::ostringstream out;
    out << empty << "|" << v;
    REQUIRE(out.str() == "|0.3.7.0");
}