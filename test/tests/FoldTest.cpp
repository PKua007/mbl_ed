//
// Created by Piotr Kubala on 21/01/2020.
//

#include <catch2/catch.hpp>

#include "utils/Fold.h"

TEST_CASE("Fold: empty") {
    std::string actual = Fold("").width(10);

    REQUIRE(actual == "");
}

TEST_CASE("Fold: normal") {
    std::string actual = Fold("012345 78 01234 6 8 01 23").width(10);

    REQUIRE(actual == "012345 78 \n01234 6 8 \n01 23");
}

TEST_CASE("Fold: touching the border") {
    std::string actual = Fold("0123 5 789 0 234").width(10);

    REQUIRE(actual == "0123 5 789\n0 234");
}

TEST_CASE("Fold: long") {
    std::string actual = Fold("012345 012345678901 34").width(10);

    REQUIRE(actual == "012345 \n0123456789\n01 34");
}

TEST_CASE("Fold: very long") {
    std::string actual = Fold("0123456789012345678901234567890123").width(10);

    REQUIRE(actual == "0123456789\n0123456789\n0123456789\n0123");
}

TEST_CASE("Fold: newlines") {
    std::string actual = Fold("012 56\n01234 01234567\n0 2").width(10);

    REQUIRE(actual == "012 56\n01234 \n01234567\n0 2");
}

TEST_CASE("Fold: multiple spaces") {
    std::string actual = Fold("0123 56     12").width(10);

    REQUIRE(actual == "0123 56   \n 12");
}

TEST_CASE("Fold: space touching border") {
    std::string actual = Fold("0123 56   0123").width(10);

    REQUIRE(actual == "0123 56   \n0123");
}

