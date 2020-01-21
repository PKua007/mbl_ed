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

TEST_CASE("Fold: margin") {
    std::string actual = Fold("012 4 0 23\n01234 01").width(10).margin(5);

    REQUIRE(actual == "     012 4\n     0 23\n     01234\n     01");
}

TEST_CASE("Fold: max width in last line") {
    std::string actual2 = Fold("0 2 4567 0123 5 789").width(10);
    std::string actual1 = Fold("0123456 8 0123 56   ").width(10);

    REQUIRE(actual1 == "0123456 8 \n0123 56   ");
    REQUIRE(actual2 == "0 2 4567 \n0123 5 789");
}