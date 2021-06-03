//
// Created by Piotr Kubala on 09/02/2020.
//

#include <catch2/catch.hpp>

#include "core/terms/HubbardHop.h"
#include "core/HamiltonianGenerator.h"
#include "core/FockBasisGenerator.h"

TEST_CASE("HubbardHop: default single hop") {
    HamiltonianGenerator generator(FockBasisGenerator{}.generate(4, 4), true);
    HubbardHop hubbardHop(2);

    SECTION("1 hops") {
        CHECK(hubbardHop.calculate({0, 1, {}, {}}, generator) == -2);
        CHECK(hubbardHop.calculate({1, 2, {}, {}}, generator) == -2);
        CHECK(hubbardHop.calculate({2, 3, {}, {}}, generator) == -2);
        CHECK(hubbardHop.calculate({3, 0, {}, {}}, generator) == -2);
    }

    SECTION("inverted 1 hops") {
        CHECK(hubbardHop.calculate({1, 0, {}, {}}, generator) == -2);
        CHECK(hubbardHop.calculate({2, 1, {}, {}}, generator) == -2);
    }

    SECTION("incorrect hops") {
        CHECK_THROWS(hubbardHop.calculate({0, 0, {}, {}}, generator));
        CHECK_THROWS(hubbardHop.calculate({0, 2, {}, {}}, generator));
    }
}

TEST_CASE("HubbardHop: default two hops") {
    HamiltonianGenerator generator(FockBasisGenerator{}.generate(4, 4), true);
    HubbardHop hubbardHop({2, 3});  // J1 = 2, J2 = 3

    SECTION("correct hops") {
        CHECK(hubbardHop.calculate({0, 1, {}, {}}, generator) == -2);
        CHECK(hubbardHop.calculate({1, 2, {}, {}}, generator) == -2);
        CHECK(hubbardHop.calculate({1, 3, {}, {}}, generator) == -3);
        CHECK(hubbardHop.calculate({1, 0, {}, {}}, generator) == -2);
    }

    SECTION("incorrect hops") {
        CHECK_THROWS(hubbardHop.calculate({0, 0, {}, {}}, generator));
    }
}

TEST_CASE("HubbardHop: only J2") {
    HamiltonianGenerator generator(FockBasisGenerator{}.generate(4, 4), true);
    HubbardHop hubbardHop({2}, {3});  // J1 = 2, J2 = 3

    SECTION("correct hops") {
        CHECK(hubbardHop.calculate({1, 3, {}, {}}, generator) == -3);
        CHECK(hubbardHop.calculate({3, 1, {}, {}}, generator) == -3);
    }

    SECTION("incorrect hops") {
        CHECK_THROWS(hubbardHop.calculate({0, 1, {}, {}}, generator));
        CHECK_THROWS(hubbardHop.calculate({0, 3, {}, {}}, generator));
    }
}