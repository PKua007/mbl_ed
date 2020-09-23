//
// Created by Piotr Kubala on 21/02/2020.
//

#include <sstream>

#include <catch2/catch.hpp>

#include "evolution/CorrelationsTimeEntry.h"

TEST_CASE("CorrelationsTimeEntry: no entries") {
    SECTION("no fields") {
        CorrelationsTimeEntry correlationsTimeEntry{};

        REQUIRE(correlationsTimeEntry.toString() == "0 ");
    }

    SECTION("some fields") {
        CorrelationsTimeEntry correlationsTimeEntry(2, 3);

        REQUIRE_THAT(correlationsTimeEntry.toString(), Catch::StartsWith("2 0 0 0 "));
    }
}

TEST_CASE("CorrelationsTimeEntry: single observables set") {
    CorrelationsTimeEntry correlationsTimeEntry(2, 3);

    correlationsTimeEntry.addValues({3, 4, 5});

    REQUIRE(correlationsTimeEntry.toString() == "2 3 4 5 ");
}

TEST_CASE("CorrelationsTimeEntry: averaging") {
    CorrelationsTimeEntry correlationsTimeEntry(2, 3);

    correlationsTimeEntry.addValues({3, 4, 5});
    correlationsTimeEntry.addValues({5, 6, 7});

    REQUIRE(correlationsTimeEntry.toString() == "2 4 5 6 ");
}