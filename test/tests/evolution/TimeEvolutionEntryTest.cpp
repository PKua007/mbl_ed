//
// Created by Piotr Kubala on 21/02/2020.
//

#include <sstream>

#include <catch2/catch.hpp>

#include "evolution/TimeEvolutionEntry.h"

TEST_CASE("TimeEvolutionEntry: construction") {
    SECTION("default") {
        TimeEvolutionEntry correlationsTimeEntry{};

        REQUIRE(correlationsTimeEntry.toString() == "0 ");
    }

    SECTION("with fields") {
        TimeEvolutionEntry correlationsTimeEntry(2, 3);

        REQUIRE_THAT(correlationsTimeEntry.toString(), Catch::StartsWith("2 0 0 0 "));
    }

    SECTION("from vector") {
        TimeEvolutionEntry correlationsTimeEntry(2, {1, 2, 3});

        REQUIRE_THAT(correlationsTimeEntry.toString(), Catch::StartsWith("2 1 2 3 "));
    }
}

TEST_CASE("TimeEvolutionEntry: single observables set") {
    TimeEvolutionEntry correlationsTimeEntry(2, 3);

    correlationsTimeEntry.addValues({3, 4, 5});

    REQUIRE(correlationsTimeEntry.toString() == "2 3 4 5 ");
}

TEST_CASE("TimeEvolutionEntry: averaging") {
    TimeEvolutionEntry correlationsTimeEntry(2, 3);

    correlationsTimeEntry.addValues({3, 4, 5});
    correlationsTimeEntry.addValues({5, 6, 7});

    REQUIRE(correlationsTimeEntry.toString() == "2 4 5 6 ");
}