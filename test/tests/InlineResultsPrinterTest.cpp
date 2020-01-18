//
// Created by Piotr Kubala on 18/01/2020.
//

#include <catch2/catch.hpp>

#include "InlineResultsPrinter.h"

TEST_CASE("InlineResultsPrinter: only parameters") {
    Parameters parameters;
    parameters.numberOfBosons = 5;
    parameters.phi0 = "7";

    InlineResultsPrinter irp(parameters, {}, {}, {"numberOfBosons", "phi0"});

    REQUIRE(irp.getHeader() == "numberOfBosons phi0 ");
    REQUIRE(irp.getFields() == "5 7 ");
}

TEST_CASE("InlineResultsPrinter: only inline results") {
    InlineResultsPrinter irp(Parameters{}, {"param1", "param2"}, {"r1", "r2"}, {});

    REQUIRE(irp.getHeader() == "param1 param2 ");
    REQUIRE(irp.getFields() == "r1 r2 ");
}

TEST_CASE("InlineResultsPrinter: both parameters and inline results") {
    Parameters parameters;
    parameters.numberOfBosons = 5;
    InlineResultsPrinter irp(parameters, {"param1"}, {"r1"}, {"numberOfBosons"});

    REQUIRE(irp.getHeader() == "numberOfBosons param1 ");
    REQUIRE(irp.getFields() == "5 r1 ");
}

TEST_CASE("InlineResultsPrinter: unknown parameter should throw") {
    REQUIRE_THROWS(InlineResultsPrinter(Parameters{}, {}, {}, {"killme"}));
}