//
// Created by pkua on 08.11.2019.
//

#include <Catch2/catch.hpp>

#include "GapRatioCalculator.h"

TEST_CASE("GapRatioCalculator: single energy set") {
    std::vector<double> energies = {0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0};

    GapRatioCalculator ratioCalculator(0.5, 0.4);
    ratioCalculator.addEigenenergies(energies);

    REQUIRE(ratioCalculator.calculateMean().value == Approx(11./18));
}

TEST_CASE("GapRatioCalculator: normalization") {
    std::vector<double> energies = {1, 11, 41, 51, 61, 81, 91, 101};

    GapRatioCalculator ratioCalculator(0.5, 0.4);
    ratioCalculator.addEigenenergies(energies);

    REQUIRE(ratioCalculator.calculateMean().value == Approx(11./18));
}

TEST_CASE("GapRatioCalculator: calculating mean") {
    std::vector<double> energies1 = {0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0};
    std::vector<double> energies2 = {0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0};

    GapRatioCalculator ratioCalculator(0.5, 0.4);
    ratioCalculator.addEigenenergies(energies1);
    ratioCalculator.addEigenenergies(energies2);

    REQUIRE(ratioCalculator.calculateMean().value == Approx(23./30));
}