//
// Created by pkua on 08.11.2019.
//

#include <catch2/catch.hpp>

#include "analyzer/tasks/GapRatioCalculator.h"

TEST_CASE("GapRadioCalculator: names") {
    GapRatioCalculator ratioCalculator(0.5, 0.1);

    REQUIRE(ratioCalculator.getName() == "mean_gap_ratio");
    REQUIRE(ratioCalculator.getResultHeader() == std::vector<std::string>{"mean gap ratio", "mean gap ratio error"});
}

TEST_CASE("GapRatioCalculator: single energy set") {
    GapRatioCalculator ratioCalculator(0.5, 0.4);

    ratioCalculator.analyze({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0});

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.6111", "0.2003"});
}

TEST_CASE("GapRatioCalculator: normalization") {
    GapRatioCalculator ratioCalculator(0.5, 0.4);

    ratioCalculator.analyze({1, 11, 41, 51, 61, 81, 91, 101});

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.6111", "0.2003"});
}

TEST_CASE("GapRatioCalculator: calculating mean") {
    GapRatioCalculator ratioCalculator(0.5, 0.4);

    ratioCalculator.analyze({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0});
    ratioCalculator.analyze({0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0});

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.7667", "0.1453"});
}