//
// Created by pkua on 08.11.2019.
//

#include <catch2/catch.hpp>

#include "simulation/FockBaseGenerator.h"
#include "analyzer/tasks/MeanGapRatio.h"

TEST_CASE("GapRadioCalculator: names") {
    MeanGapRatio ratioCalculator(0.5, 0.1);

    REQUIRE(ratioCalculator.getName() == "mgr");
    REQUIRE(ratioCalculator.getResultHeader() == std::vector<std::string>{"meanGapRatio", "meanGapRatioError"});
}

TEST_CASE("MeanGapRatio: single energy set (number)") {
    MeanGapRatio ratioCalculator(0.5, 0.4);
    std::ostringstream logger;

    // So we are averaging for: 0.4, 0.5, 0.6
    ratioCalculator.analyze(Eigensystem({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.6111", "0.2003"});
}

TEST_CASE("MeanGapRatio: single energy set (around vector)") {
    auto base = std::shared_ptr<FockBase>(FockBaseGenerator{}.generate(7, 2));
    MeanGapRatio ratioCalculator(FockBase::Vector{4, 3}, 0.4);
    std::ostringstream logger;

    // {4, 3} has index 3, so we are averagong for: 0.4, 0.5, 0.6
    ratioCalculator.analyze(Eigensystem({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0}, base), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.6111", "0.2003"});
}

TEST_CASE("MeanGapRatio: normalization") {
    MeanGapRatio ratioCalculator(0.5, 0.4);
    std::ostringstream logger;

    ratioCalculator.analyze(Eigensystem({1, 11, 41, 51, 61, 81, 91, 101}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.6111", "0.2003"});
}

TEST_CASE("MeanGapRatio: calculating mean") {
    MeanGapRatio ratioCalculator(0.5, 0.4);
    std::ostringstream logger;

    ratioCalculator.analyze(Eigensystem({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0}), logger);
    ratioCalculator.analyze(Eigensystem({0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.7667", "0.1453"});
}