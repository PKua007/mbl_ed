//
// Created by pkua on 08.11.2019.
//

#include <catch2/catch.hpp>

#include "core/FockBasisGenerator.h"
#include "analyzer/tasks/MeanGapRatio.h"

TEST_CASE("GapRadioCalculator: names") {
    MeanGapRatio ratioCalculator(BandExtractor::EpsilonRange(0.5, 0.1));

    REQUIRE(ratioCalculator.getName() == "mgr");
    REQUIRE(ratioCalculator.getResultHeader() == std::vector<std::string>{"meanGapRatio", "meanGapRatioError"});
}

TEST_CASE("MeanGapRatio: single energy set") {
    MeanGapRatio ratioCalculator(BandExtractor::EpsilonRange(0.5, 0.4));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // So we are averaging for: 0.4, 0.5, 0.6
    ratioCalculator.analyze(Eigensystem({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.611111", "0"});
}

TEST_CASE("MeanGapRatio: normalization") {
    MeanGapRatio ratioCalculator(BandExtractor::EpsilonRange(0.5, 0.4));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    ratioCalculator.analyze(Eigensystem({1, 11, 41, 51, 61, 81, 91, 101}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.611111", "0"});
}

TEST_CASE("MeanGapRatio: calculating mean") {
    MeanGapRatio ratioCalculator(BandExtractor::EpsilonRange(0.5, 0.4));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // Here: 0.4, 0.5, 0.6 - mgr 0.6(1)
    ratioCalculator.analyze(Eigensystem({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0}), logger);
    // Here: 0.4, 0.6 - mgr 0.80(5)
    ratioCalculator.analyze(Eigensystem({0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.8056", "0.1944"});
}