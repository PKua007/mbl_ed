//
// Created by pkua on 08.11.2019.
//

#include <catch2/catch.hpp>

#include "core/FockBasisGenerator.h"
#include "analyzer/tasks/MeanGapRatio.h"

TEST_CASE("GapRadioCalculator: names") {
    MeanGapRatio ratioCalculator(EigenvaluesExtractor::EpsilonRange(0.5, 0.1));

    REQUIRE(ratioCalculator.getName() == "mgr");
    REQUIRE(ratioCalculator.getResultHeader() == std::vector<std::string>{"meanGapRatio", "meanGapRatioError"});
}

TEST_CASE("MeanGapRatio: single energy set (epsilon range)") {
    MeanGapRatio ratioCalculator(EigenvaluesExtractor::EpsilonRange(0.5, 0.4));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // So we are averaging for: 0.4, 0.5, 0.6
    ratioCalculator.analyze(Eigensystem({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.611111", "0"});
}

TEST_CASE("MeanGapRatio: single energy set (cdf range)") {
    MeanGapRatio ratioCalculator(EigenvaluesExtractor::CDFRange(0.5, 0.375));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // So we are averaging for: 0.5, 0.6, 0.7
    ratioCalculator.analyze(Eigensystem({0, 0.2, 0.5, 0.6, 0.7, 0.9, 0.95, 1.0}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.611111", "0"});
}

TEST_CASE("MeanGapRatio: single energy set (around vector)") {
    auto base = std::shared_ptr<FockBasis>(FockBasisGenerator{}.generate(7, 2));
    MeanGapRatio ratioCalculator(EigenvaluesExtractor::VectorRange(FockBasis::Vector{5, 2}, 0.3));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // Here, we want our middle vector |52> to correspond to eigenenergy 0.5, around which we will calculate mgr,
    // so eigenbasis permutes |52>, |43> and |34>
    // Margin 0.3 means that 0.4, 0.5, 0.6 are taken into account
    arma::vec eigval = {0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0};
    arma::mat eigvec = {{1, 0, 0, 0, 0, 0, 0, 0},
                        {0, 1, 0, 0, 0, 0, 0, 0},
                        {0, 0, 0, 1, 0, 0, 0, 0},
                        {0, 0, 0, 0, 1, 0, 0, 0},
                        {0, 0, 1, 0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0, 1, 0, 0},
                        {0, 0, 0, 0, 0, 0, 1, 0},
                        {0, 0, 0, 0, 0, 0, 0, 1}};

    ratioCalculator.analyze(Eigensystem(eigval, eigvec, base), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.611111", "0"});
}

TEST_CASE("MeanGapRatio: normalization") {
    MeanGapRatio ratioCalculator(EigenvaluesExtractor::EpsilonRange(0.5, 0.4));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    ratioCalculator.analyze(Eigensystem({1, 11, 41, 51, 61, 81, 91, 101}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.611111", "0"});
}

TEST_CASE("MeanGapRatio: calculating mean") {
    MeanGapRatio ratioCalculator(EigenvaluesExtractor::EpsilonRange(0.5, 0.4));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // Here: 0.4, 0.5, 0.6 - mgr 0.6(1)
    ratioCalculator.analyze(Eigensystem({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0}), logger);
    // Here: 0.4, 0.6 - mgr 0.80(5)
    ratioCalculator.analyze(Eigensystem({0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"0.8056", "0.1944"});
}