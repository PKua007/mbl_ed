//
// Created by Piotr Kubala on 28/03/2021.
//

#include <catch2/catch.hpp>

#include "analyzer/BandExtractor.h"
#include "core/FockBasisGenerator.h"

TEST_CASE("BandExtractor: single energy set (epsilon range)") {
    BandExtractor extractor(BandExtractor::EpsilonRange(0.5, 0.4));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // 0.4, 0.5, 0.6
    auto indices = extractor.getBandIndices(Eigensystem({0, 0.1, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0}), logger);

    REQUIRE(indices == std::vector<std::size_t>{2, 3, 4});
}

TEST_CASE("MeanGapRatio: single energy set (cdf range)") {
    BandExtractor extractor(BandExtractor::CDFRange(0.5, 0.375));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    //  0.5, 0.6, 0.7
    auto indices = extractor.getBandIndices(Eigensystem({0, 0.2, 0.5, 0.6, 0.7, 0.9, 0.95, 1.0}), logger);

    REQUIRE(indices == std::vector<std::size_t>{2, 3, 4});
}

TEST_CASE("MeanGapRatio: single energy set (around vector)") {
    auto base = std::shared_ptr<FockBasis>(FockBasisGenerator{}.generate(7, 2));
    BandExtractor extractor(BandExtractor::VectorRange(FockBasis::Vector{5, 2}, 0.3));
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

    auto indices = extractor.getBandIndices(Eigensystem(eigval, eigvec, base), logger);

    REQUIRE(indices == std::vector<std::size_t>{2, 3, 4});
}