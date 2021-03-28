//
// Created by Piotr Kubala on 06/03/2021.
//

#include <catch2/catch.hpp>

#include "analyzer/tasks/ParticipationEntropy.h"
#include "core/terms/QuasiperiodicDisorder.h"

TEST_CASE("ParticipationEntropy: names") {
    ParticipationEntropy participationEntropy(2, BandExtractor::EpsilonRange(0.5, 0.1));

    REQUIRE(participationEntropy.getName() == "pe");
    REQUIRE(participationEntropy.getResultHeader() == std::vector<std::string>{"participationEntropy",
                                                                               "participationEntropyError"});
}

TEST_CASE("ParticipationEntropy: single energy set") {
    // We are interested in 0.2, 0.4, 0.6
    double sqrt1_5 = std::sqrt(1./5);
    Eigensystem eigensystem(
            {       0, 0.2,        0.4,       0.6,     1.0},

            {{sqrt1_5,   0,          0,   sqrt1_5, sqrt1_5},
             {sqrt1_5,   1,  M_SQRT1_2,         0, sqrt1_5},
             {sqrt1_5,   0,          0,         0, sqrt1_5},
             {sqrt1_5,   0,          0, 2*sqrt1_5, sqrt1_5},
             {sqrt1_5,   0, -M_SQRT1_2,         0, sqrt1_5}});
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    SECTION("q = 1") {
        ParticipationEntropy participationEntropy(1, BandExtractor::EpsilonRange(0.4, 0.7));

        participationEntropy.analyze(eigensystem, logger);

        // Mathematica:
        // S_2(epsilon = 0.2) = 0
        // S_2(epsilon = 0.4) = Log[2]
        // S_2(epsilon = 0.6) = 4/5 Log[5/4] + Log[5]/5
        REQUIRE(participationEntropy.getResultFields() == std::vector<std::string>{"0.39785", "0"});
    }

    SECTION("q = 2") {
        ParticipationEntropy participationEntropy(2, BandExtractor::EpsilonRange(0.4, 0.7));

        participationEntropy.analyze(eigensystem, logger);

        // Mathematica:
        // S_2(epsilon = 0.2) = 0
        // S_2(epsilon = 0.4) = Log[2]
        // S_2(epsilon = 0.6) = Log[25/17]
        REQUIRE(participationEntropy.getResultFields() == std::vector<std::string>{"0.359603", "0"});
    }
}

TEST_CASE("ParticipationEntropy: calculating mean") {
    ParticipationEntropy participationEntropy(2, BandExtractor::EpsilonRange(0.5, 0.3));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // Here: energy 0.4 with pe = 0
    participationEntropy.analyze(Eigensystem({0, 0.4, 0.9, 1}, arma::eye(4, 4)), logger);
    // Here: energies 0.5, 0.6 with avg pe = Log[2]/2
    participationEntropy.analyze(Eigensystem({         0,        0.5, 0.6, 1},
                                             {{M_SQRT1_2, -M_SQRT1_2,   0, 0},
                                              {M_SQRT1_2,  M_SQRT1_2,   0, 0},
                                              {        0,          0,   1, 0},
                                              {        0,          0,   0, 1}}), logger);

    REQUIRE(participationEntropy.getResultFields() == std::vector<std::string>{"0.1733", "0.1733"});
}

TEST_CASE("ParticipationEntropy: requires eigenvectors") {
    std::ostringstream loggerStream;
    Logger logger(loggerStream);
    auto pe = ParticipationEntropy(2, BandExtractor::EpsilonRange(0.5, 0.1));
    REQUIRE_THROWS_WITH(pe.analyze(Eigensystem({0, 1, 2, 3}), logger), Catch::Contains("hasEigenvectors"));
}