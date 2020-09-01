//
// Created by Piotr Kubala on 30/01/2020.
//

#include <catch2/catch.hpp>

#include "analyzer/tasks/InverseParticipationRatio.h"

TEST_CASE("InverseParticipationRatio: names") {
    InverseParticipationRatio ratioCalculator(0.5, 0.1);

    REQUIRE(ratioCalculator.getName() == "ipr");
}

TEST_CASE("InverseParticipationRatio: single energy set") {
    InverseParticipationRatio ratioCalculator(0.5, 0.7);
    std::ostringstream logger;

    // We are interested in 0.4, 0.5, 0.6, 0.8
    Eigensystem eigensystem(
            { 0, 0.1, 0.4,   0.5,   0.6,        0.8, 0.9, 1.0},

            {{1,   1,   0,  2./3,     0,          0,   1,   1},
             {1,   1,   0,  1./3,  1./2,          0,   1,   1},
             {1,   1,   1,  1./3,     0,  M_SQRT1_2,   1,   1},
             {1,   1,   0, -1./3, -1./2,          0,   1,   1},
             {1,   1,   0,     0,  1./2,          0,   1,   1},
             {1,   1,   0,  1./3,  1./2,          0,   1,   1},
             {1,   1,   0,     0,     0,          0,   1,   1},
             {1,   1,   0, -1./3,     0, -M_SQRT1_2,   1,   1}});
    ratioCalculator.analyze(eigensystem, logger);

    // ( (1)^4 )^(-1)  =  1
    // ( (2/3)^4 + 5*(1/3)^4 )^(-1)  =  3.85714285714
    // ( 4*(1/2)^4 )^(-1)  =  4
    // ( 2*(1/sqrt(2))^4 )^(-1)  =  2
    std::ostringstream out;
    ratioCalculator.storeResult(out);
    REQUIRE(out.str() == "0.4 1\n0.5 3.85714\n0.6 4\n0.8 2\n");
}

TEST_CASE("InverseParticipationRatio: multiple simulations") {
    InverseParticipationRatio ratioCalculator(0.5, 0.2);
    std::ostringstream logger;

    ratioCalculator.analyze(Eigensystem(
            { 0, 0.4, 1},

            {{1,   0, 1},
             {1,   1, 1},
             {1,   0, 1}}), logger);
    ratioCalculator.analyze(Eigensystem(
            { 0,        0.5, 1},

            {{1, -M_SQRT1_2, 1},
             {1,  M_SQRT1_2, 1},
             {1,          0, 1}}), logger);

    // ( (1)^4 )^(-1)  =  1
    // ( 2*(1/sqrt(2))^4 )^(-1)  =  2
    std::ostringstream out;
    ratioCalculator.storeResult(out);
    REQUIRE(out.str() == "0.4 1\n0.5 2\n");
}

TEST_CASE("InverseParticipationRatio: requires eigenvectors") {
    std::ostringstream logger;
    REQUIRE_THROWS_WITH(InverseParticipationRatio(0.5, 0.1).analyze(Eigensystem({0, 1, 2, 3}), logger),
                        Catch::Contains("hasEigenvectors"));
}