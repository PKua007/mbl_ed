//
// Created by Piotr Kubala on 23/01/2020.
//

#include <catch2/catch.hpp>

#include "analyzer/tasks/MeanInverseParticipationRatio.h"

TEST_CASE("MeanInverseParticipationRatio: names") {
    MeanInverseParticipationRatio ratioCalculator(0.5, 0.1);

    REQUIRE(ratioCalculator.getName() == "ipr");
    REQUIRE(ratioCalculator.getResultHeader() == std::vector<std::string>{"inverseParticipationRatio",
                                                                          "inverseParticipationRatioError"});
}

TEST_CASE("MeanInverseParticipationRatio: single energy set") {
    MeanInverseParticipationRatio ratioCalculator(0.5, 0.7);

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
    ratioCalculator.analyze(eigensystem);

    // (   ( (1)^4 )^(-1)
    //   + ( (2/3)^4 + 5*(1/3)^4 )^(-1)
    //   + ( 4*(1/2)^4 )^(-1)
    //   + ( 2*(1/sqrt(2))^4 )^(-1)     ) / 4  =  2.71428571429
    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"2.7143", "0.7308"});
}

TEST_CASE("MeanInverseParticipationRatio: calculating mean") {
    MeanInverseParticipationRatio ratioCalculator(0.5, 0.2);

    ratioCalculator.analyze(Eigensystem(
            { 0, 0.4, 1},

            {{1,   0, 1},
             {1,   1, 1},
             {1,   0, 1}}));
    ratioCalculator.analyze(Eigensystem(
            { 0,        0.5, 1},

            {{1, -M_SQRT1_2, 1},
             {1,  M_SQRT1_2, 1},
             {1,          0, 1}}));

    // (   ( (1)^4 )^(-1)
    //   + ( 2*(1/sqrt(2))^4 )^(-1) ) / 2  =  1.5
    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"1.5000", "0.5000"});
}

TEST_CASE("MeanInverseParticipationRatio: requires eigenvectors") {
    REQUIRE_THROWS_WITH(MeanInverseParticipationRatio(0.5, 0.1).analyze(Eigensystem({0, 1, 2, 3})),
                        Catch::Contains("hasEigenvectors"));
}