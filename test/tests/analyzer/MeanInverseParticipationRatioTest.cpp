//
// Created by Piotr Kubala on 23/01/2020.
//

#include <catch2/catch.hpp>

#include "analyzer/tasks/MeanInverseParticipationRatio.h"

#include "core/FockBasisGenerator.h"
#include "core/terms/HubbardHop.h"
#include "core/terms/HubbardOnsite.h"
#include "core/terms/QuasiperiodicDisorder.h"
#include "core/averaging_models/UniformPhi0AveragingModel.h"

TEST_CASE("MeanInverseParticipationRatio: names") {
    MeanInverseParticipationRatio ratioCalculator(0.5, 0.1);

    REQUIRE(ratioCalculator.getName() == "mipr");
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
    std::ostringstream loggerStream;
    Logger logger(loggerStream);
    ratioCalculator.analyze(eigensystem, logger);

    // (   ( (1)^4 )^(-1)
    //   + ( (2/3)^4 + 5*(1/3)^4 )^(-1)
    //   + ( 4*(1/2)^4 )^(-1)
    //   + ( 2*(1/sqrt(2))^4 )^(-1)     ) / 4  =  2.71428571429
    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"2.71429", "0"});
}

TEST_CASE("MeanInverseParticipationRatio: calculating mean") {
    MeanInverseParticipationRatio ratioCalculator(0.5, 0.3);
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // Here: energy 0.4 with ipr = 1
    ratioCalculator.analyze(Eigensystem({0, 0.4, 0.9, 1}, arma::eye(4, 4)), logger);
    // Here: energies 0.5, 0.6 with avg ipr = 1.5
    ratioCalculator.analyze(Eigensystem({         0,        0.5, 0.6, 1},
                                        {{M_SQRT1_2, -M_SQRT1_2,   0, 0},
                                         {M_SQRT1_2,  M_SQRT1_2,   0, 0},
                                         {        0,          0,   1, 0},
                                         {        0,          0,   0, 1}}), logger);

    REQUIRE(ratioCalculator.getResultFields() == std::vector<std::string>{"1.2500", "0.2500"});
}

TEST_CASE("MeanInverseParticipationRatio: requires eigenvectors") {
    std::ostringstream loggerStream;
    Logger logger(loggerStream);
    REQUIRE_THROWS_WITH(MeanInverseParticipationRatio(0.5, 0.1).analyze(Eigensystem({0, 1, 2, 3}), logger),
                        Catch::Contains("hasEigenvectors"));
}