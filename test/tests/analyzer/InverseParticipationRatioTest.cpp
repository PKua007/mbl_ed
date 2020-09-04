//
// Created by Piotr Kubala on 30/01/2020.
//

#include <catch2/catch.hpp>

#include "analyzer/tasks/InverseParticipationRatio.h"

#include "core/FockBaseGenerator.h"
#include "core/terms/HubbardHop.h"
#include "core/terms/HubbardOnsite.h"
#include "core/terms/QuasiperiodicDisorder.h"
#include "core/averaging_models/UniformPhi0AveragingModel.h"

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

TEST_CASE("InverseParticipationRatio: storing, restoring and cleaning") {
    auto fockBase = std::shared_ptr(FockBaseGenerator().generate(4, 4));
    HamiltonianGenerator hamiltonianGenerator(fockBase, false);
    hamiltonianGenerator.addHoppingTerm(std::make_unique<HubbardHop>(1));
    hamiltonianGenerator.addDiagonalTerm(std::make_unique<HubbardOnsite>(1));
    hamiltonianGenerator.addDiagonalTerm(std::make_unique<QuasiperiodicDisorder>(1, 0.3, 0));
    UniformPhi0AveragingModel averagingModel{};
    RND rnd{};
    averagingModel.setupHamiltonianGenerator(hamiltonianGenerator, rnd, 1, 5);
    Eigensystem eigensystem1 = hamiltonianGenerator.calculateEigensystem(true);
    averagingModel.setupHamiltonianGenerator(hamiltonianGenerator, rnd, 2, 5);
    Eigensystem eigensystem2 = hamiltonianGenerator.calculateEigensystem(true);
    InverseParticipationRatio inverseParticipationRatio(0.5, 0.05);
    std::ostringstream logger;

    SECTION("clearing") {
        inverseParticipationRatio.analyze(eigensystem1, logger);
        std::ostringstream result1;
        inverseParticipationRatio.storeResult(result1);

        inverseParticipationRatio.analyze(eigensystem2, logger);
        inverseParticipationRatio.clear();
        inverseParticipationRatio.analyze(eigensystem1, logger);
        std::ostringstream result2;
        inverseParticipationRatio.storeResult(result2);

        REQUIRE(result1.str() == result2.str());
        std::cout << logger.str() << std::endl;
    }

    SECTION("storing and joining restored") {
        inverseParticipationRatio.analyze(eigensystem1, logger);
        inverseParticipationRatio.analyze(eigensystem2, logger);
        std::ostringstream normalResult;
        inverseParticipationRatio.storeResult(normalResult);

        inverseParticipationRatio.clear();
        inverseParticipationRatio.analyze(eigensystem2, logger);
        std::stringstream simulation2;
        inverseParticipationRatio.storeState(simulation2);

        inverseParticipationRatio.clear();
        inverseParticipationRatio.analyze(eigensystem1, logger);
        inverseParticipationRatio.joinRestoredState(simulation2);
        std::ostringstream restoredResult;
        inverseParticipationRatio.storeResult(restoredResult);

        REQUIRE(normalResult.str() == restoredResult.str());
    }
}