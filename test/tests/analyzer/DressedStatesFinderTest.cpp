//
// Created by Piotr Kubala on 19/07/2020.
//

#include <catch2/catch.hpp>

#include "analyzer/tasks/DressedStatesFinder.h"
#include "core/FockBasisGenerator.h"

#include "core/FockBasisGenerator.h"
#include "core/terms/HubbardHop.h"
#include "core/terms/HubbardOnsite.h"
#include "core/terms/QuasiperiodicDisorder.h"
#include "core/averaging_models/UniformPhi0AveragingModel.h"

TEST_CASE("DressedStatesFinder: single diagonalization") {
    FockBasisGenerator generator;
    auto base = generator.generate(4, 2);
    Eigensystem eigensystem({0, 0.25, 0.5, 0.75, 1},
                            {{1,   0,    0,       0,        0},
                             {0, 0.6, -0.8,       0,        0},
                             {0, 0.8,  0.6,       0,        0},
                             {0,   0,    0, M_SQRT2, -M_SQRT2},
                             {0,   0,    0, M_SQRT2,  M_SQRT2}},
                            std::move(base));
    DressedStatesFinder finder(0.75, BandExtractor::EpsilonRange(0.5, 0.8));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // To sum up: we will investigate eigenvalues 0.25, 0.5, 0.75 and only first 2 will be ok (above 0.75 threshold)
    // 0.25 -> 0.8 on 3rd coeff, so 2.2
    // 0.5 -> -0.8 on 2nd coeff, so 3.1   <= here is negative - negative ones are also ok
    finder.analyze(eigensystem, logger);

    std::ostringstream out;
    finder.storeResult(out);
    REQUIRE(out.str() == "0 2.2 0.25 0.8\n0 3.1 0.5 -0.8\n");
}

TEST_CASE("DressedStatesFinder: 2 diagonalizations") {
    FockBasisGenerator generator;
    auto base = std::shared_ptr<FockBasis>(generator.generate(2, 2));
    Eigensystem eigensystem1({0, 0.5, 1},
                             {{1, 0, 0},
                              {0, 1, 0},
                              {0, 0, 1}},
                             base);
    Eigensystem eigensystem2({0, 0.6, 1},
                             {{0, -1, 0},
                              {1,  0, 0},
                              {0,  0, 1}},
                             base);
    DressedStatesFinder finder(0.8, BandExtractor::EpsilonRange(0.5, 0.3));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    finder.analyze(eigensystem1, logger);
    finder.analyze(eigensystem2, logger);

    std::ostringstream out;
    finder.storeResult(out);
    REQUIRE(out.str() == "0 1.1 0.5 1\n1 2.0 0.6 -1\n");
}