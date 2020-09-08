//
// Created by Piotr Kubala on 31/07/2020.
//

#include <catch2/catch.hpp>

#include "analyzer/tasks/BulkMeanGapRatio.h"

TEST_CASE("BulkMeanGapRatio: names") {
    BulkMeanGapRatio bmgr(4);

    REQUIRE(bmgr.getName() == "mgrs");
}

TEST_CASE("BulkMeanGapRatio: single set") {
    // bins: [0, 0.5], [0.5, 1.0]
    BulkMeanGapRatio bmgr(2);
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // for [0, 0.5] => mean{ (0.05 - 0)/(0.25 - 0.05) = 0.25, (0.25 - 0.05)/(0.75 - 0.25) = 0.4 } = 0.325
    // for [0.5, 1] => (1.0 - 0.75)/(0.75 - 0.25) = 0.5
    bmgr.analyze(Eigensystem({0, 0.05, 0.25, 0.75, 1.0}), logger);

    std::stringstream out;
    bmgr.storeResult(out);
    REQUIRE(out.str() == "0 0.32500 0.07500\n0.5 0.5 0\n");
}

TEST_CASE("BulkMeanGapRatio: 2 sets") {
    // bins: [0, 0.5], [0.5, 1.0]
    BulkMeanGapRatio bmgr(2);
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    // Two eigensystems with "opposite gaps", averaging in both bins to
    // mean{ (0.8 - 0.4)/(0.4 - 0) = 1, (1.0 - 0.8)/(0.8 - 0.4) = 0.5 } = 0.75
    bmgr.analyze(Eigensystem({0, 0.4, 0.8, 1.0}), logger);
    bmgr.analyze(Eigensystem({0, 0.2, 0.6, 1.0}), logger);

    std::stringstream out;
    bmgr.storeResult(out);
    REQUIRE(out.str() == "0 0.7500 0.2500\n0.5 0.7500 0.2500\n");
}