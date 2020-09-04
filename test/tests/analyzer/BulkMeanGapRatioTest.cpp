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
    std::ostringstream logger;

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
    std::ostringstream logger;

    // Two eigensystems with "opposite gaps", averaging in both bins to
    // mean{ (0.8 - 0.4)/(0.4 - 0) = 1, (1.0 - 0.8)/(0.8 - 0.4) = 0.5 } = 0.75
    bmgr.analyze(Eigensystem({0, 0.4, 0.8, 1.0}), logger);
    bmgr.analyze(Eigensystem({0, 0.2, 0.6, 1.0}), logger);

    std::stringstream out;
    bmgr.storeResult(out);
    REQUIRE(out.str() == "0 0.7500 0.2500\n0.5 0.7500 0.2500\n");
}

TEST_CASE("BulkMeanGapRatio: storing, restoring and cleaning") {
    BulkMeanGapRatio bulkMeanGapRatio(3);
    Eigensystem eigensystem0({0, 0.1, 0.15, 0.16, 0.2, 0.25, 0.40, 0.47, 0.58, 0.71, 0.76, 0.8, 0.94, 1.0});
    Eigensystem eigensystem1({0, 0.05, 0.25, 0.27, 0.28, 0.32, 0.45, 0.58, 0.59, 0.6, 0.75, 0.85, 0.91, 1.0});
    std::ostringstream logger;

    SECTION("clearing") {
        bulkMeanGapRatio.analyze(eigensystem0, logger);
        std::ostringstream result1;
        bulkMeanGapRatio.storeResult(result1);

        bulkMeanGapRatio.analyze(eigensystem1, logger);
        bulkMeanGapRatio.clear();
        bulkMeanGapRatio.analyze(eigensystem0, logger);
        std::ostringstream result2;
        bulkMeanGapRatio.storeResult(result2);

        REQUIRE(result1.str() == result2.str());
    }

    SECTION("storing and joining restored") {
        bulkMeanGapRatio.analyze(eigensystem0, logger);
        bulkMeanGapRatio.analyze(eigensystem1, logger);
        std::ostringstream normalResult;
        bulkMeanGapRatio.storeResult(normalResult);

        bulkMeanGapRatio.clear();
        bulkMeanGapRatio.analyze(eigensystem1, logger);
        std::stringstream simulation1;
        bulkMeanGapRatio.storeState(simulation1);

        bulkMeanGapRatio.clear();
        bulkMeanGapRatio.analyze(eigensystem0, logger);
        bulkMeanGapRatio.joinRestoredState(simulation1);
        std::ostringstream restoredResult;
        bulkMeanGapRatio.storeResult(restoredResult);

        REQUIRE(normalResult.str() == restoredResult.str());
    }
}