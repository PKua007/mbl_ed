//
// Created by Piotr Kubala on 17/01/2020.
//

#include <catch2/catch.hpp>
#include <sstream>

#include "analyzer/tasks/CDF.h"

TEST_CASE("CDF: name") {
    CDF cdf(10);

    REQUIRE(cdf.getName() == "cdf");
}

TEST_CASE("CDF: single set") {
    SECTION("already normalized") {
        CDF cdf(6);
        std::ostringstream logger;
        cdf.analyze(Eigensystem({0, 0.05, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8, 0.9, 1.0}), logger);

        std::ostringstream out;
        cdf.storeResult(out);

        std::ostringstream expected;
        expected << "0 0.1" << std::endl << "0.2 0.4" << std::endl << "0.4 0.5" << std::endl << "0.6 0.6" << std::endl;
        expected << "0.8 0.8" << std::endl << "1 1" << std::endl;
        REQUIRE(out.str() == expected.str());
    }

    SECTION("normalization") {
        CDF cdf(6);
        std::ostringstream logger;
        cdf.analyze(Eigensystem({1, 1.5, 2, 3, 5, 6, 8, 9, 10, 11}), logger);

        std::ostringstream out;
        cdf.storeResult(out);

        std::ostringstream expected;
        expected << "0 0.1" << std::endl << "0.2 0.4" << std::endl << "0.4 0.5" << std::endl << "0.6 0.6" << std::endl;
        expected << "0.8 0.8" << std::endl << "1 1" << std::endl;
        REQUIRE(out.str() == expected.str());
    }
}

TEST_CASE("CDF: two sets") {
    CDF cdf(6);
    std::ostringstream logger;
    cdf.analyze(Eigensystem({0, 0.05, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8, 0.9, 1.0}), logger);
    cdf.analyze(Eigensystem({0, 0.1, 0.3, 0.4, 0.5, 0.6, 0.8, 0.85, 0.95, 1.0}), logger);

    std::ostringstream out;
    cdf.storeResult(out);

    std::ostringstream expected;
    expected << "0 0.1" << std::endl << "0.2 0.3" << std::endl << "0.4 0.45" << std::endl << "0.6 0.6" << std::endl;
    expected << "0.8 0.75" << std::endl << "1 1" << std::endl;
    REQUIRE(out.str() == expected.str());
}

TEST_CASE("CDF: corner cases") {
    SECTION("only 2 energies") {
        CDF cdf(3);
        std::ostringstream logger;
        cdf.analyze(Eigensystem({0, 1}), logger);

        std::ostringstream out;
        cdf.storeResult(out);

        std::ostringstream expected;
        expected << "0 0.5" << std::endl << "0.5 0.5" << std::endl << "1 1" << std::endl;
        REQUIRE(out.str() == expected.str());
    }

    SECTION("only 2 bins") {
        CDF cdf(2);
        std::ostringstream logger;
        cdf.analyze(Eigensystem({0, 0.1, 0.8, 1}), logger);

        std::ostringstream out;
        cdf.storeResult(out);

        std::ostringstream expected;
        expected << "0 0.25" << std::endl << "1 1" << std::endl;
        REQUIRE(out.str() == expected.str());
    }

    SECTION("0 or 1 bins should throw") {
        REQUIRE_THROWS(CDF(0));
        REQUIRE_THROWS(CDF(1));
    }
}

TEST_CASE("CDF: storing, restoring and cleaning") {
    CDF cdf(3);
    Eigensystem eigensystem0({0, 0.1, 0.15, 0.16, 0.2, 0.25, 0.40, 0.47, 0.58, 0.71, 0.76, 0.8, 0.94, 1.0});
    Eigensystem eigensystem1({0, 0.05, 0.25, 0.27, 0.28, 0.32, 0.45, 0.58, 0.59, 0.6, 0.75, 0.85, 0.91, 1.0});
    std::ostringstream logger;

    SECTION("clearing") {
        cdf.analyze(eigensystem0, logger);
        std::ostringstream result1;
        cdf.storeResult(result1);

        cdf.analyze(eigensystem1, logger);
        cdf.clear();
        cdf.analyze(eigensystem0, logger);
        std::ostringstream result2;
        cdf.storeResult(result2);

        REQUIRE(result1.str() == result2.str());
    }

    SECTION("storing and joining restored") {
        cdf.analyze(eigensystem0, logger);
        cdf.analyze(eigensystem1, logger);
        std::ostringstream normalResult;
        cdf.storeResult(normalResult);

        cdf.clear();
        cdf.analyze(eigensystem1, logger);
        std::stringstream simulation1;
        cdf.storeState(simulation1);

        cdf.clear();
        cdf.analyze(eigensystem0, logger);
        cdf.joinRestoredState(simulation1);
        std::ostringstream restoredResult;
        cdf.storeResult(restoredResult);

        REQUIRE(normalResult.str() == restoredResult.str());
    }
}