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
        cdf.analyze({0, 0.05, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8, 0.9, 1.0});

        std::ostringstream out;
        cdf.storeResult(out);

        std::ostringstream expected;
        expected << "0 0.1" << std::endl << "0.2 0.4" << std::endl << "0.4 0.5" << std::endl << "0.6 0.6" << std::endl;
        expected << "0.8 0.8" << std::endl << "1 1" << std::endl;
        REQUIRE(out.str() == expected.str());
    }

    SECTION("normalization") {
        CDF cdf(6);
        cdf.analyze({1, 1.5, 2, 3, 5, 6, 8, 9, 10, 11});

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
    cdf.analyze({0, 0.05, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8, 0.9, 1.0});
    cdf.analyze({0, 0.1, 0.3, 0.4, 0.5, 0.6, 0.8, 0.85, 0.95, 1.0});

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
        cdf.analyze({0, 1});

        std::ostringstream out;
        cdf.storeResult(out);

        std::ostringstream expected;
        expected << "0 0.5" << std::endl << "0.5 0.5" << std::endl << "1 1" << std::endl;
        REQUIRE(out.str() == expected.str());
    }

    SECTION("only 2 bins") {
        CDF cdf(2);
        cdf.analyze({0, 0.1, 0.8, 1});

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