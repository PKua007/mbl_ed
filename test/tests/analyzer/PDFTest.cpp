//
// Created by Piotr Kubala on 01/11/2020.
//

#include <catch2/catch.hpp>
#include <sstream>

#include "analyzer/tasks/PDF.h"

TEST_CASE("PDF: name") {
    PDF pdf(10);

    REQUIRE(pdf.getName() == "pdf");
}

TEST_CASE("PDF: single set") {
    PDF pdf(5);
    std::ostringstream loggerStream;
    Logger logger(loggerStream);
    pdf.analyze(Eigensystem({0, 0.05, 0.1, 0.21, 0.41, 0.5, 0.7, 0.81, 0.9, 1.0}), logger);

    std::ostringstream out;
    pdf.storeResult(out);

    std::ostringstream expected;
    expected << "0.1 0.3" << std::endl << "0.3 0.1" << std::endl << "0.5 0.2" << std::endl << "0.7 0.1" << std::endl;
    expected << "0.9 0.3" << std::endl;
    REQUIRE(out.str() == expected.str());
}

TEST_CASE("PDF: two sets") {
    PDF pdf(5);
    std::ostringstream loggerStream;
    Logger logger(loggerStream);
    pdf.analyze(Eigensystem({0, 0.05, 0.11, 0.21, 0.41, 0.5, 0.7, 0.81, 0.9, 1.0}), logger);
    pdf.analyze(Eigensystem({0, 0.1, 0.3, 0.41, 0.5, 0.61, 0.81, 0.85, 0.95, 1.0}), logger);

    std::ostringstream out;
    pdf.storeResult(out);

    std::ostringstream expected;
    expected << "0.1 0.25" << std::endl << "0.3 0.1" << std::endl << "0.5 0.2" << std::endl << "0.7 0.1" << std::endl;
    expected << "0.9 0.35" << std::endl;
    REQUIRE(out.str() == expected.str());
}

TEST_CASE("PDF: corner cases") {
    SECTION("only 2 bins") {
        PDF cdf(2);
        std::ostringstream loggerStream;
        Logger logger(loggerStream);
        cdf.analyze(Eigensystem({0, 0.1, 0.8, 1}), logger);

        std::ostringstream out;
        cdf.storeResult(out);

        std::ostringstream expected;
        expected << "0.25 0.5" << std::endl << "0.75 0.5" << std::endl;
        REQUIRE(out.str() == expected.str());
    }

    SECTION("0 or 1 bins should throw") {
        REQUIRE_THROWS(PDF(0));
        REQUIRE_THROWS(PDF(1));
    }
}