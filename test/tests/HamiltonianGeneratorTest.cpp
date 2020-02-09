//
// Created by pkua on 01.11.2019.
//

#include <catch2/catch.hpp>

#include "mocks/DiagonalTermMock.h"
#include "mocks/HoppingTermMock.h"

#include "matchers/ArmaApproxEqualMatcher.h"

#include "simulation/FockBaseGenerator.h"
#include "simulation/HamiltonianGenerator.h"
#include "utils/Assertions.h"

using namespace trompeloeil;

namespace {
    class HopBetween {
    private:
        std::size_t a{};
        std::size_t b{};

    public:
        HopBetween(std::size_t a, std::size_t b) : a{a}, b{b} {}

        bool operator==(HopBetween other) {
            return (this->a == other.a && this->b == other.b) || (this->a == other.b && this->b == other.a);
        }
    };
}

TEST_CASE("HamiltonianGenerator: non-periodic hop") {
    auto hopping = std::make_unique<HoppingTermMock>();
    REQUIRE_CALL(*hopping, calculate(_, _, _, _, _))
        .WITH(HopBetween(_3, _4) == HopBetween(0, 1))
        .RETURN(-1)
        .TIMES(AT_LEAST(1));
    REQUIRE_CALL(*hopping, calculate(_, _, _, _, _))
        .WITH(HopBetween(_3, _4) == HopBetween(1, 2))
        .RETURN(-1)
        .TIMES(AT_LEAST(1));
    FockBaseGenerator baseGenerator;
    auto fockBase = baseGenerator.generate(3, 2);
    HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), false);
    hamiltonianGenerator.addHoppingTerm(std::move(hopping));

    arma::mat result = hamiltonianGenerator.generate();

    arma::mat expected = {{ 0,       -M_SQRT2,  0,  0,        0,        0},
                          {-M_SQRT2,  0,       -1, -M_SQRT2,  0,        0},
                          { 0,       -1,        0,  0,       -1,        0},
                          { 0,       -M_SQRT2,  0,  0,       -M_SQRT2,  0},
                          { 0,        0,       -1, -M_SQRT2,  0,       -M_SQRT2},
                          { 0,        0,        0,  0,       -M_SQRT2,  0}};
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("HamiltonianGenerator: diagonal") {
    auto diagonal = std::make_unique<DiagonalTermMock>();
    ALLOW_CALL(*diagonal, calculate(_, _))
        .RETURN(_2.getFockBase().findIndex(_1).value());
    FockBaseGenerator baseGenerator;
    auto fockBase = baseGenerator.generate(3, 2);
    HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), true);
    hamiltonianGenerator.addDiagonalTerm(std::move(diagonal));

    arma::mat result = hamiltonianGenerator.generate();

    arma::mat expected = arma::diagmat(arma::vec{0, 1, 2, 3, 4, 5});
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("HamiltonianGenerator: PBC") {
    SECTION("Periodic BC - periodic hopping should be present") {
        auto hopping = std::make_unique<HoppingTermMock>();
        REQUIRE_CALL(*hopping, calculate(_, _, _, _, _))
                .WITH(HopBetween(_3, _4) == HopBetween(0, 1))
                .RETURN(-1)
                .TIMES(AT_LEAST(1));
        REQUIRE_CALL(*hopping, calculate(_, _, _, _, _))
                .WITH(HopBetween(_3, _4) == HopBetween(1, 2))
                .RETURN(-1)
                .TIMES(AT_LEAST(1));
        REQUIRE_CALL(*hopping, calculate(_, _, _, _, _))
                .WITH(HopBetween(_3, _4) == HopBetween(0, 2))
                .RETURN(-1)
                .TIMES(AT_LEAST(1));
        FockBaseGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(3, 2);
        HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), true);
        hamiltonianGenerator.addHoppingTerm(std::move(hopping));

        REQUIRE_NOTHROW(hamiltonianGenerator.generate());
    }

    SECTION("Non periodic BC - periodic hopping should not be present") {
        auto hopping = std::make_unique<HoppingTermMock>();
        REQUIRE_CALL(*hopping, calculate(_, _, _, _, _))
                .WITH(HopBetween(_3, _4) == HopBetween(0, 1))
                .RETURN(-1)
                .TIMES(AT_LEAST(1));
        REQUIRE_CALL(*hopping, calculate(_, _, _, _, _))
                .WITH(HopBetween(_3, _4) == HopBetween(1, 2))
                .RETURN(-1)
                .TIMES(AT_LEAST(1));
        FockBaseGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(3, 2);
        HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), false);
        hamiltonianGenerator.addHoppingTerm(std::move(hopping));

        REQUIRE_NOTHROW(hamiltonianGenerator.generate());
    }
}

TEST_CASE("HamiltonianGenerator: site distance") {
    FockBaseGenerator baseGenerator;
    auto evenBase = baseGenerator.generate(6, 1);
    auto oddBase = baseGenerator.generate(7, 1);

    SECTION("Periodic BC - site distance calculated normally") {
        SECTION("even size") {
            HamiltonianGenerator hamiltonianGenerator(std::move(evenBase), true);

            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 1) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(4, 5) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 3) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 4) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 5) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 5) == 1);
        }

        SECTION("odd size") {
            HamiltonianGenerator hamiltonianGenerator(std::move(oddBase), true);

            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 1) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(5, 6) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 3) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 4) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 5) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 6) == 1);
        }
    }

    SECTION("Non-periodic BC - site distance calculated periodically") {
        SECTION("even size") {
            HamiltonianGenerator hamiltonianGenerator(std::move(evenBase), false);

            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 1) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(4, 5) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 3) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 4) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 5) == 4);
            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 5) == 5);
        }

        SECTION("odd size") {
            HamiltonianGenerator hamiltonianGenerator(std::move(oddBase), false);

            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 1) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(5, 6) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 3) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 4) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 5) == 4);
            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 6) == 6);
        }
    }
}