//
// Created by pkua on 01.11.2019.
//

#include <catch2/catch.hpp>

#include "mocks/DiagonalTermMock.h"
#include "mocks/HoppingTermMock.h"
#include "mocks/DoubleHoppingTermMock.h"

#include "matchers/ArmaApproxEqualCatchMatcher.h"

#include "core/FockBasisGenerator.h"
#include "core/HamiltonianGenerator.h"
#include "utils/Assertions.h"

using namespace trompeloeil;

namespace {
    class HopBetween {
    private:
        std::size_t a{};
        std::size_t b{};

    public:
        HopBetween(std::size_t a, std::size_t b) : a{a}, b{b} { }
        HopBetween(const HopData &hopData) : a{hopData.fromSite}, b{hopData.toSite} { }

        bool operator==(HopBetween other) {
            return (this->a == other.a && this->b == other.b) || (this->a == other.b && this->b == other.a);
        }
    };
}

TEST_CASE("HamiltonianGenerator: non-periodic hop") {
    auto hopping = std::make_unique<HoppingTermMock>();
    ALLOW_CALL(*hopping, calculate(_, _))
         .WITH(HopBetween(_1) == HopBetween(0, 1))
        .RETURN(1);
    ALLOW_CALL(*hopping, calculate(_, _))
        .WITH(HopBetween(_1) == HopBetween(1, 2))
        .RETURN(1);
    FockBasisGenerator baseGenerator;
    auto fockBase = baseGenerator.generate(2, 3);
    HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), false);
    hamiltonianGenerator.addHoppingTerm(std::move(hopping));

    arma::mat result = arma::mat(hamiltonianGenerator.generate());

    arma::mat expected = {{0,       M_SQRT2, 0, 0,       0,       0},
                          {M_SQRT2, 0,       1, M_SQRT2, 0,       0},
                          {0,       1,       0, 0,       1,       0},
                          {0,       M_SQRT2, 0, 0,       M_SQRT2, 0},
                          {0,       0,       1, M_SQRT2, 0,       M_SQRT2},
                          {0,       0,       0, 0,       M_SQRT2, 0}};
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("HamiltonianGenerator: double hop 2 on 2") {
    auto hopping = std::make_unique<DoubleHoppingTermMock>();
    ALLOW_CALL(*hopping, calculate(_, _, _))
            .WITH(_3.getSiteDistance(_1.fromSite, _1.toSite) == 1 && _3.getSiteDistance(_2.fromSite, _2.toSite) == 1)
            .RETURN(1);
    FockBasisGenerator baseGenerator;
    auto fockBase = baseGenerator.generate(2, 2);
    HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), false);
    hamiltonianGenerator.addDoubleHoppingTerm(std::move(hopping));

    arma::mat result = arma::mat(hamiltonianGenerator.generate());

    arma::mat expected = {{2, 0, 2},
                          {0, 4, 0},
                          {2, 0, 2}};
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("HamiltonianGenerator: double hop 2 on 3") {
    auto hopping = std::make_unique<DoubleHoppingTermMock>();
    ALLOW_CALL(*hopping, calculate(_, _, _))
            .WITH(_3.getSiteDistance(_1.fromSite, _1.toSite) == 1 && _3.getSiteDistance(_2.fromSite, _2.toSite) == 1)
            .RETURN(1);
    FockBasisGenerator baseGenerator;
    auto fockBase = baseGenerator.generate(2, 3);
    HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), false);
    hamiltonianGenerator.addDoubleHoppingTerm(std::move(hopping));

    arma::mat result = arma::mat(hamiltonianGenerator.generate());

    arma::mat expected = {{      2, 0,   M_SQRT2,         2, 0,       0},
                          {      0, 5,         0,         0, 3,       0},
                          {M_SQRT2, 0,         2, 2*M_SQRT2, 0, M_SQRT2},
                          {      2, 0, 2*M_SQRT2,         4, 0,       2},
                          {      0, 3,         0,         0, 5,       0},
                          {      0, 0,   M_SQRT2,         2, 0,       2}};
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("HamiltonianGenerator: diagonal") {
    auto diagonal = std::make_unique<DiagonalTermMock>();
    ALLOW_CALL(*diagonal, calculate(_, _))
        .RETURN(_2.getFockBasis()->findIndex(_1).value());
    FockBasisGenerator baseGenerator;
    auto fockBase = baseGenerator.generate(2, 3);
    HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), true);
    hamiltonianGenerator.addDiagonalTerm(std::move(diagonal));

    arma::mat result = arma::mat(hamiltonianGenerator.generate());

    arma::mat expected = arma::diagmat(arma::vec{0, 1, 2, 3, 4, 5});
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("HamiltonianGenerator: diagonalization") {
    SECTION("diagonal") {
        // Simple {2, 1} on diagonal - it should only get the order {1, 2} right
        auto diagonal = std::make_unique<DiagonalTermMock>();
        ALLOW_CALL(*diagonal, calculate(FockBasis::Vector{1, 0}, _))
            .RETURN(2);
        ALLOW_CALL(*diagonal, calculate(FockBasis::Vector{0, 1}, _))
            .RETURN(1);
        FockBasisGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(1, 2);
        HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), true);
        hamiltonianGenerator.addDiagonalTerm(std::move(diagonal));

        Eigensystem result = hamiltonianGenerator.calculateEigensystem(true);

        REQUIRE_THAT(result.getEigenenergies(), IsApproxEqual(arma::vec{1, 2}, 1e-8));
        REQUIRE_THAT(result.getEigenstates(), IsApproxEqual(arma::mat{{0, 1}, {1, 0}}, 1e-8));
    }

    SECTION("off-diagonal") {
        // Simple zeros on diagonal and ones off diagonal
        auto hopping = std::make_unique<HoppingTermMock>();
        ALLOW_CALL(*hopping, calculate(_, _))
                .WITH(HopBetween(_1) == HopBetween(0, 1))
                .RETURN(1);
        FockBasisGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(1, 2);
        HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), false);
        hamiltonianGenerator.addHoppingTerm(std::move(hopping));

        Eigensystem result = hamiltonianGenerator.calculateEigensystem(true);

        REQUIRE_THAT(result.getEigenenergies(), IsApproxEqual(arma::vec{-1, 1}, 1e-8));
        // Allow for +-eigenvectors
        REQUIRE_THAT(result.getEigenstates().col(0), IsApproxEqual(arma::vec{-M_SQRT1_2, M_SQRT1_2}, 1e-8)
                                                     || IsApproxEqual(arma::vec{M_SQRT1_2, -M_SQRT1_2}, 1e-8));
        REQUIRE_THAT(result.getEigenstates().col(1), IsApproxEqual(arma::vec{M_SQRT1_2, M_SQRT1_2}, 1e-8)
                                                     || IsApproxEqual(arma::vec{-M_SQRT1_2, -M_SQRT1_2}, 1e-8));
    }
}

TEST_CASE("HamiltonianGenerator: PBC") {
    SECTION("Periodic BC - periodic hopping should be present") {
        auto hopping = std::make_unique<HoppingTermMock>();
        ALLOW_CALL(*hopping, calculate(_, _))
                .WITH(HopBetween(_1) == HopBetween(0, 1))
                .RETURN(-1);
        ALLOW_CALL(*hopping, calculate(_, _))
                .WITH(HopBetween(_1) == HopBetween(1, 2))
                .RETURN(-1);
        ALLOW_CALL(*hopping, calculate(_, _))
                .WITH(HopBetween(_1) == HopBetween(0, 2))
                .RETURN(-1);
        FockBasisGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(2, 3);
        HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), true);
        hamiltonianGenerator.addHoppingTerm(std::move(hopping));

        REQUIRE_NOTHROW(hamiltonianGenerator.generate());
    }

    SECTION("Non periodic BC - periodic hopping should not be present") {
        auto hopping = std::make_unique<HoppingTermMock>();
        ALLOW_CALL(*hopping, calculate(_, _))
                .WITH(HopBetween(_1) == HopBetween(0, 1))
                .RETURN(-1);
        ALLOW_CALL(*hopping, calculate(_, _))
                .WITH(HopBetween(_1) == HopBetween(1, 2))
                .RETURN(-1);
        FockBasisGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(2, 3);
        HamiltonianGenerator hamiltonianGenerator(std::move(fockBase), false);
        hamiltonianGenerator.addHoppingTerm(std::move(hopping));

        REQUIRE_NOTHROW(hamiltonianGenerator.generate());
    }
}

TEST_CASE("HamiltonianGenerator: site distance") {
    FockBasisGenerator baseGenerator;
    auto evenBase = baseGenerator.generate(1, 6);
    auto oddBase = baseGenerator.generate(1, 7);

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