//
// Created by Piotr Kubala on 22/01/2020.
//

#include <catch2/catch.hpp>

#include "matchers/ArmaApproxEqualCatchMatcher.h"

#include "core/Eigensystem.h"

TEST_CASE("Eigensystem: empty") {
    Eigensystem eigensystem;

    CHECK(eigensystem.empty());
    CHECK(eigensystem.size() == 0);
    CHECK_FALSE(eigensystem.hasEigenvectors());
    CHECK_THROWS(eigensystem.getEigenstates());
    CHECK_THAT(eigensystem.getEigenenergies(), IsApproxEqual(arma::vec{}, 1e-15));
    CHECK_THROWS(eigensystem.getEigenstate(0));
}

TEST_CASE("Eigensystem: only eigenvalues") {
    Eigensystem eigensystem({0, 0.5, 1});

    CHECK_FALSE(eigensystem.empty());
    REQUIRE(eigensystem.size() == 3);
    CHECK_FALSE(eigensystem.hasEigenvectors());
    CHECK_THROWS(eigensystem.getEigenstates());
    CHECK_THAT(eigensystem.getEigenenergies(), IsApproxEqual(arma::vec{0, 0.5, 1}, 1e-15));
    CHECK_THROWS(eigensystem.getEigenstate(0));
}

TEST_CASE("Eigensystem: eigenvalues and eigenvectors") {
    arma::vec eigenenergies{0, 0.5, 1};
    arma::mat eigenstates{{0, 1, -0.5},
                          {-M_SQRT1_2, 0, 0.5},
                          {M_SQRT1_2, 0, M_SQRT1_2}};
    Eigensystem eigensystem(eigenenergies, eigenstates);

    CHECK_FALSE(eigensystem.empty());
    REQUIRE(eigensystem.size() == 3);
    REQUIRE(eigensystem.hasEigenvectors());
    CHECK_THAT(eigensystem.getEigenstates(), IsApproxEqual(eigenstates, 1e-16));
    CHECK_THAT(eigensystem.getEigenenergies(), IsApproxEqual(eigenenergies, 1e-15));
    CHECK_THAT(eigensystem.getEigenstate(0), IsApproxEqual(arma::vec{0, -M_SQRT1_2, M_SQRT1_2}, 1e-15));
    CHECK_THAT(eigensystem.getEigenstate(1), IsApproxEqual(arma::vec{1, 0, 0}, 1e-15));
    CHECK_THAT(eigensystem.getEigenstate(2), IsApproxEqual(arma::vec{-0.5, 0.5, M_SQRT1_2}, 1e-15));
}

TEST_CASE("Eigensystem: eigenvalues sorting") {
    arma::vec eigenenergies{0.5, 1, 0};
    arma::mat eigenstates{{         0, 1,      -0.5},
                          {-M_SQRT1_2, 0,       0.5},
                          { M_SQRT1_2, 0, M_SQRT1_2}};
    Eigensystem eigensystem(eigenenergies, eigenstates);

    arma::mat expectedEigenstates{{     -0.5,          0, 1},
                                  {      0.5, -M_SQRT1_2, 0},
                                  {M_SQRT1_2,  M_SQRT1_2, 0}};
    REQUIRE_THAT(eigensystem.getEigenstates(), IsApproxEqual(expectedEigenstates, 1e-16));
}

TEST_CASE("Eigensystem: eigenvectors normalization") {
    arma::vec eigenenergies{0, 0.5, 1};
    arma::mat eigenstates{{ 0, 3,      -1},
                          {-1, 0,       1},
                          { 1, 0, M_SQRT2}};
    Eigensystem eigensystem(eigenenergies, eigenstates);

    arma::mat expectedEigenstates{{0,          1, -0.5},
                                  {-M_SQRT1_2, 0, 0.5},
                                  { M_SQRT1_2, 0, M_SQRT1_2}};
    REQUIRE(eigensystem.size() == 3);
    REQUIRE_THAT(eigensystem.getEigenstates(), IsApproxEqual(expectedEigenstates, 1e-15));
}

TEST_CASE("Eigensystem: orthonormality check") {
    SECTION("not orthonormal") {
        Eigensystem eigensystem({1, 2, 3, 4}, {{1, 1,  1,  1},
                                               {1, -1, 0,  0},
                                               {1, 1,  -2, 0},
                                               {1, 1,  1,  -3}});

        REQUIRE_FALSE(eigensystem.isOrthonormal());
    }

    SECTION("orthonormal") {
        const double sq2 = std::sqrt(2);
        const double sq6 = std::sqrt(6);
        const double sq12 = std::sqrt(12);
        Eigensystem eigensystem({1, 2, 3, 4}, {{   0.5,    0.5,   0.5,      0.5},
                                               { 1/sq2, -1/sq2,     0,        0},
                                               { 1/sq6,  1/sq6, -2/sq6,       0},
                                               {1/sq12, 1/sq12, 1/sq12, -3/sq12}});

        REQUIRE(eigensystem.isOrthonormal());
    }

    SECTION("no eigenvectors") {
        Eigensystem eigensystem({1, 2, 3});

        REQUIRE_THROWS(eigensystem.isOrthonormal());
    }
}

TEST_CASE("Eigensystem: store/restore") {
    SECTION("only eigenenergies") {
        Eigensystem toStore({0, 0.25, 0.5, 1});
        Eigensystem toRestore({0, 0.1, 0.8, 1});

        std::stringstream inout;
        toStore.store(inout);
        toRestore.restore(inout);

        CHECK(toStore == toRestore);
        CHECK_FALSE(toRestore.hasEigenvectors());
    }

    SECTION("eigenenergies and eigenstates") {
        Eigensystem toStore({0, 0.25, 0.5, 1},
                            {{0, -1,         0,          0},
                             {1,  0,         0,          0},
                             {0,  0, M_SQRT1_2, -M_SQRT1_2},
                             {0,  0, M_SQRT1_2,  M_SQRT1_2}});
        Eigensystem toRestore({0, 0.1, 0.8, 1});

        std::stringstream inoutEnergies, inoutStates;
        toStore.store(inoutEnergies, inoutStates);
        toRestore.restore(inoutEnergies, inoutStates);

        REQUIRE(toStore == toRestore);
        REQUIRE(toRestore.hasEigenvectors());
    }
}

TEST_CASE("Eigensystem: normalized energies") {
    Eigensystem eigensystem({1, 2, 3});

    REQUIRE_THAT(eigensystem.getNormalizedEigenenergies(), IsApproxEqual(arma::vec{0, 0.5, 1}, 1e-15));
}

TEST_CASE("Eigensystem: unmatching sizes") {
    CHECK_THROWS(Eigensystem({1, 2}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}));
    CHECK_THROWS(Eigensystem({1, 2}, {{1, 2}, {4, 5}, {7, 8}}));
}

TEST_CASE("Eigensystem: zero vector") {
    CHECK_THROWS(Eigensystem({0, 1}, {{0, 2}, {0, 1}}));
}

TEST_CASE("Eigensystem: cannot normalize equal eigenvalues") {
    CHECK_THROWS(Eigensystem({2, 2, 2}).getNormalizedEigenenergies());
}

TEST_CASE("Eigensystem: get indices of eigenenergies in band") {
    Eigensystem eigensystem({0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0});

    SECTION("range inside [0, 1]") {
        auto indices = eigensystem.getIndicesOfNormalizedEnergiesInBand(0.5, 0.5);

        REQUIRE(indices == std::vector<std::size_t>{2, 3, 4, 5});
    }

    SECTION("range sticking outside [0, 1]") {
        SECTION("below 0") {
            auto indices = eigensystem.getIndicesOfNormalizedEnergiesInBand(0.1, 0.3);

            REQUIRE(indices == std::vector<std::size_t>{0, 1});
        }

        SECTION("over 1") {
            auto indices = eigensystem.getIndicesOfNormalizedEnergiesInBand(0.95, 0.2);

            REQUIRE(indices == std::vector<std::size_t>{7, 8});
        }
    }
}

TEST_CASE("Eigensystem: get indices of eigenenergies in band - incorrect") {
    Eigensystem eigensystem({0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0});

    CHECK_THROWS(eigensystem.getIndicesOfNormalizedEnergiesInBand(-0.3, 0.1));
    CHECK_THROWS(eigensystem.getIndicesOfNormalizedEnergiesInBand(1.2, 0.1));
    CHECK_THROWS(eigensystem.getIndicesOfNormalizedEnergiesInBand(0.4, 0));
    CHECK_THROWS(eigensystem.getIndicesOfNormalizedEnergiesInBand(0.4, -0.1));
}

TEST_CASE("Eigensystem: FockBase") {
    auto base = std::make_shared<FockBasis>();
    base->add({2, 0});
    base->add({1, 1});
    base->add({0, 2});

    SECTION("valid construction") {
        CHECK_NOTHROW(Eigensystem({0, 0.5, 1}, base));
        CHECK(Eigensystem({0, 0.5, 1}, base).hasFockBasis());
        CHECK_FALSE(Eigensystem({0, 0.5, 1}).hasFockBasis());
        CHECK_NOTHROW(Eigensystem({0, 0.5, 1}, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, base));
    }

    SECTION("invalid construction") {
        CHECK_THROWS(Eigensystem({0, 1}, base));
        CHECK_THROWS(Eigensystem({0, 1}, {{1, 0, }, {0, 1}}, base));
    }

    SECTION("restore") {
        Eigensystem eigensystem({0, 1});
        Eigensystem saved({0, 0.5, 1});
        std::stringstream stream;
        saved.store(stream);

        CHECK_NOTHROW(eigensystem.restore(stream, base));
    }
}