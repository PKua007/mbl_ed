//
// Created by Piotr Kubala on 22/01/2020.
//

#include <catch2/catch.hpp>

#include "matchers/ArmaApproxEqualMatcher.h"

#include "simulation/Eigensystem.h"

TEST_CASE("Eigensystem: empty") {
    Eigensystem eigensystem;

    REQUIRE(eigensystem.empty());
    REQUIRE(eigensystem.size() == 0);
    REQUIRE_FALSE(eigensystem.hasEigenvectors());
    REQUIRE_THAT(eigensystem.getEigenstates(), IsApproxEqual(arma::mat{}, 1e-15));
    REQUIRE_THAT(eigensystem.getEigenenergies(), IsApproxEqual(arma::vec{}, 1e-15));
    REQUIRE_THROWS(eigensystem.getEigenstate(0));
}

TEST_CASE("Eigensystem: only eigenvalues") {
    Eigensystem eigensystem({0, 0.5, 1});

    REQUIRE_FALSE(eigensystem.empty());
    REQUIRE(eigensystem.size() == 3);
    REQUIRE_FALSE(eigensystem.hasEigenvectors());
    REQUIRE_THAT(eigensystem.getEigenstates(), IsApproxEqual(arma::mat(3, 3, arma::fill::zeros), 1e-15));
    REQUIRE_THAT(eigensystem.getEigenenergies(), IsApproxEqual(arma::vec{0, 0.5, 1}, 1e-15));
    REQUIRE_THAT(eigensystem.getEigenstate(0), IsApproxEqual(arma::vec{0, 0, 0}, 1e-15));
    REQUIRE_THAT(eigensystem.getEigenstate(1), IsApproxEqual(arma::vec{0, 0, 0}, 1e-15));
    REQUIRE_THAT(eigensystem.getEigenstate(2), IsApproxEqual(arma::vec{0, 0, 0}, 1e-15));
}

TEST_CASE("Eigensystem: eigenvalues and eigenvectors") {
    arma::vec eigenenergies{0, 0.5, 1};
    arma::mat eigenstates{{0, 1, -0.5},
                          {-M_SQRT1_2, 0, 0.5},
                          {M_SQRT1_2, 0, M_SQRT1_2}};
    Eigensystem eigensystem(eigenenergies, eigenstates);

    REQUIRE_FALSE(eigensystem.empty());
    REQUIRE(eigensystem.size() == 3);
    REQUIRE(eigensystem.hasEigenvectors());
    REQUIRE_THAT(eigensystem.getEigenstates(), IsApproxEqual(eigenstates, 1e-16));
    REQUIRE_THAT(eigensystem.getEigenenergies(), IsApproxEqual(eigenenergies, 1e-15));
    REQUIRE_THAT(eigensystem.getEigenstate(0), IsApproxEqual(arma::vec{0, -M_SQRT1_2, M_SQRT1_2}, 1e-15));
    REQUIRE_THAT(eigensystem.getEigenstate(1), IsApproxEqual(arma::vec{1, 0, 0}, 1e-15));
    REQUIRE_THAT(eigensystem.getEigenstate(2), IsApproxEqual(arma::vec{-0.5, 0.5, M_SQRT1_2}, 1e-15));
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

TEST_CASE("Eigensystem: store/restore") {
    Eigensystem toStore({0, 0.25, 0.5, 1});
    Eigensystem toRestore({0, 0.1, 0.8, 1});

    std::stringstream inout;
    toStore.store(inout);
    toRestore.restore(inout);

    REQUIRE(toStore == toRestore);
}

TEST_CASE("Eigensystem: normalized energies") {
    Eigensystem eigensystem({1, 2, 3});

    REQUIRE_THAT(eigensystem.getNormalizedEigenenergies(), IsApproxEqual(arma::vec{0, 0.5, 1}, 1e-15));
}

TEST_CASE("Eigensystem: unmatching sizes") {
    REQUIRE_THROWS(Eigensystem({1, 2}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}));
    REQUIRE_THROWS(Eigensystem({1, 2}, {{1, 2}, {4, 5}, {7, 8}}));
}

TEST_CASE("Eigensystem: zero vector") {
    REQUIRE_THROWS(Eigensystem({0, 1}, {{0, 2}, {0, 1}}));
}

TEST_CASE("Eigensystem: cannot normalize equal eigenvalues") {
    REQUIRE_THROWS(Eigensystem({2, 2, 2}).getNormalizedEigenenergies());
}