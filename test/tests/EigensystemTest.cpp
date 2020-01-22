//
// Created by Piotr Kubala on 22/01/2020.
//

#include <catch2/catch.hpp>

#include "simulation/Eigensystem.h"

TEST_CASE("Eigensystem: empty") {
    Eigensystem eigensystem;

    REQUIRE(eigensystem.empty() == true);
    REQUIRE(eigensystem.size() == 0);
    REQUIRE(eigensystem.getEigenenergies() == std::vector<double>{});
    REQUIRE(eigensystem.getEigenstates() == std::vector<std::vector<double>>{});
    REQUIRE(eigensystem.hasEigenvectors() == false);
    REQUIRE(eigensystem.isComplete() == true);
}

TEST_CASE("Eigensystem: no eigenvectors") {
    Eigensystem eigensystem;
    eigensystem.addEntry(0, {});
    eigensystem.addEntry(1, {});
    eigensystem.addEntry(2, {});

    REQUIRE(eigensystem.empty() == false);
    REQUIRE(eigensystem.size() == 3);
    REQUIRE(eigensystem.getEigenenergies() == std::vector<double>{0, 1, 2});
    REQUIRE(eigensystem.getEigenstates() == std::vector<std::vector<double>>{{}, {}, {}});
    REQUIRE(eigensystem.hasEigenvectors() == false);
    REQUIRE(eigensystem.isComplete() == true);
}

TEST_CASE("Eigensystem: has eigenvectors") {
    Eigensystem eigensystem;
    eigensystem.addEntry(0, {1, 2, 3});
    eigensystem.addEntry(1, {4, 5, 6});
    eigensystem.addEntry(2, {7, 8, 9});

    REQUIRE(eigensystem.empty() == false);
    REQUIRE(eigensystem.size() == 3);
    REQUIRE(eigensystem.getEigenenergies() == std::vector<double>{0, 1, 2});
    REQUIRE(eigensystem.getEigenstates() == std::vector<std::vector<double>>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    REQUIRE(eigensystem.hasEigenvectors() == true);
    REQUIRE(eigensystem.isComplete() == true);
}

TEST_CASE("Eigensystem: incomplete") {
    Eigensystem eigensystem;
    eigensystem.addEntry(0, {1, 2, 3});
    eigensystem.addEntry(1, {4, 5, 6});

    REQUIRE(eigensystem.empty() == false);
    REQUIRE(eigensystem.size() == 2);
    REQUIRE(eigensystem.getEigenenergies() == std::vector<double>{0, 1});
    REQUIRE(eigensystem.getEigenstates() == std::vector<std::vector<double>>{{1, 2, 3}, {4, 5, 6}});
    REQUIRE(eigensystem.hasEigenvectors() == true);
    REQUIRE(eigensystem.isComplete() == false);
}

TEST_CASE("Eigensystem: throw if too many entries") {
    Eigensystem eigensystem;
    eigensystem.addEntry(0, {1, 2, 3});
    eigensystem.addEntry(1, {4, 5, 6});
    eigensystem.addEntry(2, {7, 8, 9});

    REQUIRE_THROWS(eigensystem.addEntry(3, {10, 11, 12}));
}

TEST_CASE("Eigensystem: throw if non-matching vector sizes") {
    Eigensystem eigensystem1;
    eigensystem1.addEntry(0, {1, 2, 3});
    Eigensystem eigensystem2;
    eigensystem2.addEntry(1, {});

    REQUIRE_THROWS(eigensystem1.addEntry(2, {1, 2}));
    REQUIRE_THROWS(eigensystem2.addEntry(2, {1, 2}));
}