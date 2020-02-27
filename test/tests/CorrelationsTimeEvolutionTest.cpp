//
// Created by Piotr Kubala on 19/02/2020.
//


#include <catch2/catch.hpp>

#include "analyzer/tasks/correlations_time_evolution/CorrelationsTimeEvolution.h"
#include "simulation/FockBaseGenerator.h"

/*

#include <chrono>
#include <catch2/trompeloeil.hpp>

#include "simulation/HamiltonianGenerator.h"
#include "simulation/terms/HubbardHop.h"
#include "mocks/DiagonalTermMock.h"

TEST_CASE("CorrelationsTimeEvolution: benchmark") {
    auto base = FockBaseGenerator{}.generate(6, 6);
    HamiltonianGenerator hamiltonianGenerator(std::move(base), false);
    auto diagonal = std::make_unique<DiagonalTermMock>();
    using trompeloeil::_;
    ALLOW_CALL(*diagonal, calculate(_, _))
        .RETURN(_2.getFockBase()->findIndex(_1).value());
    hamiltonianGenerator.addDiagonalTerm(std::move(diagonal));
    hamiltonianGenerator.addHoppingTerm(std::make_unique<HubbardHop>(1));
    arma::mat hamiltonian = hamiltonianGenerator.generate();
    arma::mat eigenvectors;
    arma::vec eigenvalues;
    arma::eig_sym(eigenvalues, eigenvectors, hamiltonian);
    Eigensystem eigensystem(eigenvalues, eigenvectors, hamiltonianGenerator.getFockBase());

    CorrelationsTimeEvolution correlationsTimeEvolution(0, 1, 11, CorrelationsTimeEvolution::Linear, 1,
                                                        {{1, 1, 1, 1, 1, 1, 1, 1}, {2, 0, 2, 0, 2, 0, 2, 0}});

    // Warmup
    for (std::size_t i = 0; i < 2; i++)
        correlationsTimeEvolution.analyze(eigensystem);

    // Actual test
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < 5; i++)
        correlationsTimeEvolution.analyze(eigensystem);
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
}

*/

TEST_CASE("CorrelationsTimeEvolution: header") {
    auto fockBase = FockBaseGenerator{}.generate(1, 5);
    Eigensystem eigensystem({1, 1, 1, 1, 1}, arma::eye(5, 5), std::move(fockBase));
    CorrelationsTimeEvolution evolution(2, 2, 1, {{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}});

    evolution.analyze(eigensystem);
    std::stringstream out;
    evolution.storeResult(out);

    std::string line;
    std::getline(out, line);
    REQUIRE(line == "1.0.0.0.0_t G_1 G_2 G_3 G_4 G_1 G_2 rho_0 rho_1 rho_2 rho_3 rho_4 "
                    "0.1.0.0.0_t G_1 G_2 G_3 G_4 G_1 G_2 rho_0 rho_1 rho_2 rho_3 rho_4 ");
}

TEST_CASE("CorrelationsTimeEvolution: times") {
    auto fockBase = FockBaseGenerator{}.generate(1, 2);
    Eigensystem eigensystem({1, 1}, arma::eye(2, 2), std::move(fockBase));
    CorrelationsTimeEvolution evolution(0.5, 3, 0, {{1, 0}});

    evolution.analyze(eigensystem);
    std::stringstream out;
    evolution.storeResult(out);

    std::string line;
    std::getline(out, line);
    double t, G_1, bG_1, rho_1, rho2;
    out >> t >> G_1 >> bG_1 >> rho_1 >> rho2;
    REQUIRE(t == 0);
    out >> t >> G_1 >> bG_1 >> rho_1 >> rho2;
    REQUIRE(t == 0.25);
    out >> t >> G_1 >> bG_1 >> rho_1 >> rho2;
    REQUIRE(t == 0.50);
}

TEST_CASE("CorrelationsTimeEvolution: throw on non-matching eigensystems") {
    auto fockBase1 = FockBaseGenerator{}.generate(1, 2);
    Eigensystem eigensystem1({1, 1}, arma::eye(2, 2), std::move(fockBase1));
    auto fockBase2 = FockBaseGenerator{}.generate(1, 3);
    Eigensystem eigensystem2({1, 1, 1}, arma::eye(3, 3), std::move(fockBase1));
    CorrelationsTimeEvolution evolution(2, 2, 0, {{1, 0}});
    evolution.analyze(eigensystem1);

    REQUIRE_THROWS(evolution.analyze(eigensystem2));
}

TEST_CASE("CorrelationsTimeEvolution: throw on wrong initial vectors") {
    auto fockBase = FockBaseGenerator{}.generate(1, 2);
    Eigensystem eigensystem({1, 1}, arma::eye(2, 2), std::move(fockBase));
    CorrelationsTimeEvolution corr1(2, 2, 0, {{1, 0, 0}});
    CorrelationsTimeEvolution corr2(2, 2, 0, {{0, 2}});

    REQUIRE_THROWS(corr1.analyze(eigensystem));
    REQUIRE_THROWS(corr2.analyze(eigensystem));
}