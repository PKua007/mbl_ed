//
// Created by Piotr Kubala on 19/02/2020.
//


#include <catch2/catch.hpp>
#include <chrono>
#include <catch2/trompeloeil.hpp>

#include "mocks/DiagonalTermMock.h"

#include "analyzer/tasks/EDCorrelationsTimeEvolution.h"
#include "evolution/EDEvolver.h"
#include "core/FockBaseGenerator.h"
#include "core/HamiltonianGenerator.h"
#include "core/terms/HubbardHop.h"

/*TEST_CASE("EDCorrelationsTimeEvolution: benchmark") {
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

    EDCorrelationsTimeEvolution correlationsTimeEvolution(1, 11, 1,
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
}*/

TEST_CASE("CorrelationsTimeEvolution: header") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(1, 5));
    Eigensystem eigensystem({1, 1, 1, 1, 1}, arma::eye(5, 5), fockBase);
    EDCorrelationsTimeEvolution evolution({{{2, 1}}, 5, 1, fockBase,
                                          {FockBase::Vector{1, 0, 0, 0, 0}, FockBase::Vector{0, 1, 0, 0, 0}}});
    std::ostringstream logger;

    evolution.analyze(eigensystem, logger);
    std::stringstream out;
    evolution.storeResult(out);

    std::string line;
    std::getline(out, line);
    REQUIRE(line == "1.0.0.0.0_t Gm0_1 Gm0_2 Gm0_3 Gm0_4 Gm1_1 Gm1_2 rho_0 rho_1 rho_2 rho_3 rho_4 n_0 n_1 n_2 n_3 n_4 "
                    "0.1.0.0.0_t Gm0_1 Gm0_2 Gm0_3 Gm0_4 Gm1_1 Gm1_2 rho_0 rho_1 rho_2 rho_3 rho_4 n_0 n_1 n_2 n_3 n_4 ");
}

TEST_CASE("CorrelationsTimeEvolution: times") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(1, 2));
    Eigensystem eigensystem({1, 1}, arma::eye(2, 2), fockBase);
    EDCorrelationsTimeEvolution evolution({{{1, 2}, {3, 1}}, 2, 0, fockBase, {FockBase::Vector{1, 0}}});
    std::ostringstream logger;

    evolution.analyze(eigensystem, logger);
    std::stringstream out;
    evolution.storeResult(out);

    std::string line;
    std::getline(out, line);
    double t, G_1, bG_1, rho_0, rho_1, n_0, n_1;
    out >> t >> G_1 >> bG_1 >> rho_0 >> rho_1 >> n_0 >> n_1;
    REQUIRE(t == 0);
    out >> t >> G_1 >> bG_1 >> rho_0 >> rho_1 >> n_0 >> n_1;
    REQUIRE(t == 0.5);
    out >> t >> G_1 >> bG_1 >> rho_0 >> rho_1 >> n_0 >> n_1;
    REQUIRE(t == 1);
    out >> t >> G_1 >> bG_1 >> rho_0 >> rho_1 >> n_0 >> n_1;
    REQUIRE(t == 3);
}

TEST_CASE("CorrelationsTimeEvolution: external vectors") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(1, 2));
    Eigensystem eigensystem({1, 1}, arma::eye(2, 2), fockBase);
    using ExternalVector = CorrelationsTimeEvolutionParameters::ExternalVector;
    CorrelationsTimeEvolution evolution({{{1, 1}}, 2, 0, fockBase,
                                         {FockBase::Vector{1, 0}, ExternalVector{"external"}}
                                        });
    std::ostringstream logger;
    EDEvolver evolver(eigensystem);

    SECTION("correct number of external vectors") {
        evolution.addEvolution(evolver, logger, {{M_SQRT1_2, M_SQRT1_2}});
        std::stringstream out;
        evolution.storeResult(out);

        std::string line;
        std::getline(out, line);
        // Check header
        REQUIRE(line == "1.0_t Gm0_1 Gm0_1 rho_0 rho_1 n_0 n_1 external_t Gm0_1 Gm0_1 rho_0 rho_1 n_0 n_1 ");
        // Make sure that initial values of <n_0> and <n_1> for external vector are correct
        double t, G_1, bG_1, rho_0, rho_1, Fock_n_0, Fock_n_1, external_n_0, external_n_1;
        out >> t >> G_1 >> bG_1 >> rho_0 >> rho_1 >> Fock_n_0 >> Fock_n_1;
        out >> t >> G_1 >> bG_1 >> rho_0 >> rho_1 >> external_n_0 >> external_n_1;
        REQUIRE(external_n_0 == Approx(0.5));
        REQUIRE(external_n_1 == Approx(0.5));
    }

    SECTION("throw on incorrect number of external vectors") {
        REQUIRE_THROWS(evolution.addEvolution(evolver, logger, {}));
        REQUIRE_THROWS(evolution.addEvolution(evolver, logger, {{1, 0}, {0, 1}}));
    }
}

TEST_CASE("CorrelationsTimeEvolution: throw on non-matching eigensystems") {
    auto fockBase1 = std::shared_ptr(FockBaseGenerator{}.generate(1, 2));
    Eigensystem eigensystem1({1, 1}, arma::eye(2, 2), fockBase1);
    auto fockBase2 = FockBaseGenerator{}.generate(1, 3);
    Eigensystem eigensystem2({1, 1, 1}, arma::eye(3, 3), std::move(fockBase2));
    EDCorrelationsTimeEvolution evolution({{{2, 1}}, 2, 0, fockBase1, {FockBase::Vector{1, 0}}});
    std::ostringstream logger;
    evolution.analyze(eigensystem1, logger);

    REQUIRE_THROWS(evolution.analyze(eigensystem2, logger));
}

TEST_CASE("CorrelationsTimeEvolution: throw on wrong initial vectors") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(1, 2));
    Eigensystem eigensystem({1, 1}, arma::eye(2, 2), fockBase);
    EDCorrelationsTimeEvolution corr1({{{2, 1}}, 2, 0, fockBase, {FockBase::Vector{1, 0, 0}}});
    EDCorrelationsTimeEvolution corr2({{{2, 1}}, 2, 0, fockBase, {FockBase::Vector{0, 2}}});
    std::ostringstream logger;

    REQUIRE_THROWS(corr1.analyze(eigensystem, logger));
    REQUIRE_THROWS(corr2.analyze(eigensystem, logger));
}