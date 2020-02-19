//
// Created by Piotr Kubala on 19/02/2020.
//

#include <catch2/catch.hpp>
#include <catch2/trompeloeil.hpp>

#include "analyzer/tasks/CorrelationsTimeEvolution.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/FockBaseGenerator.h"
#include "simulation/terms/HubbardHop.h"
#include "mocks/DiagonalTermMock.h"

TEST_CASE("CorrelationsTimeEvolution") {
    auto base = FockBaseGenerator{}.generate(4, 4);
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

    CorrelationsTimeEvolution correlationsTimeEvolution(11, 20, 10, CorrelationsTimeEvolution::Linear, 1,
                                                        {{1, 1, 1, 1}, {2, 0, 2, 0}});
    correlationsTimeEvolution.analyze(eigensystem);

    std::ostringstream out;
    correlationsTimeEvolution.storeResult(out);
    std::cout << out.str();
}