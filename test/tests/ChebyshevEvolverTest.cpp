//
// Created by pkua on 02.05.2020.
//

#include <cmath>
#include <catch2/catch.hpp>

#include "simulation/FockBaseGenerator.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/terms/HubbardHop.h"
#include "simulation/terms/HubbardOnsite.h"

#include "evolution/EDEvolver.h"
#include "evolution/ChebyshevEvolver.h"


TEST_CASE("ChebycshevEvolver test") {
    FockBaseGenerator generator;
    auto basis = std::shared_ptr(generator.generate(3, 3));
    HamiltonianGenerator hamiltonianGenerator(basis, false);
    hamiltonianGenerator.addHoppingTerm(std::make_unique<HubbardHop>(1));
    hamiltonianGenerator.addDiagonalTerm(std::make_unique<HubbardOnsite>(2));
    auto H = hamiltonianGenerator.generate();

    arma::cx_vec psi0(basis->size(), arma::fill::zeros);
    psi0[0] = 1;
    EDEvolver edEvolver;
    edEvolver.prepareFor(H, psi0, 2, 2);
    edEvolver.evolve();
    ChebyshevEvolver chebyshevEvolver;
    chebyshevEvolver.prepareFor(H, psi0, 2, 2);
    chebyshevEvolver.evolve();

    REQUIRE(arma::norm(edEvolver.getCurrentState() - chebyshevEvolver.getCurrentState()) < 1e-8);
}