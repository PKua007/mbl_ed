//
// Created by pkua on 02.05.2020.
//

#include <cmath>

#include <catch2/catch.hpp>

#include "simulation/FockBaseGenerator.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/terms/HubbardHop.h"
#include "simulation/terms/HubbardOnsite.h"

TEST_CASE("ChebycshevEvolver test") {
    FockBaseGenerator generator;
    auto basis = generator.generate(3, 3);
    HamiltonianGenerator hamiltonianGenerator(std::move(basis), false);
    hamiltonianGenerator.addHoppingTerm(std::make_unique<HubbardHop>(1));
    hamiltonianGenerator.addDiagonalTerm(std::make_unique<HubbardOnsite>(2));

    auto H = hamiltonianGenerator.generate();

    std::cout << arma::mat(H) << std::endl;
}