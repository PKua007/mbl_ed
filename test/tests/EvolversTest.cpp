//
// Created by pkua on 02.05.2020.
//

#include <cmath>
#include <complex>

#include <catch2/catch.hpp>

#include "simulation/FockBaseGenerator.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/terms/HubbardHop.h"
#include "simulation/terms/HubbardOnsite.h"

#include "evolution/EDEvolver.h"
#include "evolution/ChebyshevEvolver.h"

using namespace std::complex_literals;

TEST_CASE("Evolvers test") {
    FockBaseGenerator generator;
    auto basis = std::shared_ptr(generator.generate(2, 3));
    HamiltonianGenerator hamiltonianGenerator(basis, false);
    hamiltonianGenerator.addHoppingTerm(std::make_unique<HubbardHop>(1));
    hamiltonianGenerator.addDiagonalTerm(std::make_unique<HubbardOnsite>(2));
    auto H = hamiltonianGenerator.generate();
    std::ostringstream logger;

    arma::cx_vec psi0(basis->size(), arma::fill::zeros);
    psi0[0] = 1;

    arma::cx_vec expected(6);
    expected[0] = 0.09916801545681185 + 0.488086379991221i;
    expected[1] = -0.1389521483608867 + 0.05922080722864662i;
    expected[2] = -0.09154129939615813 + -0.02895753433723917i;
    expected[3] = 0.2908961074241866 + -0.5668961787668898i;
    expected[4] = 0.09636154567496999 + -0.04847228249128829i;
    expected[5] = -0.4619155288962375 + -0.2981800634501822i;

    SECTION("EDEvolver") {
        arma::vec eigval;
        arma::mat eigvec;
        REQUIRE(arma::eig_sym(eigval, eigvec, arma::mat(H)));
        Eigensystem eigensystem(eigval, eigvec);
        EDEvolver edEvolver(eigensystem);
        edEvolver.prepareFor(psi0, 2, 2);
        edEvolver.evolve();

        REQUIRE(arma::norm(edEvolver.getCurrentState() - expected) < 1e-11);

        SECTION("throw on to many evolutions") {
            edEvolver.evolve(); // 1 more step is permitted for convenience
            REQUIRE_THROWS(edEvolver.evolve());
        }

        SECTION("second run - reset") {
            // Different evolution - different number of steps - but time after first step coincides with previous one
            // However not to use exactly the same, we take -psi0
            edEvolver.prepareFor(-psi0, 4, 3);
            edEvolver.evolve();

            REQUIRE(arma::norm(edEvolver.getCurrentState() - (-expected)) < 1e-11);
        }
    }

    SECTION("ChebyshevEvolver") {
        ChebyshevEvolver chebyshevEvolver(H, logger, 2);
        chebyshevEvolver.prepareFor(psi0, 2, 2);
        chebyshevEvolver.evolve();

        REQUIRE(arma::norm(chebyshevEvolver.getCurrentState() - expected) < 1e-11);

        SECTION("throw on to many evolutions") {
            chebyshevEvolver.evolve(); // 1 more step is permitted for convenience
            REQUIRE_THROWS(chebyshevEvolver.evolve());
        }

        SECTION("second run - reset") {
            // Different evolution - different number of steps - but time after first step coincides with previous one
            // However not to use exactly the same, we take -psi0
            chebyshevEvolver.prepareFor(-psi0, 4, 3);
            chebyshevEvolver.evolve();

            REQUIRE(arma::norm(chebyshevEvolver.getCurrentState() - (-expected)) < 1e-11);
        }
    }

    SECTION("ChebyshevEvolver - long times with rebuilds") {
        // We compare EDEvolver result ...
        arma::vec eigval;
        arma::mat eigvec;
        REQUIRE(arma::eig_sym(eigval, eigvec, arma::mat(H)));
        Eigensystem eigensystem(eigval, eigvec);
        EDEvolver edEvolver(eigensystem);
        edEvolver.prepareFor(psi0, 100, 2);
        edEvolver.evolve();
        // ... with out ChebyshevEvolver
        ChebyshevEvolver chebyshevEvolver(H, logger, 50);
        chebyshevEvolver.prepareFor(psi0, 100, 101);

        for (std::size_t i{}; i < 100; i++)
            chebyshevEvolver.evolve();

        REQUIRE(arma::norm(chebyshevEvolver.getCurrentState() - edEvolver.getCurrentState()) < 1e-8);
    }
}