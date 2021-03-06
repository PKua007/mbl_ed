//
// Created by Piotr Kubala on 18/08/2020.
//

#include <catch2/catch.hpp>

#include "matchers/ArmaApproxEqualCatchMatcher.h"

#include "core/QuenchCalculator.h"

TEST_CASE("QuenchCalculator: calculation") {
    arma::sp_mat finalHamiltonian(2, 2);
    finalHamiltonian(0, 0) = 1;
    finalHamiltonian(1, 1) = 2;

    arma::sp_mat initialHamiltonian1(2, 2);  // ground state {1/sqrt(2), -1/sqrt(2)}
    initialHamiltonian1(0, 0) = 1; initialHamiltonian1(0, 1) = 1;
    initialHamiltonian1(1, 0) = 1; initialHamiltonian1(1, 1) = 1;

    arma::sp_mat initialHamiltonian2(2, 2);  // ground state {-1-sqrt(2), 1} * normalization
    initialHamiltonian2(0, 0) = -1; initialHamiltonian2(0, 1) = 1;
    initialHamiltonian2(1, 0) =  1; initialHamiltonian2(1, 1) = 1;

    SECTION("single quench") {
        QuenchCalculator quenchCalculator;
        quenchCalculator.addQuench(initialHamiltonian1, finalHamiltonian);

        arma::vec expectedState = arma::vec{-M_SQRT1_2, M_SQRT1_2};
        arma::vec minusExpectedState = -expectedState;
        REQUIRE_THAT(quenchCalculator.getLastQuenchedState(),
                     IsApproxEqual(expectedState, 1e-8) || IsApproxEqual(minusExpectedState, 1e-8));
        REQUIRE(quenchCalculator.getLastQuenchEpsilon() == Approx(0.5));
        REQUIRE(quenchCalculator.getMeanEpsilonQuantumUncertainty() == Approx(0.5));

        SECTION("clearing") {
            quenchCalculator.clear();

            REQUIRE(quenchCalculator.getMeanEpsilon() == 0);
            REQUIRE(quenchCalculator.getMeanEpsilonQuantumUncertainty() == 0);
            REQUIRE(quenchCalculator.getEpsilonAveragingSampleError() == 0);
            REQUIRE(quenchCalculator.getLastQuenchedState().empty());
        }
    }

    SECTION("two quenches") {
        QuenchCalculator quenchCalculator;

        SECTION("averages - correct values") {
            quenchCalculator.addQuench(initialHamiltonian1, finalHamiltonian);
            quenchCalculator.addQuench(initialHamiltonian2, finalHamiltonian);

            REQUIRE(quenchCalculator.getMeanEpsilon() == Approx(0.3232233047033631));
            REQUIRE(quenchCalculator.getMeanEpsilonQuantumUncertainty() == Approx(0.4330127018922193));
            REQUIRE(quenchCalculator.getEpsilonAveragingSampleError() == Approx(0.25));
        }

        SECTION("storing and restoring with join") {
            QuenchCalculator secondQuenchCalculator;
            secondQuenchCalculator.addQuench(initialHamiltonian2, finalHamiltonian);
            std::stringstream memory;

            secondQuenchCalculator.storeState(memory);
            quenchCalculator.addQuench(initialHamiltonian1, finalHamiltonian);
            quenchCalculator.joinRestoredState(memory);

            REQUIRE(quenchCalculator.getMeanEpsilon() == Approx(0.3232233047033631));
            REQUIRE(quenchCalculator.getMeanEpsilonQuantumUncertainty() == Approx(0.4330127018922193));
            REQUIRE(quenchCalculator.getEpsilonAveragingSampleError() == Approx(0.25));
        }


    }
}