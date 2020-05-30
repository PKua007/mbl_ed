//
// Created by Piotr Kubala on 20/02/2020.
//

#include <catch2/catch.hpp>
#include <iterator>

#include "matchers/ArmaApproxEqualMatcher.h"
#include "matchers/VectorApproxEqualMatcher.h"

#include "evolution/OccupationEvolution.h"
#include "evolution/EDEvolver.h"
#include "simulation/FockBaseGenerator.h"

TEST_CASE("OccupationEvolution: 1 boson 4 sites") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(1, 4));
    // Just a random nontrivial orthogonal matrix as eigenstates
    const double sq2 = std::sqrt(2);
    const double sq6 = std::sqrt(6);
    const double sq12 = std::sqrt(12);
    Eigensystem eigensystem({1, 2, 3, 4}, {{   0.5,    0.5,    0.5,     0.5},
                                           { 1/sq2, -1/sq2,      0,       0},
                                           { 1/sq6,  1/sq6, -2/sq6,       0},
                                           {1/sq12, 1/sq12, 1/sq12, -3/sq12}}, fockBase);
    std::ostringstream logger;
    EDEvolver evolver(eigensystem);


    // Initial state - {0, 1, 0, 0}
    auto evolution = OccupationEvolution::perform({{2, 1}}, 1, *fockBase, evolver, logger);


    REQUIRE(evolution.size() == 2);

    // 0 time - of course nothing should change from the initial state
    // only <n_1> != 0 (= 1), all <n_i n_j> = 0 apart from <n_1 n_1> = 1
    std::vector<double> nsExpected = {0, 1, 0, 0};
    REQUIRE_THAT(evolution[0].numParticles, IsApproxEqual(nsExpected, 1e-15));
    arma::vec nnsExpected = {0, 1, 0, 0};
    REQUIRE_THAT(evolution[0].numParticlesSquared.toArma(), IsApproxEqual(arma::mat(arma::diagmat(nnsExpected)), 1e-12));

    // The values were calculated in Mathematica by creating hamiltonian in Fock basis from eigenvectors and
    // eigenvalues, computing |psi(t)> = e^(-iHt) |psi> explicitly and calculating all <psi(t)|observable|psi(t)>
    nsExpected = {0.3540367091367856, 0.2919265817264288, 0.2360244727578571, 0.11801223637892853};
    REQUIRE_THAT(evolution[1].numParticles, IsApproxEqual(nsExpected, 1e-15));
    nnsExpected = arma::vec{0.3540367091367856, 0.2919265817264288, 0.2360244727578571, 0.11801223637892853};
    REQUIRE_THAT(evolution[1].numParticlesSquared.toArma(), IsApproxEqual(arma::mat(arma::diagmat(nnsExpected)), 1e-12));
}

TEST_CASE("OccupationEvolution: 2 bosons 2 sites") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(2, 2));
    Eigensystem eigensystem({1, 2, 3}, {{ 3, 2,  6},
                                        {-6, 3,  2},
                                        { 2, 6, -3}}, fockBase);
    EDEvolver evolver(eigensystem);
    std::ostringstream logger;


    // Initial vector - {1, 1}
    auto evolution = OccupationEvolution::perform({{2, 1}}, 1, *fockBase, evolver, logger);


    REQUIRE(evolution.size() == 2);

    // all <n_i> and <n_i n_j> = 1 for t = 0
    REQUIRE(evolution[0].numParticles.size() == 2);
    REQUIRE(evolution[0].numParticles[0] == Approx(1));
    REQUIRE(evolution[0].numParticles[1] == Approx(1));
    REQUIRE_THAT(evolution[0].numParticlesSquared.toArma(), IsApproxEqual(arma::mat{{1, 1}, {1, 1}}, 1e-12));

    REQUIRE(evolution[1].numParticles.size() == 2);
    REQUIRE(evolution[1].numParticles[0] == Approx(1.0569754884490989));
    REQUIRE(evolution[1].numParticles[1] == Approx(0.9430245115509011));
    REQUIRE_THAT(evolution[1].numParticlesSquared.toArma(),
                 IsApproxEqual(
                     arma::mat{{ 1.736972669993851, 0.3769783069043470},
                               {0.3769783069043470,  1.509070716197455}},
                     1e-12
                 ));
}