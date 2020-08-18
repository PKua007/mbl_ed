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

    // The values were calculated in Mathematica by creating hamiltonian in Fock basis from eigenvectors and
    // eigenvalues, computing |psi(t)> = e^(-iHt) |psi> explicitly and calculating all <psi(t)|observable|psi(t)>
    std::vector<double> nsT0 = {0, 1, 0, 0};
    arma::vec nnsT0 = {0, 1, 0, 0};
    std::vector<double> nsT2 = {0.3540367091367856, 0.2919265817264288, 0.2360244727578571, 0.11801223637892853};
    arma::vec nnsT2 = arma::vec{0.3540367091367856, 0.2919265817264288, 0.2360244727578571, 0.11801223637892853};

    SECTION("single time segment - t = 0, 2") {
        OccupationEvolution occupationEvolution(fockBase);
        auto evolution = occupationEvolution.perform({{2, 1}}, {0, 1, 0, 0}, evolver, logger);

        REQUIRE(evolution.size() == 2);
        REQUIRE_THAT(evolution[0].numParticles, IsApproxEqual(nsT0, 1e-15));
        REQUIRE_THAT(evolution[0].numParticlesSquared.toArma(),IsApproxEqual(arma::mat(arma::diagmat(nnsT0)), 1e-15));
        REQUIRE_THAT(evolution[1].numParticles, IsApproxEqual(nsT2, 1e-12));
        REQUIRE_THAT(evolution[1].numParticlesSquared.toArma(),IsApproxEqual(arma::mat(arma::diagmat(nnsT2)), 1e-12));
    }

    SECTION("2 time segments - t = 0, 0.5, 1, 2") {
        OccupationEvolution occupationEvolution(fockBase);
        auto evolution = occupationEvolution.perform({{1, 2}, {2, 1}}, {0, 1, 0, 0}, evolver, logger);

        REQUIRE(evolution.size() == 4);
        REQUIRE_THAT(evolution[0].numParticles, IsApproxEqual(nsT0, 1e-15));
        REQUIRE_THAT(evolution[0].numParticlesSquared.toArma(),IsApproxEqual(arma::mat(arma::diagmat(nnsT0)), 1e-15));
        REQUIRE_THAT(evolution[3].numParticles, IsApproxEqual(nsT2, 1e-12));
        REQUIRE_THAT(evolution[3].numParticlesSquared.toArma(),IsApproxEqual(arma::mat(arma::diagmat(nnsT2)), 1e-12));
    }
}

TEST_CASE("OccupationEvolution: 2 bosons 2 sites") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(2, 2));
    Eigensystem eigensystem({1, 2, 3}, {{ 3, 2,  6},
                                        {-6, 3,  2},
                                        { 2, 6, -3}}, fockBase);
    EDEvolver evolver(eigensystem);
    std::ostringstream logger;


    OccupationEvolution occupationEvolution(fockBase);
    auto evolution = occupationEvolution.perform({{2, 1}}, {0, 1, 0}, evolver, logger);


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