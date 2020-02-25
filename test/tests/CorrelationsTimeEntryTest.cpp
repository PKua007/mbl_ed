//
// Created by Piotr Kubala on 21/02/2020.
//

#include <sstream>

#include <catch2/catch.hpp>

#include "analyzer/tasks/correlations_time_evolution/CorrelationsTimeEntry.h"

TEST_CASE("CorrelationsTimeEntry: basic") {
    SECTION("empty") {
        CorrelationsTimeEntry correlationsTimeEntry{};

        REQUIRE(correlationsTimeEntry.getNumberOfSites() == 0);
    }

    SECTION("non-empty") {
        CorrelationsTimeEntry correlationsTimeEntry(2, 1, 5);

        REQUIRE(correlationsTimeEntry.getNumberOfSites() == 5);
        REQUIRE(correlationsTimeEntry.getHeader() == "t x G_1 G_2 G_3 G_4 G_1 G_2 rho_0 rho_1 rho_2 rho_3 rho_4 ");
        REQUIRE_THAT(correlationsTimeEntry.toString(), Catch::StartsWith("2 "));
    }
}

TEST_CASE("CorrelationsTimeEntry: single observables set") {
    OccupationEvolution::Occupations o(5);
    o.numParticles = {1, 2, 3, 4, 5};
    o.numParticlesSquared(0, 0) = 6; o.numParticlesSquared(0, 1) = 7; o.numParticlesSquared(0, 2) = 8; o.numParticlesSquared(0, 3) = 9; o.numParticlesSquared(0, 4) = 10;
    o.numParticlesSquared(1, 1) = 11; o.numParticlesSquared(1, 2) = 12; o.numParticlesSquared(1, 3) = 13; o.numParticlesSquared(1, 4) = 14;
    o.numParticlesSquared(2, 2) = 15; o.numParticlesSquared(2, 3) = 16; o.numParticlesSquared(2, 4) = 17;
    o.numParticlesSquared(3, 3) = 18; o.numParticlesSquared(3, 4) = 19;
    o.numParticlesSquared(4, 4) = 20;
    CorrelationsTimeEntry correlationsTimeEntry(2, 1, 5);

    correlationsTimeEntry.addObservables(o);
    std::istringstream out(correlationsTimeEntry.toString());

    double t, x, G_1, G_2, G_3, G_4, bG_1, bG_2, rho_0, rho_1, rho_2, rho_3, rho_4;
    out >> t >> x >> G_1 >> G_2 >> G_3 >> G_4 >> bG_1 >> bG_2 >> rho_0 >> rho_1 >> rho_2 >> rho_3 >> rho_4;
    REQUIRE(t == 2);
    REQUIRE(x == 90);           // 2*(1*3.5 + 2*4 + 3*4.5 + 4*5) = 90
    REQUIRE(G_1 == 3.5);        // ((7 - 1*2) + (12 - 2*3) + (16 - 3*4) + (19 - 4*5)) / 4 = (5 + 6 + 4 - 1) / 4 = 3.5
    REQUIRE(G_2 == 4);          // ((8 - 1*3) + (13 - 2*4) + (17 - 3*5))/3 = (5 + 5 + 2) / 3 = 4
    REQUIRE(G_3 == 4.5);        // ((9 - 1*4) + (14 - 2*5)) / 2 = (5 + 4) / 2 = 4.5
    REQUIRE(G_4 == 5);
    REQUIRE(bG_1 == 5);         // ((12 - 2*3) + (16 - 3*4)) / 2 = (6 + 4) / 2 = 5
    REQUIRE(bG_2 == 5);
    REQUIRE(rho_0 == 5);
    REQUIRE(rho_1 == 7);
    REQUIRE(rho_2 == 6);
    REQUIRE(rho_3 == 2);
    REQUIRE(rho_4 == -5);
}

TEST_CASE("CorrelationsTimeEntry: averaging") {
    OccupationEvolution::Occupations o1(3), o2(3);
    o1.numParticles = {1, 2, 3};
    o1.numParticlesSquared(0, 0) = 4; o1.numParticlesSquared(0, 1) = 5; o1.numParticlesSquared(0, 2) = 6;
    o1.numParticlesSquared(1, 1) = 7; o1.numParticlesSquared(1, 2) = 8;
    o1.numParticlesSquared(2, 2) = 9;
    o2.numParticles = {10, 11, 12};
    o2.numParticlesSquared(0, 0) = 13; o2.numParticlesSquared(0, 1) = 14; o2.numParticlesSquared(0, 2) = 15;
    o2.numParticlesSquared(1, 1) = 16; o2.numParticlesSquared(1, 2) = 17;
    o2.numParticlesSquared(2, 2) = 18;
    CorrelationsTimeEntry correlationsTimeEntry(2, 0, 3);

    correlationsTimeEntry.addObservables(o1);
    correlationsTimeEntry.addObservables(o2);
    std::istringstream out(correlationsTimeEntry.toString());

    double t, x, G_1, G_2, bG_1, bG_2, rho_0, rho_1, rho_2;
    out >> t >> x >> G_1 >> G_2 >> bG_1 >> bG_2 >> rho_0 >> rho_1 >> rho_2;
    REQUIRE(t == 2);
    REQUIRE(x == 307);          // 1 -> 2*|1*2.5 + 2*3| = 17; 2 -> 2*|1*(-105.5) + 2*(-105)| = 631; avg = 307
    REQUIRE(G_1 == -51.5);      // 1->((5-1*2)+(8-2*3))/2=(3+2)/2=2.5; 2->((14-10*11)+(17-11*12))/2=-105.5; avg=-51.5
    REQUIRE(G_2 == -51);        // 1 -> (6 - 1*3) = 3; 2 -> (15 - 10*12) -> -105; avg = -51
    REQUIRE(bG_1 == -51.5);
    REQUIRE(bG_2 == -51);
    REQUIRE(rho_0 == -42);      // 1 -> 3; 2 -> -87; avg = -42
    REQUIRE(rho_1 == -51);      // 1 -> 3; 2 -> -105; avg = -51
    REQUIRE(rho_2 == -63);      // 1 -> (9 - 3*3) = 0; 2 -> (18 - 12*12) = ; avg = -63
}