//
// Created by Piotr Kubala on 19/02/2020.
//

#include <vector>

#include <catch2/catch.hpp>

#include "mocks/ObservableMock.h"
#include "mocks/OccupationEvolutionMock.h"

#include "evolution/TimeEvolution.h"
#include "evolution/EDEvolver.h"
#include "core/FockBaseGenerator.h"

using trompeloeil::_;

TEST_CASE("CorrelationsTimeEvolution: header and entries") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(1, 2));
    Eigensystem eigensystem({1, 1}, arma::eye(2, 2), fockBase);
    TimeEvolutionParameters params;
    params.fockBase = fockBase;
    params.timeSegmentation = {{1, 2}, {3, 1}};
    params.numberOfSites = 2;
    params.vectorsToEvolve = {FockBase::Vector{1, 0}, FockBase::Vector{0, 1}};

    auto observable1 = std::make_shared<ObservableMock>();
    auto observable2 = std::make_shared<ObservableMock>();
    ALLOW_CALL(*observable1, getHeader()).RETURN(std::vector<std::string>{"o1.1", "o1.2"});
    ALLOW_CALL(*observable2, getHeader()).RETURN(std::vector<std::string>{"o2.1", "o2.2"});
    params.storedObservables = {observable1, observable2};

    auto occupationEvolution = std::make_unique<OccupationEvolutionMock>();
    REQUIRE_CALL(*occupationEvolution, perform(params.timeSegmentation, _, _, _))
        .WITH(arma::approx_equal(_2, arma::cx_vec{1, 0}, "absdiff", 1e-8))
        .RETURN(std::vector<TimeEvolutionEntry>{
            {0,   {1, 2, 3, 4}},
            {0.5, {5, 6, 7, 8}},
            {1,   {9, 10, 11, 12}},
            {3,   {13, 14, 15, 16}},
        });
    REQUIRE_CALL(*occupationEvolution, perform(params.timeSegmentation, _, _, _))
        .WITH(arma::approx_equal(_2, arma::cx_vec{0, 1}, "absdiff", 1e-8))
        .RETURN(std::vector<TimeEvolutionEntry>{
            {0,   {17, 18, 19, 20}},
            {0.5, {21, 22, 23, 24}},
            {1,   {25, 26, 27, 28}},
            {3,   {29, 30, 31, 32}},
        });

    EDEvolver evolver(eigensystem);
    TimeEvolution evolution(params, std::move(occupationEvolution));
    std::ostringstream loggerStream;
    Logger logger(loggerStream);


    evolution.addEvolution(evolver, logger);
    std::stringstream out;
    evolution.storeResult(out);

    SECTION("header") {
        std::string line;
        std::getline(out, line);
        REQUIRE(line == "1.0_t o1.1 o1.2 o2.1 o2.2 0.1_t o1.1 o1.2 o2.1 o2.2 ");

        SECTION("entries") {
            std::getline(out, line);
            REQUIRE(line == "0 1 2 3 4 0 17 18 19 20 ");
            std::getline(out, line);
            REQUIRE(line == "0.5 5 6 7 8 0.5 21 22 23 24 ");
            std::getline(out, line);
            REQUIRE(line == "1 9 10 11 12 1 25 26 27 28 ");
            std::getline(out, line);
            REQUIRE(line == "3 13 14 15 16 3 29 30 31 32 ");
        }
    }
}

TEST_CASE("CorrelationsTimeEvolution: external vectors") {
    auto fockBase = std::shared_ptr(FockBaseGenerator{}.generate(1, 2));
    Eigensystem eigensystem({1, 1}, arma::eye(2, 2), fockBase);
    TimeEvolutionParameters params;
    params.fockBase = fockBase;
    params.timeSegmentation = {{1, 1}};
    params.numberOfSites = 2;
    params.vectorsToEvolve = {FockBase::Vector{1, 0}, TimeEvolutionParameters::ExternalVector{"external"}};

    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    SECTION("correct number of external vectors") {
        auto occupationEvolution = std::make_unique<OccupationEvolutionMock>();
        REQUIRE_CALL(*occupationEvolution, perform(params.timeSegmentation, _, _, _))
                .WITH(arma::approx_equal(_2, arma::cx_vec{1, 0}, "absdiff", 1e-8))
                .RETURN(std::vector<TimeEvolutionEntry>{{0, {}}, {1, {}}});
        REQUIRE_CALL(*occupationEvolution, perform(params.timeSegmentation, _, _, _))
                .WITH(arma::approx_equal(_2, arma::cx_vec{M_SQRT1_2, M_SQRT1_2}, "absdiff", 1e-8))
                .RETURN(std::vector<TimeEvolutionEntry>{{0, {}}, {1, {}}});

        EDEvolver evolver(eigensystem);
        TimeEvolution evolution(params, std::move(occupationEvolution));

        evolution.addEvolution(evolver, logger, {{M_SQRT1_2, M_SQRT1_2}});
        std::stringstream out;
        evolution.storeResult(out);

        std::string line;
        std::getline(out, line);
        // Check header
        REQUIRE(line == "1.0_t external_t ");
    }

    SECTION("throw on incorrect number of external vectors") {
        auto occupationEvolution = std::make_unique<OccupationEvolutionMock>();
        EDEvolver evolver(eigensystem);
        TimeEvolution evolution(params, std::move(occupationEvolution));

        REQUIRE_THROWS(evolution.addEvolution(evolver, logger, {}));
        REQUIRE_THROWS(evolution.addEvolution(evolver, logger, {{1, 0}, {0, 1}}));
    }
}