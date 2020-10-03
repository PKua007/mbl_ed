//
// Created by Piotr Kubala on 20/02/2020.
//

#include <iterator>

#include <catch2/catch.hpp>
#include <catch2/trompeloeil.hpp>

#include "matchers/ArmaApproxEqualTrompeloeilMatcher.h"
#include "matchers/FirstObservableValuesTrompeloeilMatcher.h"

#include "mocks/PrimaryObservableMock.h"
#include "mocks/SecondaryObservableMock.h"
#include "mocks/EvolverMock.h"

#include "evolution/OservablesTimeEvolution.h"

using trompeloeil::_;

TEST_CASE("OservablesTimeEvolution") {
    trompeloeil::sequence seq1;
    auto primary = std::make_shared<PrimaryObservableMock>();
    ALLOW_CALL(*primary, getHeader()).RETURN(std::vector<std::string>{"p1", "p2"});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{1, 0, 0, 0, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{1, 2});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{0, 1, 0, 0, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{3, 4});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{0, 0, 1, 0, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{5, 6});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{0, 0, 0, 1, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{7, 8});

    trompeloeil::sequence seq3;
    auto secondary = std::make_shared<SecondaryObservableMock>();
    ALLOW_CALL(*secondary, getHeader()).RETURN(std::vector<std::string>{"s1", "s2"});
    REQUIRE_CALL(*secondary, calculateForObservables(first_observable_values_eq({1, 2}))).IN_SEQUENCE(seq3);
    REQUIRE_CALL(*secondary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq3).RETURN(std::vector<double>{11, 12});
    REQUIRE_CALL(*secondary, calculateForObservables(first_observable_values_eq({3, 4}))).IN_SEQUENCE(seq3);
    REQUIRE_CALL(*secondary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq3).RETURN(std::vector<double>{13, 14});
    REQUIRE_CALL(*secondary, calculateForObservables(first_observable_values_eq({5, 6}))).IN_SEQUENCE(seq3);
    REQUIRE_CALL(*secondary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq3).RETURN(std::vector<double>{15, 16});
    REQUIRE_CALL(*secondary, calculateForObservables(first_observable_values_eq({7, 8}))).IN_SEQUENCE(seq3);
    REQUIRE_CALL(*secondary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq3).RETURN(std::vector<double>{17, 18});

    trompeloeil::sequence seq5;
    std::vector<arma::cx_vec> vectors = {{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}};
    std::size_t i{};
    EvolverMock evolver;
    REQUIRE_CALL(evolver, prepareFor(arma_eq(arma::cx_vec{1, 0, 0, 0, 0}), 1, 2ul)).IN_SEQUENCE(seq5);
    REQUIRE_CALL(evolver, getDt()).TIMES(AT_LEAST(1)).RETURN(0.5).IN_SEQUENCE(seq5);
    REQUIRE_CALL(evolver, prepareFor(arma_eq(arma::cx_vec{0, 0, 1, 0, 0}), 2, 1ul)).IN_SEQUENCE(seq5);
    REQUIRE_CALL(evolver, getDt()).TIMES(AT_LEAST(1)).RETURN(2).IN_SEQUENCE(seq5);
    ALLOW_CALL(evolver, getCurrentState()).LR_RETURN(vectors.at(i));
    ALLOW_CALL(evolver, evolve()).LR_SIDE_EFFECT(i++);

    std::ostringstream loggerStream;
    Logger logger(loggerStream);
    OservablesTimeEvolution evolution;
    evolution.setPrimaryObservables({primary});
    evolution.setSecondaryObservables({secondary});
    evolution.setStoredObservables({primary, secondary});


    auto result = evolution.perform({{1, 2}, {3, 1}}, arma::cx_vec{1, 0, 0, 0, 0}, evolver, logger);


    REQUIRE(result == std::vector<TimeEvolutionEntry>{
        {0,   {1, 2, 11, 12}},
        {0.5, {3, 4, 13, 14}},
        {1,   {5, 6, 15, 16}},
        {3,   {7, 8, 17, 18}}
    });
}