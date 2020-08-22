//
// Created by Piotr Kubala on 21/08/2020.
//

#include <catch2/catch.hpp>
#include <catch2/trompeloeil.hpp>

#include "evolution/ChebyshevEvolution.h"

using namespace trompeloeil;

namespace {
    class HamiltonianGeneratorMock {
        MAKE_CONST_MOCK0(generate, arma::sp_mat());
    };

    class AveragingModelMock {
        MAKE_MOCK4(setupHamiltonianGenerator, void(HamiltonianGeneratorMock &, RND &, std::size_t, std::size_t));
    };

    class ChebyshevEvolverMock {
    public:
        arma::sp_mat hamiltonian;
        std::ostream &logger;

        ChebyshevEvolverMock(arma::sp_mat hamiltonian, std::ostream &logger)
                : hamiltonian{std::move(hamiltonian)}, logger{logger}
        { }
    };

    class CorrelationsTimeEvolutionMock {
    public:
        MAKE_MOCK3(addEvolution, void(ChebyshevEvolverMock &, std::ostream &, const std::vector<arma::cx_vec> &));
        MAKE_CONST_MOCK1(storeResult, void(std::ostream &));
        MAKE_CONST_MOCK0(countExternalVectors, std::size_t());
    };

    class QuenchCalculatorMock {
    public:
        MAKE_MOCK2(addQuench, void(const arma::sp_mat &, const arma::sp_mat &));
        MAKE_MOCK0(getLastQuenchedState, const arma::vec &());
        MAKE_CONST_MOCK0(getLastQuenchEpsilon, double());
        MAKE_CONST_MOCK0(getLastQuenchEpsilonQuantumUncertainty, double());
        MAKE_CONST_MOCK0(getMeanEpsilon, double());
        MAKE_CONST_MOCK0(getMeanEpsilonQuantumUncertainty, double());
        MAKE_CONST_MOCK0(getEpsilonAveragingSampleError, double());
        MAKE_CONST_MOCK0(clear, void());
    };

    bool spMatApproxEqual(const arma::sp_mat &mat1, const arma::sp_mat &mat2) {
        return arma::approx_equal(arma::mat(mat1), arma::mat(mat2), "absdiff", 1e-8);
    }

    bool cxVecApproxEqualVec(const arma::cx_vec &vec1, const arma::vec &vec2) {
        arma::cx_vec vec2Complex(vec2.size());
        std::copy(vec2.begin(), vec2.end(), vec2Complex.begin());
        return arma::approx_equal(vec1, vec2Complex, "absdiff", 1e-8);
    }

    using TestChebyshevEvolution = ChebyshevEvolution<HamiltonianGeneratorMock, AveragingModelMock,
                                                      CorrelationsTimeEvolutionMock, QuenchCalculatorMock,
                                                      ChebyshevEvolverMock>;
}

TEST_CASE("ChebyshevEvolutionTest") {
    SECTION("without quench - simulations 2, 3 from 5 total") {
        SimulationsSpan simulationsSpan;
        simulationsSpan.from = 2;
        simulationsSpan.to = 4;
        simulationsSpan.totalSimulations = 5;

        arma::sp_mat hamiltonian2(1, 1), hamiltonian3(1, 1);
        hamiltonian2(0, 0) = 1;
        hamiltonian3(0, 0) = 2;

        auto rnd = std::make_unique<RND>();
        auto rndPtr = rnd.get();
        auto hamiltonianGenerator = std::make_unique<HamiltonianGeneratorMock>();
        auto hamiltonianGeneratorPtr = hamiltonianGenerator.get();
        auto averagingModel = std::make_unique<AveragingModelMock>();
        auto correlationsTimeEvolution = std::make_unique<CorrelationsTimeEvolutionMock>();

        ALLOW_CALL(*correlationsTimeEvolution, countExternalVectors()).RETURN(0);

        sequence simulationSequence;
        REQUIRE_CALL(*averagingModel, setupHamiltonianGenerator(_, _, 2ul, 5ul))
                .WITH(&_1 == hamiltonianGeneratorPtr && &_2 == rndPtr)
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*hamiltonianGenerator, generate())
                .RETURN(hamiltonian2)
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*correlationsTimeEvolution, addEvolution(_, _, _))
                .WITH(spMatApproxEqual(_1.hamiltonian, hamiltonian2) && _3.empty())
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*averagingModel, setupHamiltonianGenerator(_, _, 3ul, 5ul))
                .WITH(&_1 == hamiltonianGeneratorPtr && &_2 == rndPtr)
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*hamiltonianGenerator, generate())
                .RETURN(hamiltonian3)
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*correlationsTimeEvolution, addEvolution(_, _, _))
                .WITH(spMatApproxEqual(_1.hamiltonian, hamiltonian3) && _3.empty())
                .IN_SEQUENCE(simulationSequence);

        std::ostringstream logger;
        TestChebyshevEvolution evolution(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                         simulationsSpan, std::move(correlationsTimeEvolution));


        evolution.perform(logger);
    }

    SECTION("with quench - one simulation") {
        SimulationsSpan simulationsSpan;
        simulationsSpan.from = 0;
        simulationsSpan.to = 1;
        simulationsSpan.totalSimulations = 1;

        arma::sp_mat hamiltonian(1, 1), quenchHamiltonian(1, 1);
        hamiltonian(0, 0) = 1;
        quenchHamiltonian(0, 0) = 2;
        arma::vec quenchedState = {1};

        auto rnd = std::make_unique<RND>();
        auto rndPtr = rnd.get();
        auto quenchRnd = std::make_unique<RND>();
        auto quenchRndPtr = quenchRnd.get();
        auto hamiltonianGenerator = std::make_unique<HamiltonianGeneratorMock>();
        auto hamiltonianGeneratorPtr = hamiltonianGenerator.get();
        auto quenchHamiltonianGenerator = std::make_unique<HamiltonianGeneratorMock>();
        auto quenchHamiltonianGeneratorPtr = quenchHamiltonianGenerator.get();
        auto averagingModel = std::make_unique<AveragingModelMock>();
        auto correlationsTimeEvolution = std::make_unique<CorrelationsTimeEvolutionMock>();
        auto quenchCalculator = std::make_unique<QuenchCalculatorMock>();

        ALLOW_CALL(*quenchCalculator, getLastQuenchEpsilon()).RETURN(0);
        ALLOW_CALL(*quenchCalculator, getLastQuenchEpsilonQuantumUncertainty()).RETURN(0);
        ALLOW_CALL(*quenchCalculator, getMeanEpsilon()).RETURN(0);
        ALLOW_CALL(*quenchCalculator, getMeanEpsilonQuantumUncertainty()).RETURN(0);
        ALLOW_CALL(*quenchCalculator, getEpsilonAveragingSampleError()).RETURN(0);
        ALLOW_CALL(*correlationsTimeEvolution, countExternalVectors()).RETURN(1);

        sequence simulationSequence;
        REQUIRE_CALL(*quenchCalculator, clear())
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*averagingModel, setupHamiltonianGenerator(_, _, 0ul, 1ul))
                .WITH(&_1 == hamiltonianGeneratorPtr && &_2 == rndPtr)
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*hamiltonianGenerator, generate())
                .RETURN(hamiltonian)
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*averagingModel, setupHamiltonianGenerator(_, _, 0ul, 1ul))
                .WITH(&_1 == quenchHamiltonianGeneratorPtr && &_2 == quenchRndPtr)
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*quenchHamiltonianGenerator, generate())
                .RETURN(quenchHamiltonian)
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*quenchCalculator, addQuench(_, _))
                .WITH(spMatApproxEqual(_1, quenchHamiltonian))
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*quenchCalculator, getLastQuenchedState())
                .RETURN(quenchedState)
                .IN_SEQUENCE(simulationSequence);
        REQUIRE_CALL(*correlationsTimeEvolution, addEvolution(_, _, _))
                .WITH(spMatApproxEqual(_1.hamiltonian, hamiltonian) &&
                      _3.size() == 1 && cxVecApproxEqualVec(_3[0], quenchedState))
                .IN_SEQUENCE(simulationSequence);

        std::ostringstream logger;
        TestChebyshevEvolution evolution(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                         simulationsSpan, std::move(correlationsTimeEvolution),
                                         std::move(quenchCalculator), std::move(quenchHamiltonianGenerator),
                                         std::move(quenchRnd));


        evolution.perform(logger);
    }
}