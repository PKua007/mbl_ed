//
// Created by Piotr Kubala on 21/08/2020.
//

#include <sstream>

#include <catch2/catch.hpp>
#include <catch2/trompeloeil.hpp>

#include "simulation/ChebyshevEvolution.h"

#include "core/FockBaseGenerator.h"
#include "core/terms/HubbardHop.h"
#include "core/terms/HubbardOnsite.h"
#include "core/terms/QuasiperiodicDisorder.h"
#include "core/averaging_models/UniformPhi0AveragingModel.h"

#include "evolution/observables/OnsiteOccupations.h"

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

class CorrelationsTimeEvolutionMock : public trompeloeil::mock_interface<Restorable> {
    public:
        MAKE_MOCK3(addEvolution, void(ChebyshevEvolverMock &, std::ostream &, const std::vector<arma::cx_vec> &));
        MAKE_CONST_MOCK1(storeResult, void(std::ostream &));
        MAKE_CONST_MOCK0(countExternalVectors, std::size_t());
        IMPLEMENT_MOCK0(clear);
        IMPLEMENT_CONST_MOCK1(storeState);
        IMPLEMENT_MOCK1(joinRestoredState);
    };

    class QuenchCalculatorMock : public trompeloeil::mock_interface<Restorable> {
    public:
        MAKE_MOCK2(addQuench, void(const arma::sp_mat &, const arma::sp_mat &));
        MAKE_MOCK0(getLastQuenchedState, const arma::vec &());
        MAKE_CONST_MOCK0(getLastQuenchEpsilon, double());
        MAKE_CONST_MOCK0(getLastQuenchEpsilonQuantumUncertainty, double());
        MAKE_CONST_MOCK0(getMeanEpsilon, double());
        MAKE_CONST_MOCK0(getMeanEpsilonQuantumUncertainty, double());
        MAKE_CONST_MOCK0(getEpsilonAveragingSampleError, double());
        IMPLEMENT_MOCK0(clear);
        IMPLEMENT_CONST_MOCK1(storeState);
        IMPLEMENT_MOCK1(joinRestoredState);
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

TEST_CASE("ChebyshevEvolution: evolutions") {
    SECTION("without quench - simulations 2, 3 from 5 total") {
        arma::sp_mat hamiltonian2(1, 1);
        hamiltonian2(0, 0) = 1;

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

        std::ostringstream loggerStream;
        Logger logger(loggerStream);
        TestChebyshevEvolution evolution(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                         std::move(correlationsTimeEvolution));


        evolution.performSimulation(2, 5, logger);
    }

    SECTION("with quench - one simulation") {
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

        std::ostringstream loggerStream;
        Logger logger(loggerStream);
        TestChebyshevEvolution evolution(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                         std::move(correlationsTimeEvolution), std::move(quenchCalculator),
                                         std::move(quenchHamiltonianGenerator), std::move(quenchRnd));


        evolution.performSimulation(0, 1, logger);
    }
}

TEST_CASE("ChebyshevEvolution: clearing, storing and restoring") {
    // Actually we are doing integration test in the whole hierarchy of storing/restoring down CorrelationsTimeEvolution
    auto fockBase = std::shared_ptr(FockBaseGenerator().generate(4, 4));
    auto hamiltonianGenerator = std::make_unique<HamiltonianGenerator>(fockBase, false);
    hamiltonianGenerator->addHoppingTerm(std::make_unique<HubbardHop>(1));
    hamiltonianGenerator->addDiagonalTerm(std::make_unique<HubbardOnsite>(1));
    hamiltonianGenerator->addDiagonalTerm(std::make_unique<QuasiperiodicDisorder>(1, 0.3, 0));
    auto averagingModel = std::make_unique<UniformPhi0AveragingModel>();
    auto rnd = std::make_unique<RND>();
    TimeEvolutionParameters params;
    params.timeSegmentation = {{1, 1}};
    params.fockBase = fockBase;
    params.numberOfSites = 4;
    params.setVectorsToEvolveFromTags({"unif"});
    auto onsiteOccupations = std::make_shared<OnsiteOccupations>(4, fockBase);
    params.primaryObservables = {onsiteOccupations};
    params.storedObservables = {onsiteOccupations};
    auto occupationEvolution = std::make_unique<OservablesTimeEvolution>(
        params.primaryObservables, params.secondaryObservables, params.storedObservables
    );
    auto correlationsTimeEvolution = std::make_unique<TimeEvolution>(
        params, std::move(occupationEvolution)
    );
    std::ostringstream loggerStream;
    Logger logger(loggerStream);
    ChebyshevEvolution evolution(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                 std::move(correlationsTimeEvolution));

    SECTION("clearing") {
        evolution.performSimulation(0, 2, logger);
        std::ostringstream result1;
        evolution.storeResults(result1);

        evolution.performSimulation(1, 2, logger);
        evolution.clear();
        evolution.performSimulation(0, 2, logger);
        std::ostringstream result2;
        evolution.storeResults(result2);

        REQUIRE(result1.str() == result2.str());
    }

    SECTION("storing and joining restored") {
        evolution.performSimulation(0, 2, logger);
        evolution.performSimulation(1, 2, logger);
        std::ostringstream normalResult;
        evolution.storeResults(normalResult);

        evolution.clear();
        evolution.performSimulation(1, 2, logger);
        std::stringstream simulation1;
        evolution.storeState(simulation1);

        evolution.clear();
        evolution.performSimulation(0, 2, logger);
        evolution.joinRestoredState(simulation1);
        std::ostringstream restoredResult;
        evolution.storeResults(restoredResult);

        REQUIRE(normalResult.str() == restoredResult.str());
    }
}