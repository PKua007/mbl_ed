//
// Created by pkua on 13.11.2019.
//

#include <sstream>

#include <catch2/catch.hpp>
#include <catch2/trompeloeil.hpp>

#include "mocks/FileUtilsMock.h"

#include "core/HamiltonianGenerator.h"
#include "simulation/ExactDiagonalization.h"

using trompeloeil::_;

namespace {
    class MockHamiltonianGenerator {
    public:
        MAKE_CONST_MOCK1(calculateEigensystem, Eigensystem(bool));
        MAKE_CONST_MOCK0(getFockBase, std::shared_ptr<const FockBase>());
    };

    class MockAveragingModel {
    public:
        MAKE_CONST_MOCK4(setupHamiltonianGenerator, void(MockHamiltonianGenerator &, RND &, std::size_t, std::size_t));
    };

    class MockAnalyzer : public trompeloeil::mock_interface<Restorable> {
    public:
        MAKE_CONST_MOCK2(analyze, void(const Eigensystem &, std::ostream &));
        IMPLEMENT_CONST_MOCK1(storeState);
        IMPLEMENT_MOCK1(joinRestoredState);
        IMPLEMENT_MOCK0(clear);
    };
}

TEST_CASE("ExactDiagonlization: 3 'random' hamiltonians") {
    ExactDiagonalizationParameters params;
    params.calculateEigenvectors = true;
    params.saveEigenenergies = false;
    params.fileSignature = "";
    Eigensystem eigensystem({-1, 1, 2}, {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}});

    auto hamiltonianGenerator = std::make_unique<MockHamiltonianGenerator>();
    auto hamiltonianGeneratorPtr = hamiltonianGenerator.get();
    auto rnd = std::make_unique<RND>();
    auto rndPtr = rnd.get();
    auto averagingModel = std::make_unique<MockAveragingModel>();
    auto analyzer = std::make_unique<MockAnalyzer>();
    trompeloeil::sequence seq;
    REQUIRE_CALL(*averagingModel, setupHamiltonianGenerator(_, _, 1ul, 3ul))
        .WITH(&_1 == hamiltonianGeneratorPtr && &_2 == rndPtr)
        .IN_SEQUENCE(seq);
    REQUIRE_CALL(*hamiltonianGenerator, calculateEigensystem(true))
        .RETURN(eigensystem)
        .IN_SEQUENCE(seq);
    REQUIRE_CALL(*analyzer, analyze(eigensystem, _))
        .IN_SEQUENCE(seq);

    using TestSimulation = ExactDiagonalization<MockHamiltonianGenerator, MockAveragingModel, MockAnalyzer>;
    TestSimulation simulation(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                              std::make_unique<FileOstreamProviderMock>(), params, std::move(analyzer));
    std::ostringstream dummyLoggerStream;
    Logger dummyLogger(dummyLoggerStream);

    simulation.performSimulation(1, 3, dummyLogger);
}