//
// Created by pkua on 13.11.2019.
//

#include <catch2/catch.hpp>
#include <sstream>

#include "matchers/ArmaApproxEqualMatcher.h"

#include "simulation/HamiltonianGenerator.h"
#include "simulation/Simulation.h"

namespace {
    class MockHamiltonianGenerator {
    private:
        std::size_t simulationIndex{};
        std::vector<arma::mat> matrices{};

        friend class MockAveragingModel;

    public:
        MockHamiltonianGenerator() {
            // -1, 0, 9
            matrices.push_back({{5, 2, 4},
                                {2, 0, 2},
                                {4, 2, 3}});
            // -1, 1, 2
            matrices.push_back({{0, 0, 1},
                                {0, 2, 0},
                                {1, 0, 0}});
            // 1, 2, 3
            matrices.push_back({{1, 0, 0},
                                {0, 2, 0},
                                {0, 0, 3}});
        }

        [[nodiscard]] arma::mat generate() const {
            return this->matrices.at(simulationIndex - 1);
        }

        [[nodiscard]] std::shared_ptr<const FockBase> getFockBase() const {
            return nullptr;
        }
    };

    class MockAveragingModel {
    public:
        static void setupHamiltonianGenerator(MockHamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                              std::size_t simulationIndex, std::size_t totalSimulations)
        {
            static_cast<void>(rnd);
            REQUIRE((simulationIndex >= 1 && simulationIndex <= 3));
            REQUIRE(totalSimulations == 5);
            hamiltonianGenerator.simulationIndex = simulationIndex;
        }
    };

    class MockAnalyzer {
    private:
        std::size_t index = 1;

    public:
        ~MockAnalyzer() {
            REQUIRE(index == 4);
        }

        void analyze(const Eigensystem &eigensystem, std::ostream &logger) {
            static_cast<void>(logger);

            REQUIRE(index >= 1);
            REQUIRE(index <= 3);
            const auto &eigenenergies = eigensystem.getEigenenergies();
            if (index == 1)
                REQUIRE_THAT(eigenenergies, IsApproxEqual(arma::vec{-1, 0, 9}, 1e-8));
            else if (index == 2)
                REQUIRE_THAT(eigenenergies, IsApproxEqual(arma::vec{-1, 1, 2}, 1e-8));
            else if (index == 3)
                REQUIRE_THAT(eigenenergies, IsApproxEqual(arma::vec{1, 2, 3}, 1e-8));
            index++;
        }
    };

    struct DummyOstreamProvider : public FileOstreamProvider {
    public:
        [[nodiscard]] std::unique_ptr<std::ostream> openFile(const std::string &filename) const override {
            FAIL("tried to open file: " + filename);
            return nullptr;
        }
    };
}

TEST_CASE("Simulation: 3 'random' hamiltonians") {
    SimulationParameters params;
    params.from = 1;
    params.to = 4;
    params.totalSimulations = 5;
    params.calculateEigenvectors = true;
    params.saveEigenenergies = false;
    params.fileSignature = "";
    using TestSimulation = Simulation<MockHamiltonianGenerator, MockAveragingModel, MockAnalyzer>;
    TestSimulation simulation(std::make_unique<MockHamiltonianGenerator>(), std::make_unique<MockAveragingModel>(),
                              std::make_unique<RND>(), std::make_unique<DummyOstreamProvider>(), params);
    MockAnalyzer mockAnalyzer;
    std::ostringstream dummyLogger;

    simulation.perform(dummyLogger, mockAnalyzer);
}