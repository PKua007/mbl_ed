//
// Created by pkua on 13.11.2019.
//

#include <Catch2/catch.hpp>
#include <sstream>
#include <HamiltonianGenerator.h>

#include "Simulation.h"

namespace {
    bool approxVectorsEqual(std::vector<double> first, std::vector<double> second) {
        if (first.size() != second.size())
            return false;
        for (std::size_t i{}; i < first.size(); i++)
            if (std::abs(first[i] - second[i]) > 0.0000001)
                return false;
        return true;
    }

    class MockHamiltonianGenerator {
    private:
        std::size_t timesResampled{};
        std::vector<arma::mat> matrices{};

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
            return this->matrices.at(timesResampled - 1);
        }

        void resampleOnsiteEnergies() {
            this->timesResampled++;
        }

        std::string fileSignature() const {
            return "";
        }
    };

    class MockGapRatioCalculator {
    private:
        std::size_t index{};

    public:
        MockGapRatioCalculator(double relativeMiddleEnergy, double relativeMargin) {
            REQUIRE(relativeMiddleEnergy == 0.5);
            REQUIRE(relativeMargin == 0.99);
        }

        ~MockGapRatioCalculator() {
            REQUIRE(index == 3);
        }

        void addEigenenergies(const std::vector<double> &eigenenergies) {
            REQUIRE(index < 3);
            if (index == 0)
                REQUIRE(approxVectorsEqual(eigenenergies, std::vector<double>{-1, 0, 9}));
            else if (index == 1)
                REQUIRE(approxVectorsEqual(eigenenergies, std::vector<double>{-1, 1, 2}));
            else if (index == 2)
                REQUIRE(approxVectorsEqual(eigenenergies, std::vector<double>{1, 2, 3}));
            index++;
        }

        Quantity calculateMean() {
            return Quantity{};
        }
    };

    struct DummyOstreamProvider {
    public:
        static std::ostringstream openFile(const std::string &filename) {
            static_cast<void>(filename);
            return std::ostringstream{};
        }
    };
}


TEST_CASE("Simulation: 3 'random' hamiltonians") {
    using TestSimulation = Simulation<MockHamiltonianGenerator, MockGapRatioCalculator, DummyOstreamProvider>;
    TestSimulation simulation(std::make_unique<MockHamiltonianGenerator>(),3, 0.5, 0.99);
    std::ostringstream dummyLogger;

    simulation.perform(dummyLogger);
}