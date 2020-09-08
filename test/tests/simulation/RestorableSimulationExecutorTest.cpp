//
// Created by Piotr Kubala on 01/09/2020.
//

#include <filesystem>
#include <vector>

#include <catch2/catch.hpp>
#include <ostream>

#include "simulation/RestorableSimulationExecutor.h"
#include "simulation/RestorableHelper.h"


namespace {
    class MockRestorableSimulation : public RestorableSimulation {
    private:
        unsigned long seed{};
        std::size_t interruptOn{};

    public:
        struct Data {
            std::size_t index{};
            std::size_t total{};
            unsigned long seed{};

            bool operator==(const Data &rhs) const {
                return std::tie(index, total, seed) == std::tie(rhs.index, rhs.total, rhs.seed);
            }

            friend std::ostream &operator<<(std::ostream &os, const Data &data) {
                os << "from: " << data.index << ", total: " << data.total << ", seed: " << data.seed;
                return os;
            }
        };

        std::vector<Data> simulations;

        explicit MockRestorableSimulation(unsigned long interruptOn = std::numeric_limits<unsigned long>::max())
                : interruptOn{interruptOn}
        { }

        void storeState(std::ostream &binaryOut) const override {
            RestorableHelper::storeStateForVector(this->simulations, binaryOut);
        }

        void joinRestoredState(std::istream &binaryIn) override {
            RestorableHelper::joinRestoredStateForVector(this->simulations, binaryIn);
        }

        void performSimulation(std::size_t simulationIndex, std::size_t totalSimulations, Logger&) override {
            if (simulationIndex == this->interruptOn)
                throw std::runtime_error("interruption");
            this->simulations.emplace_back(Data{simulationIndex, totalSimulations, this->seed});
        }

        void clear() override { this->simulations.clear(); }
        void seedRandomGenerators(unsigned long seed_) override { this->seed = seed_; };
        [[nodiscard]] std::string getTagName() const override { return "simulation"; }
    };
}

TEST_CASE("RestorableSimulationExecutor") {
    std::filesystem::path testDir = "RestorableSimulationExecutor_test/";
    std::filesystem::remove_all(testDir);
    std::filesystem::create_directory(testDir);

    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    using Simulations = std::vector<MockRestorableSimulation::Data>;

    SECTION("normal simulation [1, 3] out of 4") {
        MockRestorableSimulation restorableSimulation;
        RestorableSimulationExecutor executor({1, 4, 4}, "N.8_K.8_from.1_to.4_term.value", false, testDir);
        loggerStream.clear();

        executor.performSimulations(restorableSimulation, 1234, logger);

        CHECK(restorableSimulation.simulations == Simulations{{1, 4, 1235}, {2, 4, 1235}, {3, 4, 1235}});
        CHECK(executor.shouldSaveSimulation());
        CHECK_THAT(loggerStream.str(), Catch::Contains("No state file found"));
        CHECK(std::filesystem::is_empty(testDir));
    }

    SECTION("interrupted simulation [0, 1] + [2] out of 3") {
        MockRestorableSimulation restorableSimulation1(2);
        MockRestorableSimulation restorableSimulation2;
        RestorableSimulationExecutor executor({0, 3, 3}, "N.8_K.8_from.0_to.3_term.value", false, testDir);

        CHECK_THROWS_WITH(executor.performSimulations(restorableSimulation1, 1234, logger), "interruption");

        CHECK(restorableSimulation1.simulations == Simulations{{0, 3, 1234}, {1, 3, 1234}});
        CHECK_FALSE(executor.shouldSaveSimulation());


        loggerStream.clear();

        executor.performSimulations(restorableSimulation2, 1234, logger);

        CHECK(restorableSimulation2.simulations == Simulations{{0, 3, 1234}, {1, 3, 1234}, {2, 3, 1234 + 2}});
        CHECK(executor.shouldSaveSimulation());
        CHECK_THAT(loggerStream.str(), Catch::Contains("State file found"));
        CHECK(std::filesystem::is_empty(testDir));
    }

    SECTION("split simulation") {
        SECTION("first part [2] - not finished yet") {
            MockRestorableSimulation restorableSimulation1;
            RestorableSimulationExecutor executor1({2, 3, 3}, "N.8_K.8_from.2_to.3_term.value", true, testDir);
            loggerStream.clear();

            executor1.performSimulations(restorableSimulation1, 1234, logger);

            CHECK(restorableSimulation1.simulations == Simulations{{2, 3, 1234 + 2}});
            CHECK_FALSE(executor1.shouldSaveSimulation());
            CHECK_THAT(loggerStream.str(), Catch::Contains("Some state files are missing"));

            SECTION("second part [0, 1] - finished") {
                MockRestorableSimulation restorableSimulation2;
                RestorableSimulationExecutor executor2({0, 2, 3}, "N.8_K.8_from.0_to.2_term.value", true, testDir);
                loggerStream.clear();

                executor2.performSimulations(restorableSimulation2, 1234, logger);

                CHECK(restorableSimulation2.simulations == Simulations{{0, 3, 1234}, {1, 3, 1234}, {2, 3, 1234 + 2}});
                CHECK(executor2.shouldSaveSimulation());
                CHECK_THAT(loggerStream.str(), Catch::Contains("No state files are missing"));
                CHECK_THAT(loggerStream.str(), Catch::Contains("All simulations are finished"));
                CHECK(std::filesystem::is_empty(testDir));
            }
        }

        SECTION("[0] (interrupted) + [2] - all files, but first not ready, so unfinished") {
            MockRestorableSimulation restorableSimulation1(1);
            RestorableSimulationExecutor executor1({0, 2, 3}, "N.8_K.8_from.0_to.2_term.value", true, testDir);
            loggerStream.clear();

            CHECK_THROWS_WITH(executor1.performSimulations(restorableSimulation1, 1234, logger), "interruption");

            CHECK(restorableSimulation1.simulations == Simulations{{0, 3, 1234}});
            CHECK_FALSE(executor1.shouldSaveSimulation());


            MockRestorableSimulation restorableSimulation2;
            RestorableSimulationExecutor executor2({2, 3, 3}, "N.8_K.8_from.2_to.3_term.value", true, testDir);
            loggerStream.clear();

            executor2.performSimulations(restorableSimulation2, 1234, logger);

            CHECK(restorableSimulation2.simulations.empty());
            CHECK_FALSE(executor2.shouldSaveSimulation());
            CHECK_THAT(loggerStream.str(), Catch::Contains("No state files are missing"));
            CHECK_THAT(loggerStream.str(), Catch::Contains("Some simulations must have been interrupted"));

            SECTION("redoing [0, 1] - finishing") {
                MockRestorableSimulation restorableSimulation3;
                RestorableSimulationExecutor executor3({0, 2, 3}, "N.8_K.8_from.0_to.2_term.value", true, testDir);
                loggerStream.clear();

                executor3.performSimulations(restorableSimulation3, 1234, logger);

                CHECK(restorableSimulation3.simulations == Simulations{{0, 3, 1234}, {1, 3, 1234 + 1}, {2, 3, 1234+2}});
                CHECK(executor3.shouldSaveSimulation());
                CHECK_THAT(loggerStream.str(), Catch::Contains("No state files are missing"));
                CHECK_THAT(loggerStream.str(), Catch::Contains("All simulations are finished"));
                CHECK(std::filesystem::is_empty(testDir));
            }
        }

        SECTION("[0, 1] + [0, 2] - not finished, because messed up ranges") {
            MockRestorableSimulation restorableSimulation1;
            RestorableSimulationExecutor executor1({0, 2, 3}, "N.8_K.8_from.0_to.2_term.value", true, testDir);
            loggerStream.clear();

            executor1.performSimulations(restorableSimulation1, 1234, logger);

            CHECK(restorableSimulation1.simulations == Simulations{{0, 3, 1234}, {1, 3, 1234}});
            CHECK_FALSE(executor1.shouldSaveSimulation());
            CHECK_THAT(loggerStream.str(), Catch::Contains("Some state files are missing"));


            MockRestorableSimulation restorableSimulation2;
            RestorableSimulationExecutor executor2({0, 3, 3}, "N.8_K.8_from.0_to.3_term.value", true, testDir);
            loggerStream.clear();

            executor2.performSimulations(restorableSimulation2, 1234, logger);

            CHECK(restorableSimulation2.simulations == Simulations{{0, 3, 1234}, {1, 3, 1234}, {2, 3, 1234}});
            CHECK_FALSE(executor2.shouldSaveSimulation());
            CHECK_THAT(loggerStream.str(), Catch::Contains("State files have broken ranges"));
        }
    }

    std::filesystem::remove_all(testDir);
}