//
// Created by Piotr Kubala on 31/08/2020.
//

#ifndef MBL_ED_RESTORABLESIMULATIONEXECUTOR_H
#define MBL_ED_RESTORABLESIMULATIONEXECUTOR_H

#include <memory>
#include <filesystem>

#include "RestorableSimulation.h"

class RestorableSimulationExecutor {
private:
    struct SimulationStatus {
        bool finished{};
        std::size_t nextSimulationIndex{};
    };

    struct StateFileData {
        std::size_t from{};
        std::size_t to{};
        std::filesystem::path path;

        friend bool operator<(const StateFileData &lhs, const StateFileData &rhs) {
            return std::tie(lhs.from, lhs.to) < std::tie(rhs.from, rhs.to);
        }
    };

    enum StateFilesCoverage {
        COMPLETE,
        INCOMPLETE,
        BROKEN
    };

    SimulationsSpan simulationsSpan;
    std::string fileSignature;
    bool splitWorkload{};
    bool shouldSaveSimulation_{};
    std::filesystem::path workingDirectory;

    void storeSimulations(const RestorableSimulation &simulation, std::ostream &binaryOut,
                          const SimulationStatus &simulationStatus) const;
    void superviseSimulationsSplit(RestorableSimulation &simulation, const std::string &stateFilename,
                                   std::ostream &logger);

    void doPerformSimulations(RestorableSimulation &simulation, const SimulationsSpan &actualSpan,
                              const std::string &stateFilename, std::ostream &logger) const;

    [[nodiscard]] SimulationStatus tryRestoringSimulation(RestorableSimulation &simulation,
                                                          const std::string &stateFilename, std::ostream &logger) const;

    [[nodiscard]] bool joinAllRestoredSimulations(RestorableSimulation &simulation,
                                                  const std::vector<StateFileData> &stateFileDatas) const;

    [[nodiscard]] std::string prepareStateFilenamePattern(const std::string &stateFilename) const;
    [[nodiscard]] std::vector<StateFileData> discoverStateFiles(const std::string &stateFilename) const;
    [[nodiscard]] StateFilesCoverage checkStateFilesCoverage(const std::vector<StateFileData> &stateFileDatas) const;
    SimulationStatus joinRestoredSimulations(RestorableSimulation &simulation, std::istream &binaryIn) const;

public:
    explicit RestorableSimulationExecutor(const SimulationsSpan &simulationsSpan, std::string fileSignature,
                                          bool splitWorkload, std::filesystem::path workingDirectory = "./");

    void performSimulations(RestorableSimulation &simulation, unsigned long seed, std::ostream &logger);
    [[nodiscard]] bool shouldSaveSimulation() const { return this->shouldSaveSimulation_; }
};


#endif //MBL_ED_RESTORABLESIMULATIONEXECUTOR_H
