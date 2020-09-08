//
// Created by Piotr Kubala on 31/08/2020.
//

#ifndef MBL_ED_RESTORABLESIMULATIONEXECUTOR_H
#define MBL_ED_RESTORABLESIMULATIONEXECUTOR_H

#include <memory>
#include <filesystem>

#include "RestorableSimulation.h"
#include "utils/Logger.h"

/**
 * @brief A class performing a span of simulations, which can restore their state after an interruption and managing
 * simulations split between multiple processes.
 */
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
                                   Logger &logger);

    void doPerformSimulations(RestorableSimulation &simulation, const SimulationsSpan &actualSpan,
                              const std::string &stateFilename, Logger &logger) const;

    [[nodiscard]] SimulationStatus tryRestoringSimulation(RestorableSimulation &simulation,
                                                          const std::string &stateFilename, Logger &logger) const;

    [[nodiscard]] bool joinAllRestoredSimulations(RestorableSimulation &simulation,
                                                  const std::vector<StateFileData> &stateFileDatas) const;

    [[nodiscard]] std::string prepareStateFilenamePattern(const std::string &stateFilename) const;
    [[nodiscard]] std::vector<StateFileData> discoverStateFiles(const std::string &stateFilename) const;
    [[nodiscard]] StateFilesCoverage checkStateFilesCoverage(const std::vector<StateFileData> &stateFileDatas) const;
    SimulationStatus joinRestoredSimulations(RestorableSimulation &simulation, std::istream &binaryIn) const;

public:
    /**
     * @brief Construct the class.
     * @param simulationsSpan The span of simulations that should be performed.
     * @param fileSignature The signature of file used for storing state files.
     * @param splitWorkload If true, the class will work in workload split mode. See performSimulations.
     * @param workingDirectory The directory in which simulation state files will be stored/read from.
     */
    explicit RestorableSimulationExecutor(const SimulationsSpan &simulationsSpan, std::string fileSignature,
                                          bool splitWorkload, std::filesystem::path workingDirectory = "./");

    /**
     * @brief Intelligently perform a span of simulations.
     * @details <p> The method basicaly performs a span of simulations given by SimlationsSpan from the constructor.
     * This is the short story.
     * <p> The long story: first of all, it checks, whether there exists a state file corresponding to this simulation
     * span. If so, it concludes, that this span was already being performed and had been interrupted. So it loads
     * already done simulations and performs the rest. After each simulation the state is stored.
     * Simulation::seedRandomGenerators is invoked with @a seed increased by the first simulation index to be performed.
     * Note, that it would yield different results for interrupted vs not interrupted simulations, but it guarantees
     * that simulations won't be repeated.
     * <p> Later behaviour is determined by @a splitWorkload flag from the constructor. If false, the work is done,
     * state file is deleted, the method returns and shouldSaveSimulation() will return true. If the flag is true:
     * <ul>
     * <li> State file is not deleted.
     * <li> @a workingDirectory is searched for other state files. If they cover the whole
     * [0, SimulationsSpan::total - 1] range, they are all restored into (cleared) @a simulation in a proper,
     * ascending order. shouldSaveSimulation() will now return true. Moreover, all state files will be deleted.
     * <li> If any state files are missing, they contain interrupted simulations or they have broken, overlapping
     * ranges, @a simulation is cleared and shouldSaveSimulation() will return false. State files are not deleted -
     * either the last running process will finish the job or user intervention is required for interrupted simulations
     * or broken ranges.
     * </ul>
     */
    void performSimulations(RestorableSimulation &simulation, unsigned long seed, Logger &logger);

    /**
     * @brief After invoking performSimulations(), it indicates, weather results should be saved of they are not ready.
     * @details Unless @a splitWorkload from the constructor is true, it will always return true for non-interrupted
     * simulation.
     */
    [[nodiscard]] bool shouldSaveSimulation() const { return this->shouldSaveSimulation_; }
};


#endif //MBL_ED_RESTORABLESIMULATIONEXECUTOR_H
