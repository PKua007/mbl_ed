//
// Created by Piotr Kubala on 31/08/2020.
//

#include <utility>
#include <ostream>
#include <fstream>
#include <regex>

#include "RestorableSimulationExecutor.h"
#include "utils/Assertions.h"
#include "utils/Utils.h"

RestorableSimulationExecutor::RestorableSimulationExecutor(const SimulationsSpan &simulationsSpan,
                                                           std::string fileSignature, bool splitWorkload,
                                                           bool secureSimulationState,
                                                           std::filesystem::path workingDirectory)
        : simulationsSpan{simulationsSpan}, fileSignature{std::move(fileSignature)}, storeSimulations{secureSimulationState},
          splitWorkload{splitWorkload}, workingDirectory{std::move(workingDirectory)}
{
    Expects(simulationsSpan.total > 0);
    Expects(simulationsSpan.from < simulationsSpan.to);
    Expects(simulationsSpan.to <= simulationsSpan.total);
    if (splitWorkload)
        Expects(secureSimulationState);
}

void RestorableSimulationExecutor::performSimulations(RestorableSimulation &simulation, unsigned long seed,
                                                      Logger &logger)
{
    simulation.clear();
    this->shouldSaveSimulation_ = false;

    std::string stateFilename = this->fileSignature + "_state_" + simulation.getTagName() + ".bin";
    SimulationStatus simulationStatus = this->tryRestoringSimulation(simulation, stateFilename, logger);

    SimulationsSpan actualSpan;
    actualSpan = this->simulationsSpan;
    actualSpan.from = simulationStatus.nextSimulationIndex;
    simulation.seedRandomGenerators(seed + actualSpan.from);
    this->doPerformSimulations(simulation, actualSpan, stateFilename, logger);

    if (!this->splitWorkload) {
        this->shouldSaveSimulation_ = true;
        std::filesystem::remove(this->workingDirectory / stateFilename);
    } else {
        this->superviseSimulationsSplit(simulation, stateFilename, logger);
    }
}

RestorableSimulationExecutor::SimulationStatus
RestorableSimulationExecutor::tryRestoringSimulation(RestorableSimulation &simulation, const std::string &stateFilename,
                                                     Logger &logger) const
{
    SimulationStatus simulationStatus;
    std::ifstream stateFile(this->workingDirectory / stateFilename, std::ios::in | std::ios::binary);
    if (stateFile.is_open()) {
        simulationStatus = this->joinRestoredSimulations(simulation, stateFile);
        Assert(simulationStatus.nextSimulationIndex > this->simulationsSpan.from);
        Assert(simulationStatus.nextSimulationIndex <= this->simulationsSpan.to);

        logger.warn() << "State file found, restored simulations [" << this->simulationsSpan.from << ", ";
        logger << (simulationStatus.nextSimulationIndex - 1) << "], performing [";
        logger << simulationStatus.nextSimulationIndex << ", " << simulationsSpan.to << ") out of ";
        logger << simulationsSpan.total << std::endl;
    } else {
        logger.info() << "No state file found, starting simulations from scratch, [" << this->simulationsSpan.from;
        logger << ", " << this->simulationsSpan.to << ") out of " << this->simulationsSpan.total << std::endl;

        simulationStatus.finished = false;
        simulationStatus.nextSimulationIndex = this->simulationsSpan.from;
    }
    stateFile.close();

    return simulationStatus;
}

void RestorableSimulationExecutor::doPerformSimulations(RestorableSimulation &simulation,
                                                        const SimulationsSpan &actualSpan,
                                                        const std::string &stateFilename, Logger &logger) const
{
    for (std::size_t i = actualSpan.from; i < actualSpan.to; i++) {
        simulation.performSimulation(i, actualSpan.total, logger);

        if (this->storeSimulations) {
            std::ofstream storeFile(this->workingDirectory / stateFilename, std::ios::out | std::ios::binary);
            SimulationStatus simulationStatus;
            simulationStatus.finished = (i == actualSpan.to - 1);
            simulationStatus.nextSimulationIndex = i + 1;
            this->doStoreSimulations(simulation, storeFile, simulationStatus);
        }
    }
}

void RestorableSimulationExecutor::superviseSimulationsSplit(RestorableSimulation &simulation,
                                                             const std::string &stateFilename, Logger &logger)
{
    this->shouldSaveSimulation_ = false;

    logger.info() << "Workload split, searching for state files... " << std::flush;
    std::vector<StateFileData> stateFileDatas = this->discoverStateFiles(stateFilename);
    std::sort(stateFileDatas.begin(), stateFileDatas.end());
    logger << "found " << stateFileDatas.size() << " files." << std::endl;
    if (stateFileDatas.empty())
        return;

    switch (this->checkStateFilesCoverage(stateFileDatas)) {
        case INCOMPLETE:
            logger.info() << "Some state files are missing. Waiting for other tasks to finish." << std::endl;
            return;
        case BROKEN:
            logger.error() << "State files have broken ranges. Fix before reruning simulations." << std::endl;
            return;
        case COMPLETE:
            break;
    }
    logger.info() << "No state files are missing. Checking if all are finished." << std::endl;

    bool allSimulationsFinished = this->joinAllRestoredSimulations(simulation, stateFileDatas);
    this->shouldSaveSimulation_ = allSimulationsFinished;
    if (allSimulationsFinished) {
        logger.info() << "All simulations are finished. Removing state files." << std::endl;
        for (const auto &stateFileData : stateFileDatas)
            std::filesystem::remove(stateFileData.path);
    } else {
        logger.error() << "Some simulations must have been interrupted. Rerun the whole batch." << std::endl;
    }
}

bool RestorableSimulationExecutor::joinAllRestoredSimulations(RestorableSimulation &simulation,
                                                              const std::vector<StateFileData> &stateFileDatas) const
{
    simulation.clear();
    for (const auto &stateFileData : stateFileDatas) {
        std::ifstream file(stateFileData.path, std::ios::in | std::ios::binary);
        Assert(file);
        auto [finished, size] = joinRestoredSimulations(simulation, file);
        if (!finished) {
            simulation.clear();
            return false;

        }
    }
    return true;
}

RestorableSimulationExecutor::StateFilesCoverage RestorableSimulationExecutor
    ::checkStateFilesCoverage(const std::vector<RestorableSimulationExecutor::StateFileData> &stateFileDatas) const
{
    for (std::size_t i = 0; i < stateFileDatas.size() - 1; i++)
        if (stateFileDatas[i].to > stateFileDatas[i + 1].from)
            return StateFilesCoverage::BROKEN;

    if (stateFileDatas.front().from != 0)
        return StateFilesCoverage::INCOMPLETE;
    if (stateFileDatas.back().to != this->simulationsSpan.total)
        return StateFilesCoverage::INCOMPLETE;
    for (std::size_t i = 0; i < stateFileDatas.size() - 1; i++)
        if (stateFileDatas[i].to != stateFileDatas[i + 1].from)
            return StateFilesCoverage::INCOMPLETE;

    return StateFilesCoverage::COMPLETE;
}

std::vector<RestorableSimulationExecutor::StateFileData>
RestorableSimulationExecutor::discoverStateFiles(const std::string &stateFilename) const
{
    std::regex matcher(this->prepareStateFilenamePattern(stateFilename));

    std::vector<StateFileData> stateFileDatas;
    for (const auto &entry : std::filesystem::directory_iterator(workingDirectory)) {
        const std::filesystem::path &filePath = entry.path();
        std::smatch matches;
        std::string name = filePath.filename();
        if (std::regex_search(name, matches, matcher)) {
            Assert(matches.size() == 3);
            std::size_t from = std::stoul(matches[1].str());
            std::size_t to = std::stoul(matches[2].str());
            Assert(to > from);
            stateFileDatas.emplace_back(StateFileData{from, to, filePath});
        }
    }
    return stateFileDatas;
}

/**
 * @brief Given a @a stateFilename for any range, returns a regex pattern which captured @a from and @a to values.
 * @details Example:
 *    N.6_K.6_from.5_to.7_term.value
 * will be translated to
 *    N\.6_K\.6_from\.([0-9]+)_to\.([0-9]+)_term\.value
 */
std::string RestorableSimulationExecutor::prepareStateFilenamePattern(const std::string &stateFilename) const {
    std::regex specialChars(R"([-[\]{}()*+?.,\^$|#\s])");
    std::string pattern = std::regex_replace(stateFilename, specialChars, R"(\$&)");
    pattern = std::regex_replace(pattern, std::regex(R"(_from\\\.([0-9]+))"), R"(_from\.([0-9]+))");
    pattern = std::regex_replace(pattern, std::regex(R"(_to\\\.([0-9]+))"), R"(_to\.([0-9]+))");
    pattern = "^" + pattern + "$";
    return pattern;
}

void RestorableSimulationExecutor::doStoreSimulations(const RestorableSimulation &simulation, std::ostream &binaryOut,
                                                      const SimulationStatus &simulationStatus) const
{
    binaryOut.write(reinterpret_cast<const char*>(&simulationStatus), sizeof(simulationStatus));
    Assert(binaryOut);
    simulation.storeState(binaryOut);
}

RestorableSimulationExecutor::SimulationStatus
RestorableSimulationExecutor::joinRestoredSimulations(RestorableSimulation &simulation, std::istream &binaryIn) const
{
    SimulationStatus simulationStatus{};
    binaryIn.read(reinterpret_cast<char*>(&simulationStatus), sizeof(simulationStatus));
    Assert(binaryIn);
    simulation.joinRestoredState(binaryIn);
    return simulationStatus;
}
