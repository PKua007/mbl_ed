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
                                                           std::filesystem::path workingDirectory)
        : simulationsSpan{simulationsSpan}, fileSignature{std::move(fileSignature)}, splitWorkload{splitWorkload},
          workingDirectory{std::move(workingDirectory)}
{
    Expects(simulationsSpan.total > 0);
    Expects(simulationsSpan.from < simulationsSpan.to);
    Expects(simulationsSpan.to <= simulationsSpan.total);
}

void RestorableSimulationExecutor::performSimulations(RestorableSimulation &simulation, unsigned long seed,
                                                      std::ostream &logger)
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
                                                     std::ostream &logger) const
{
    SimulationStatus simulationStatus;
    std::ifstream stateFile(this->workingDirectory / stateFilename, std::ios::in | std::ios::binary);
    if (stateFile.is_open()) {
        simulationStatus = joinRestoredSimulations(simulation, stateFile);
        Assert(simulationStatus.nextSimulationIndex > this->simulationsSpan.from);
        Assert(simulationStatus.nextSimulationIndex <= this->simulationsSpan.to);

        logger << "[RestorableSimulationExecutor::tryRestoringSimulation] State file found, restored simulations [";
        logger << this->simulationsSpan.from << ", " << (simulationStatus.nextSimulationIndex - 1);
        logger << "], performing [" << simulationStatus.nextSimulationIndex << ", " << simulationsSpan.to;
        logger << ") out of " << simulationsSpan.total << std::endl;
    } else {
        logger << "[RestorableSimulationExecutor::tryRestoringSimulation] No state file found, starting simulations from ";
        logger << "scratch, [" << this->simulationsSpan.from << ", " << this->simulationsSpan.to << ") out of ";
        logger << this->simulationsSpan.total << std::endl;

        simulationStatus.finished = false;
        simulationStatus.nextSimulationIndex = this->simulationsSpan.from;
    }
    stateFile.close();

    return simulationStatus;
}

void RestorableSimulationExecutor::doPerformSimulations(RestorableSimulation &simulation,
                                                        const SimulationsSpan &actualSpan,
                                                        const std::string &stateFilename, std::ostream &logger) const
{
    for (std::size_t i = actualSpan.from; i < actualSpan.to; i++) {
        simulation.performSimulation(i, actualSpan.total, logger);
        std::ofstream storeFile(workingDirectory / stateFilename, std::ios::out | std::ios::binary);
        SimulationStatus simulationStatus;
        simulationStatus.finished = (i == actualSpan.to - 1);
        simulationStatus.nextSimulationIndex = i + 1;
        storeSimulations(simulation, storeFile, simulationStatus);
    }
}

void RestorableSimulationExecutor::superviseSimulationsSplit(RestorableSimulation &simulation,
                                                             const std::string &stateFilename, std::ostream &logger)
{
    this->shouldSaveSimulation_ = false;

    logger << "[RestorableSimulationExecutor::superviseSimulationsSplit] Workload split, searching for state files... ";
    logger << std::flush;
    std::vector<StateFileData> stateFileDatas = this->discoverStateFiles(stateFilename);
    std::sort(stateFileDatas.begin(), stateFileDatas.end());
    logger << "found " << stateFileDatas.size() << " files." << std::endl;
    if (stateFileDatas.empty())
        return;

    switch (this->checkStateFilesCoverage(stateFileDatas)) {
        case INCOMPLETE:
            logger << "[RestorableSimulationExecutor::superviseSimulationsSplit] Some state files are missing. ";
            logger << "Waiting for other tasks to finish." << std::endl;
            return;
        case BROKEN:
            logger << "[RestorableSimulationExecutor::superviseSimulationsSplit] State files have broken ranges. ";
            logger << "Fix before reruning simulations." << std::endl;
            return;
        case COMPLETE:
            break;
    }
    logger << "[RestorableSimulationExecutor::superviseSimulationsSplit] No state files are missing. ";
    logger << "Checking if all are finished." << std::endl;

    bool allSimulationsFinished = this->joinAllRestoredSimulations(simulation, stateFileDatas);
    this->shouldSaveSimulation_ = allSimulationsFinished;
    if (allSimulationsFinished) {
        logger << "[RestorableSimulationExecutor::superviseSimulationsSplit] All simulations are finished. ";
        logger << "Removing state files." << std::endl;
        for (const auto &stateFileData : stateFileDatas)
            std::filesystem::remove(stateFileData.path);
    } else {
        logger << "[RestorableSimulationExecutor::superviseSimulationsSplit] Some simulations must have been ";
        logger << "interrupted. Rerun the whole batch." << std::endl;
    }
}

bool RestorableSimulationExecutor::joinAllRestoredSimulations(RestorableSimulation &simulation,
                                                              const std::vector<StateFileData> &stateFileDatas) const {
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
    std::regex matcher(prepareStateFilenamePattern(stateFilename));

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

std::string RestorableSimulationExecutor::prepareStateFilenamePattern(const std::string &stateFilename) const {
    std::regex specialChars(R"([-[\]{}()*+?.,\^$|#\s])");
    std::string pattern = std::regex_replace(stateFilename, specialChars, R"(\$&)");
    pattern = std::regex_replace(pattern, std::regex(R"(_from\\\.([0-9]+))"), R"(_from\.([0-9]+))");
    pattern = std::regex_replace(pattern, std::regex(R"(_to\\\.([0-9]+))"), R"(_to\.([0-9]+))");
    pattern = "^" + pattern + "$";
    return pattern;
}

void RestorableSimulationExecutor::storeSimulations(const RestorableSimulation &simulation, std::ostream &binaryOut,
                                                    const SimulationStatus &status) const
{
    binaryOut.write(reinterpret_cast<const char*>(&status), sizeof(status));
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
