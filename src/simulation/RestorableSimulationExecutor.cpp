//
// Created by Piotr Kubala on 31/08/2020.
//

#include <utility>
#include <ostream>
#include <fstream>
#include <regex>
#include <numeric>

#include <ZipIterator.hpp>

#include "RestorableSimulationExecutor.h"
#include "utils/Assertions.h"
#include "utils/Utils.h"

RestorableSimulationExecutor::RestorableSimulationExecutor(const SimulationsSpan &simulationsSpan,
                                                           std::string fileSignature, bool workloadSplit,
                                                           std::filesystem::path workingDirectory)
        : simulationsSpan{simulationsSpan}, fileSignature{std::move(fileSignature)}, workloadSplit{workloadSplit},
          workingDirectory{std::move(workingDirectory)}
{
    Expects(simulationsSpan.total > 0);
    Expects(simulationsSpan.from < simulationsSpan.to);
    Expects(simulationsSpan.to <= simulationsSpan.total);
}

void RestorableSimulationExecutor::performSimulations(RestorableSimulation &simulation, unsigned long seed,
                                                      std::ostream &logger)
{
    std::string filename = this->fileSignature + "_state_" + simulation.getTagName() + ".dat";
    std::ifstream restoreFile(this->workingDirectory / filename, std::ios::in | std::ios::binary);
    SimulationsSpan actualSpan;

    simulation.clear();

    if (restoreFile.is_open()) {
        auto [finished, simulationIndex] = this->joinRestoredSimulations(simulation, restoreFile);
        if (finished) {
            logger << "That simulation already done. Aborting" << std::endl;
            this->shouldSave_ = false;
            return;
        }

        Assert(simulationIndex >= this->simulationsSpan.from);
        Assert(simulationIndex <= this->simulationsSpan.to - 2);
        logger << "[RestorableSimulationExecutor::performSimulations] State file found, restored simulations [";
        logger << this->simulationsSpan.from << ", " << simulationIndex << "], performing [";
        logger << (simulationIndex + 1) << ", " << this->simulationsSpan.to << ") out of ";
        logger << this->simulationsSpan.total << std::endl;
        actualSpan = this->simulationsSpan;
        actualSpan.from = simulationIndex + 1;
    } else {
        logger << "[RestorableSimulationExecutor::performSimulations] No state file found, starting simulations from ";
        logger << "scratch, [" << this->simulationsSpan.from << ", " << this->simulationsSpan.to << ") out of ";
        logger << this->simulationsSpan.total << std::endl;
        actualSpan = this->simulationsSpan;
    }

    simulation.seedRandomGenerators(seed + actualSpan.from);

    for (std::size_t i = actualSpan.from; i < actualSpan.to; i++) {
        simulation.performSimulation(i, actualSpan.total, logger);
        std::ofstream storeFile(this->workingDirectory / filename, std::ios::out | std::ios::binary);
        bool finished = (i == actualSpan.to - 1);
        this->storeSimulations(simulation, storeFile, i, finished);
    }

    if (!this->workloadSplit) {
        logger << "Workload not split, finishing" << std::endl;
        this->shouldSave_ = true;
        std::filesystem::remove(this->workingDirectory / filename);
    } else {
        logger << "Workload split, will search for files" << std::endl;
        this->shouldSave_ = false;

        std::regex specialChars(R"([-[\]{}()*+?.,\^$|#\s])");
        std::string pattern = std::regex_replace(filename, specialChars, R"(\$&)");
        pattern = std::regex_replace(pattern, std::regex(R"(_from\\\.([0-9]+))"), R"(_from\.([0-9]+))");
        pattern = std::regex_replace(pattern, std::regex(R"(_to\\\.([0-9]+))"), R"(_to\.([0-9]+))");
        pattern = "^" + pattern + "$";
        std::regex matcher(pattern);

        std::vector<std::filesystem::path> files;
        std::vector<std::size_t> froms;
        std::vector<std::size_t> tos;

        for (const auto &entry : std::filesystem::directory_iterator(this->workingDirectory)) {
            const std::filesystem::path &filePath = entry.path();
            std::smatch matches;
            std::string name = filePath.filename();
            if (std::regex_search(name, matches, matcher)) {
                files.push_back(filePath);
                Assert(matches.size() == 3);
                std::size_t from = std::stoul(matches[1].str());
                std::size_t to = std::stoul(matches[2].str());
                Assert(to > from);
                froms.push_back(from);
                tos.push_back(to);
            }
        }

        if (!froms.empty()) {
            logger << "Found " << froms.size() << "files" << std::endl;

            auto zipped = Zip(froms, tos, files);
            std::sort(zipped.begin(), zipped.end());

            bool correct = true;
            if (froms.front() != 0)
                correct = false;
            if (tos.back() != this->simulationsSpan.total)
                correct = false;
            for (std::size_t i = 0; i < froms.size() - 1; i++)
                if (tos[i] != froms[i + 1])
                    correct = false;

            if (correct) {
                logger << "All files are there. Checking if all finished" << std::endl;

                simulation.clear();

                bool allFinished = true;
                for (const auto &path : files) {
                    std::ifstream file(path, std::ios::in | std::ios::binary);
                    Assert(file);
                    auto [finished, size] = this->joinRestoredSimulations(simulation, file);
                    if (!finished) {
                        allFinished = false;
                        break;
                    }
                }

                logger << "Finished? " << allFinished << std::endl;

                this->shouldSave_ = allFinished;

                if (allFinished) {
                    logger << "Removing state files" << std::endl;
                    for (const auto &file : files)
                        std::filesystem::remove(file);
                } else {
                    simulation.clear();
                }
            }
        } else {
            logger << "Found " << froms.size() << "no files" << std::endl;
        }
    }
}

void RestorableSimulationExecutor::storeSimulations(const RestorableSimulation &simulation, std::ostream &binaryOut,
                                                    std::size_t simulationIndex, bool finished) const
{
    binaryOut.write(reinterpret_cast<const char*>(&finished), sizeof(finished));
    binaryOut.write(reinterpret_cast<const char*>(&simulationIndex), sizeof(simulationIndex));
    Assert(binaryOut);
    simulation.storeState(binaryOut);
}

std::pair<bool, std::size_t> RestorableSimulationExecutor::joinRestoredSimulations(RestorableSimulation &simulation,
                                                                                   std::istream &binaryIn) const
{
    bool finished{};
    std::size_t simulationIndex{};
    binaryIn.read(reinterpret_cast<char*>(&finished), sizeof(finished));
    binaryIn.read(reinterpret_cast<char*>(&simulationIndex), sizeof(simulationIndex));
    Assert(binaryIn);
    simulation.joinRestoredState(binaryIn);
    return std::make_pair(finished, simulationIndex);
}
