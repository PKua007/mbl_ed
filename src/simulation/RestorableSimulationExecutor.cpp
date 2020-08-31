//
// Created by Piotr Kubala on 31/08/2020.
//

#include <utility>
#include <ostream>
#include <istream>
#include <fstream>
#include <filesystem>

#include "RestorableSimulationExecutor.h"
#include "utils/Assertions.h"

RestorableSimulationExecutor::RestorableSimulationExecutor(const SimulationsSpan &simulationsSpan,
                                                           std::string fileSignature)
        : simulationsSpan{simulationsSpan}, fileSignature{std::move(fileSignature)}
{
    Expects(simulationsSpan.total > 0);
    Expects(simulationsSpan.from < simulationsSpan.to);
    Expects(simulationsSpan.to <= simulationsSpan.total);
}

void RestorableSimulationExecutor::performSimulations(RestorableSimulation &simulation, unsigned long seed,
                                                      std::ostream &logger) const
{
    std::string filename = fileSignature + "_state_" + simulation.getTagName() + ".dat";
    std::ifstream restoreFile(filename, std::ios::in | std::ios::binary);
    SimulationsSpan actualSpan;

    if (restoreFile.is_open()) {
        std::size_t simulationIndex = this->restoreSimulations(simulation, restoreFile);
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
    simulation.clear();

    for (std::size_t i = actualSpan.from; i < actualSpan.to; i++) {
        simulation.performSimulation(i, actualSpan.total, logger);
        std::ofstream storeFile(filename, std::ios::out | std::ios::binary);
        bool finished = (i == actualSpan.to - 1);
        this->storeSimulations(simulation, storeFile, i, finished);
    }

    std::filesystem::remove(filename);
}

void RestorableSimulationExecutor::storeSimulations(const RestorableSimulation &simulation, std::ostream &binaryOut,
                                                    std::size_t simulationIndex, bool finished) const
{
    binaryOut.write(reinterpret_cast<const char*>(&finished), sizeof(finished));
    binaryOut.write(reinterpret_cast<const char*>(&simulationIndex), sizeof(simulationIndex));
    Assert(binaryOut);
    simulation.storeState(binaryOut);
}

std::size_t RestorableSimulationExecutor::restoreSimulations(RestorableSimulation &simulation,
                                                             std::istream &binaryIn) const
{
    bool finished{};
    std::size_t simulationIndex{};
    binaryIn.read(reinterpret_cast<char*>(&finished), sizeof(finished));
    binaryIn.read(reinterpret_cast<char*>(&simulationIndex), sizeof(simulationIndex));
    Assert(binaryIn);
    simulation.restoreState(binaryIn);
    return simulationIndex;
}
