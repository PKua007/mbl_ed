//
// Created by Piotr Kubala on 31/08/2020.
//

#include "RestorableSimulationExecutor.h"

#include "utils/Assertions.h"

RestorableSimulationExecutor::RestorableSimulationExecutor(const SimulationsSpan &simulationsSpan)
        : simulationsSpan{simulationsSpan}
{
    Expects(simulationsSpan.total > 0);
    Expects(simulationsSpan.from < simulationsSpan.to);
    Expects(simulationsSpan.to <= simulationsSpan.total);
}

void RestorableSimulationExecutor::performSimulations(RestorableSimulation &simulation, unsigned long seed,
                                                      std::ostream &logger) const
{
    simulation.clear();
    simulation.seedRandomGenerators(seed + this->simulationsSpan.from);
    simulation.performSimulations(this->simulationsSpan, logger);
}
