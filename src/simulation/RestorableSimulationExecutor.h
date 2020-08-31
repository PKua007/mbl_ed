//
// Created by Piotr Kubala on 31/08/2020.
//

#ifndef MBL_ED_RESTORABLESIMULATIONEXECUTOR_H
#define MBL_ED_RESTORABLESIMULATIONEXECUTOR_H

#include <memory>

#include "RestorableSimulation.h"

class RestorableSimulationExecutor {
private:
    SimulationsSpan simulationsSpan;
    std::string fileSignature;

    void storeSimulations(const RestorableSimulation &simulation, std::ostream &binaryOut,
                          std::size_t simulationIndex, bool finished) const;
    std::size_t restoreSimulations(RestorableSimulation &simulation, std::istream &binaryIn) const;

public:
    explicit RestorableSimulationExecutor(const SimulationsSpan &simulationsSpan,
                                          std::string fileSignature);

    void performSimulations(RestorableSimulation &simulation, unsigned long seed, std::ostream &logger) const;
};


#endif //MBL_ED_RESTORABLESIMULATIONEXECUTOR_H
