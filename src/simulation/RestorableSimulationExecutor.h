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

public:
    explicit RestorableSimulationExecutor(const SimulationsSpan &simulationsSpan);

    void performSimulations(RestorableSimulation &simulation, unsigned long seed, std::ostream &logger) const;
};


#endif //MBL_ED_RESTORABLESIMULATIONEXECUTOR_H
