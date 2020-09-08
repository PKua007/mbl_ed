//
// Created by Piotr Kubala on 31/08/2020.
//

#ifndef MBL_ED_RESTORABLESIMULATION_H
#define MBL_ED_RESTORABLESIMULATION_H

#include <string>

#include "Restorable.h"
#include "SimulationsSpan.h"
#include "utils/Logger.h"

/**
 * @brief An interface for a simulation accepted by RestorableSimulationExecutor.
 * @details It extends Restorable interface and also provides additional interface methods to control the simulations.
 */
class RestorableSimulation : public Restorable {
public:
    /**
     * @brief Seeds all random generators present in this simulation with @a seed.
     */
    virtual void seedRandomGenerators(unsigned long seed) = 0;

    /**
     * @brief Performs simulation of a given @a simulationIndex out of @a totalSimulations.
     * @details Subsequent simulations should be stored and be accessible via Restorable storing/restorin/clearing
     * methods. The method should to be called always with the same value of @a totalSimulations and subsequent
     * values of @a totalSimulations, with @a simulationsIndex < @a totalSimulations. However calling clear() should
     * allow another set ot potentially different simmulations.
     */
    virtual void performSimulation(std::size_t simulationIndex, std::size_t totalSimulations, Logger &logger) = 0;

    /**
     * @brief Return discernable class tag name used for state file name.
     */
    [[nodiscard]] virtual std::string getTagName() const = 0;
};


#endif //MBL_ED_RESTORABLESIMULATION_H
