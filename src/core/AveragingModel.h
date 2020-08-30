//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_AVERAGINGMODEL_H
#define MBL_ED_AVERAGINGMODEL_H

#include <cmath>

#include "core/terms/QuasiperiodicDisorder.h"
#include "core/terms/LookupCavityZ2.h"
#include "core/terms/LookupCavityYZ.h"
#include "core/terms/LookupCavityY2.h"
#include "core/terms/CavityLongInteraction.h"
#include "core/terms/QuasiperiodicDisorder.h"
#include "core/terms/OnsiteDisorder.h"

#include "utils/Assertions.h"
#include "HamiltonianGenerator.h"
#include "RND.h"

/**
 * @brief A class responsible for a specific reparametrising the HamiltonianGenerator for each new simulation.
 * @details Implementations can do stuff like sampling new random disorder terms, introducing random (or uniformly
 * spaces) phases, and so on.
 */
class AveragingModel {
public:
    /**
     * @brief Prepares @a hamiltonianGenerator for a new simulation.
     * @details @a rnd is used for resampling random stuff, while @a simulationIndex and @a numberOfSimulations
     * probably form uniformly spreading some parameter among all simulations.
     */
    virtual void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                           std::size_t simulationIndex, std::size_t numberOfSimulations) = 0;
};

#endif //MBL_ED_AVERAGINGMODEL_H
