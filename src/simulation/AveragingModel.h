//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_AVERAGINGMODEL_H
#define MBL_ED_AVERAGINGMODEL_H

#include <cmath>
#include <simulation/terms/QuasiperiodicDisorder.h>

#include "simulation/terms/LookupCavityZ2.h"
#include "simulation/terms/LookupCavityYZ.h"
#include "simulation/terms/LookupCavityY2.h"
#include "simulation/terms/CavityLongInteraction.h"
#include "simulation/terms/QuasiperiodicDisorder.h"
#include "simulation/terms/OnsiteDisorder.h"

#include "utils/Assertions.h"
#include "HamiltonianGenerator.h"
#include "RND.h"

class AveragingModel {
public:
    virtual void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                      std::size_t simulationIndex, std::size_t numberOfSimulations) = 0;
};

#endif //MBL_ED_AVERAGINGMODEL_H
