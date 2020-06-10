//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_UNIFORMPHI0AVERAGINGMODEL_H
#define MBL_ED_UNIFORMPHI0AVERAGINGMODEL_H


#include "simulation/AveragingModel.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/RND.h"

/**
 * @brief Averaging model, used in Simulation, which distributes evenly phi0 values across simulations based on
 * @a simulationIndex and @a numberOfSimulations.
 * @details Onsite disorder is also resampled.
 */
class UniformPhi0AveragingModel : public AveragingModel {
public:
    /**
     * @brief The method looks for CavityLongInteraction and OnsiteDisorder terms in the @a hamiltonianGenerator.
     */
    void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd, std::size_t simulationIndex,
                                   std::size_t numberOfSimulations) override;
};


#endif //MBL_ED_UNIFORMPHI0AVERAGINGMODEL_H
