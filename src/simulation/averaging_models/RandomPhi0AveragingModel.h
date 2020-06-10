//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_RANDOMPHI0AVERAGINGMODEL_H
#define MBL_ED_RANDOMPHI0AVERAGINGMODEL_H


#include "simulation/AveragingModel.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/RND.h"

/**
 * @brief Averaging model, used in Simulation, which picks random phi0 values.
 * @details Onsite disorder is also resampled.
 */
class RandomPhi0AveragingModel : public AveragingModel {
public:
    /**
     * @brief The method looks for CavityLongInteraction and OnsiteDisorder terms in the @a hamiltonianGenerator.
     */
    void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd, std::size_t simulationIndex,
                                   std::size_t numberOfSimulations) override;
};


#endif //MBL_ED_RANDOMPHI0AVERAGINGMODEL_H
