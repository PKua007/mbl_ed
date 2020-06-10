//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_CAVITYCONSTANTSAVERAGINGMODEL_H
#define MBL_ED_CAVITYCONSTANTSAVERAGINGMODEL_H

#include "simulation/AveragingModel.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/RND.h"

/**
 * @brief Averaging model, used in Simulation which picks subsequent realisations from CavityConstants.
 * @details It also resamples random onsite disorder.
 * @tparam DisorderGenerator_t disorder generator used
 */
class CavityConstantsAveragingModel : public AveragingModel {
public:
    /**
     * @brief The method looks for LookupCavityZ2, LookupCavityYZ and OnsiteDisorder terms in the
     * @a hamiltonianGenerator.
     */
    void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd, std::size_t simulationIndex,
                                   std::size_t numberOfSimulations) override;
};

#endif //MBL_ED_CAVITYCONSTANTSAVERAGINGMODEL_H
