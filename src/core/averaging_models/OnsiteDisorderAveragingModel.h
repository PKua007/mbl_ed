//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_ONSITEDISORDERAVERAGINGMODEL_H
#define MBL_ED_ONSITEDISORDERAVERAGINGMODEL_H

#include "core/AveragingModel.h"
#include "core/HamiltonianGenerator.h"
#include "core/RND.h"

/**
 * @brief Averaging model, used in Simulation, which samples new onsite disorder terms for each simulation.
 */
class OnsiteDisorderAveragingModel : public AveragingModel {
public:
    /**
     * @brief The method looks for OnsiteDisorder term in the @a hamiltonianGenerator.
     */
    void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd, std::size_t simulationIndex,
                                   std::size_t numberOfSimulations) override;
};

#endif //MBL_ED_ONSITEDISORDERAVERAGINGMODEL_H
