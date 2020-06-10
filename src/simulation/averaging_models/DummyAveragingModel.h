//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_DUMMYAVERAGINGMODEL_H
#define MBL_ED_DUMMYAVERAGINGMODEL_H

#include "simulation/AveragingModel.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/RND.h"

/**
 * @brief Averaging model, which does nothing.
 */
class DummyAveragingModel : public AveragingModel {
public:
    void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd, std::size_t simulationIndex,
                                   std::size_t numberOfSimulations) override
    {
        static_cast<void>(simulationIndex);
        static_cast<void>(numberOfSimulations);
        static_cast<void>(rnd);
        static_cast<void>(hamiltonianGenerator);
    }
};


#endif //MBL_ED_DUMMYAVERAGINGMODEL_H
