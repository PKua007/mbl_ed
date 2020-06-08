//
// Created by pkua on 08.06.2020.
//

#include "OnsiteDisorderAveragingModel.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/RND.h"
#include "simulation/terms/OnsiteDisorder.h"

void OnsiteDisorderAveragingModel::setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                                             unsigned long simulationIndex,
                                                             unsigned long numberOfSimulations) {
    static_cast<void>(simulationIndex);
    static_cast<void>(numberOfSimulations);

    bool changed{};
    for (auto &diagonalTerm : hamiltonianGenerator.getDiagonalTerms()) {
        try {
            auto &onsiteDisorder = dynamic_cast<OnsiteDisorder&>(*diagonalTerm);
            onsiteDisorder.resampleOnsiteEnergies(rnd);
            changed = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }
    }

    if (!changed)
        throw std::runtime_error("OnsiteDisorderAveragingModel: lack of any term to average on");
}