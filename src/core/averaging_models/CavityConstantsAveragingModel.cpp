//
// Created by pkua on 08.06.2020.
//

#include "CavityConstantsAveragingModel.h"
#include "core/HamiltonianGenerator.h"
#include "core/RND.h"
#include "core/terms/LookupCavityZ2.h"
#include "core/terms/OnsiteDisorder.h"
#include "core/terms/LookupCavityYZ.h"
#include "core/terms/LookupCavityY2.h"

void CavityConstantsAveragingModel::setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                                              unsigned long simulationIndex,
                                                              unsigned long numberOfSimulations)
{
    Expects(numberOfSimulations > 0);
    Expects(simulationIndex < numberOfSimulations);

    bool termFound{};
    for (auto &diagonalTerm : hamiltonianGenerator.getDiagonalTerms()) {
        try {
            auto &lookupCavityZ2 = dynamic_cast<LookupCavityZ2&>(*diagonalTerm);
            lookupCavityZ2.changeRealisation(simulationIndex);
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }

        try {
            auto &onsiteDisorder = dynamic_cast<OnsiteDisorder&>(*diagonalTerm);
            onsiteDisorder.resampleOnsiteEnergies(rnd);
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }
    }

    for (auto &hoppingTerm : hamiltonianGenerator.getHoppingTerms()) {
        try {
            auto &lookupCavityYZ = dynamic_cast<LookupCavityYZ&>(*hoppingTerm);
            lookupCavityYZ.changeRealisation(simulationIndex);
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }
    }

    for (auto &doubleHoppingTerm : hamiltonianGenerator.getDoubleHoppingTerms()) {
        try {
            auto &lookupCavityY2 = dynamic_cast<LookupCavityY2&>(*doubleHoppingTerm);
            lookupCavityY2.changeRealisation(simulationIndex);
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }
    }

    if (!termFound)
        throw std::runtime_error("CavityConstantsAveragingModel: lack of any term to average on");
}